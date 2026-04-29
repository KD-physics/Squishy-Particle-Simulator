"""
test_D0_single_drop_repro.py — Reproduce paper Benchmark E (falling droplet)
using the new System / ParticleSpec API.

Paper parameters (emulsion_paper_figures.py):
  N=120, gamma=1.0, K_area=50 (κ=0.02), C=500, rho_d=1.0, k_reg=10.0
  alpha_damp=5.0, dt = 0.1*R0/(c_cap*N) ≈ 8.33e-4
  Bo=0.10, H_drop=5.0 R0, t_max=80 τ0

Expected results from paper:
  - Trajectory error vs free-fall < 3e-4 R0
  - Sag ratio (R_top - R_bot)/R0 ≈ 0.065 at Bo=0.10
  - Area conservation |ΔA/A| < 0.3%
  - No perimeter collapse

Node sweep: N = 120 → 60 → 36 to find stability floor.

Usage:
    python src/validation/test_D0_single_drop_repro.py
"""

import sys, os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Wall
from src.epd.system import System

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--nodes', type=int, nargs='+', default=[120, 60, 36],
                    help='Node counts to test (default: 120 60 36)')
_args = parser.parse_args()

os.makedirs('results', exist_ok=True)

print("=" * 60)
print("Test D0: single falling droplet — paper Benchmark E repro")
print("=" * 60)

# ── Paper parameters (emulsion_paper_figures.py) ──────────────────────────────
GAMMA     = 1.0
KAPPA     = 0.02     # area compressibility; C̃=500 default now correct in ParticleSpec
ALPHA     = 5.0
RHO_D     = 1.0
BO        = 0.10     # paper Benchmark E
G_SIM     = BO       # with gamma=rho_d=R0=1, Bo = g

H_DROP    = 5.0      # droplet CM starts H_DROP R0 above floor contact point
T_MAX_TAU = 80.0     # run 80 τ0 (τ0 = sqrt(rho_d R0^3/gamma) = 1)

# ── Node sweep ────────────────────────────────────────────────────────────────
NODE_COUNTS = _args.nodes

results = {}

for N in NODE_COUNTS:
    print(f"\n{'─'*50}")
    print(f"N = {N} nodes")
    print(f"{'─'*50}")

    # dt follows paper formula: 0.1 * R0 / (c_cap * N), c_cap=1
    # We pass it explicitly so System doesn't override it
    dt_paper = 0.1 / N
    tau0     = 1.0   # sqrt(rho_d * R0^3 / gamma) = 1
    t_max    = T_MAX_TAU * tau0
    n_steps  = int(t_max / dt_paper) + 1

    print(f"  dt = {dt_paper:.2e}   n_steps = {n_steps:,}   t_max = {t_max:.1f} τ0")

    # ── Build system ──────────────────────────────────────────────────────────
    # Box: wide enough that the floor wall covers the droplet
    LX = 20.0
    # Droplet starts at y = r_c + H_DROP ≈ 6 above floor at y=0
    # Use LY large enough to contain it
    y_cm0 = H_DROP + 1.5   # approximate; r_c ≈ 1.5 × L0/2
    LY    = y_cm0 + 5.0

    sys0 = System(LX, LY, periodic_x=False, periodic_y=False,
                  g=G_SIM)

    # Floor only (no side walls — single droplet, no lateral confinement)
    floor = Wall((0.0, 0.0), (LX, 0.0), normal=(0.0, 1.0))
    floor.set_render(color='#333333', linewidth=3.0, alpha=0.9)
    sys0.add_object(floor)

    spec = ParticleSpec(count=1, type='emulsion',
                        gamma=GAMMA, kappa=KAPPA,
                        alpha_damp=ALPHA,
                        N_nodes=N)
    sys0.add_particles(spec)

    d = spec.derived
    print(f"  γ={d['gamma']:.2f}  κ={d['kappa']:.3f}  K_area={d['K_area']:.1f}"
          f"  C={d['C']:.0f}  (C̃={d['C']/GAMMA:.0f})  α={d['alpha']:.1f}")

    # Place droplet at y_cm0 above floor, centred in x
    x_center = LX / 2.0

    # Initialize and then manually position
    sys0.initialize(phi_target=0.80, seed=0, verbose=False,
                    relax_only=True, n_relax_init=0)

    # Override dt to match paper formula exactly
    import src.simulation.tf_sim as tf_sim
    NP_DTYPE = np.float64
    sys0._dt     = dt_paper
    sys0._dt_tf  = tf.constant(NP_DTYPE(dt_paper))

    # Reposition droplet CM to (x_center, y_cm0)
    state = sys0.state
    x_cm_np = state['x_cm'].numpy()
    dx = x_center - float(x_cm_np[0, 0])
    dy = y_cm0    - float(x_cm_np[0, 1])
    import src.simulation.tf_sim as tfsim

    # Shift all node positions and CM
    x_all_np  = state['x_all'].numpy()
    x_cm_np[0] += [dx, dy]
    x_all_np[0] += [dx, dy]

    DTYPE = tf.float64
    new_state = dict(state)
    new_state['x_all'] = tf.constant(x_all_np, dtype=DTYPE)
    new_state['x_cm']  = tf.constant(x_cm_np,  dtype=DTYPE)
    sys0._state = new_state

    actual_y_cm0 = float(sys0.state['x_cm'].numpy()[0, 1])
    print(f"  Droplet CM at y = {actual_y_cm0:.3f}   (target {y_cm0:.3f})")

    # ── Run ───────────────────────────────────────────────────────────────────
    sample_every = max(1, n_steps // 80)   # ~80 frames
    print(f"  Running {n_steps:,} steps (sample every {sample_every})…")
    sys0.run(n_steps, sample_every=sample_every, verbose=False)

    # ── Metrics ───────────────────────────────────────────────────────────────
    frames = sys0.frames

    # Trajectory error vs free-fall (before floor contact)
    # Approximate impact: y_cm0 - 0.5*g*t^2 = r_c → t_impact ~ sqrt(2*(y_cm0-1)/g)
    t_impact = np.sqrt(2.0 * max(0.1, actual_y_cm0 - 1.0) / G_SIM)
    ts   = np.array([fr['t'] for fr in frames])
    ycms = np.array([fr['x_cm'][0, 1] for fr in frames])
    mask = (ts > 0) & (ts < 0.9 * t_impact)
    traj_err = 0.0
    if mask.sum() > 3:
        y_an = actual_y_cm0 - 0.5 * G_SIM * ts[mask]**2
        traj_err = float(np.max(np.abs(ycms[mask] - y_an)) / H_DROP)
    print(f"  Trajectory error = {traj_err:.2e} R0  (paper: <3e-4)")

    # Final shape metrics (last frame)
    x_final = frames[-1]['x_all'][0]   # (N, 2)
    x_cm_f  = frames[-1]['x_cm'][0]    # (2,)
    R_top = float(x_final[:, 1].max() - x_cm_f[1])
    R_bot = float(x_cm_f[1] - x_final[:, 1].min())
    sag   = (R_top - R_bot)   # /R0 with R0=1
    print(f"  Sag (R_top - R_bot)/R0 = {sag:.4f}  (paper Bo=0.10: ≈ 0.065)")

    # Area conservation (last frame)
    # Shoelace area from final x_all
    x_fin = x_final
    x_next = np.roll(x_fin, -1, axis=0)
    A_final = 0.5 * abs(np.sum(x_fin[:,0]*x_next[:,1] - x_next[:,0]*x_fin[:,1]))
    A0_circ = np.pi   # pi R0^2 with R0=1
    dA = abs(A_final / A0_circ - 1.0)
    print(f"  Area error |ΔA/A| = {dA:.4f}  (paper: <0.3%)")

    # Minimum circularity across all frames (collapse indicator)
    circ_min = 1.0
    for fr in frames:
        xn = fr['x_all'][0]
        xnext = np.roll(xn, -1, axis=0)
        edges = np.linalg.norm(xnext - xn, axis=1)
        L_perim = edges.sum()
        A_fr = 0.5 * abs(np.sum(xn[:,0]*xnext[:,1] - xnext[:,0]*xn[:,1]))
        circ = 4.0 * np.pi * A_fr / (L_perim**2)
        circ_min = min(circ_min, circ)
    print(f"  Min circularity = {circ_min:.4f}  (collapse if << 1)")

    results[N] = dict(traj_err=traj_err, sag=sag, dA=dA, circ_min=circ_min)

    # ── Movie ─────────────────────────────────────────────────────────────────
    out = f'results/test_D0_N{N}.gif'
    sys0.make_movie(
        out, fps=10,
        xlim=(x_center - 3, x_center + 3),
        ylim=(-1.5, actual_y_cm0 + 2.0),
        title=f'Single drop  N={N}  Bo={BO}  κ={KAPPA}  C̃=500',
    )
    print(f"  Movie → {out}")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Summary (Bo={BO}, κ={KAPPA}, C̃=500)")
print(f"{'='*60}")
print(f"{'N':>6}  {'traj_err':>10}  {'sag':>8}  {'dA':>8}  {'circ_min':>10}")
print(f"{'paper':>6}  {'<3e-4':>10}  {'≈0.065':>8}  {'<0.3%':>8}  {'≈1.00':>10}")
print(f"{'─'*60}")
for N, r in results.items():
    print(f"{N:>6}  {r['traj_err']:>10.2e}  {r['sag']:>8.4f}"
          f"  {r['dA']:>8.4f}  {r['circ_min']:>10.4f}")
