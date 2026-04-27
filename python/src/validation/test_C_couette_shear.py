"""
test_C_couette_shear.py — Test C Part 2: Couette cell shear with rotating inner wall

Geometry: same annular cell as Part 1 (R_inner:R_outer = 1:5)
          50 elastic particles, N_nodes=36, nu=0.5
          Swell to phi=0.80, then set inner-wall particles as driven + frozen.

Driven particles: all with CM within 0.8 * 2*R0_mean of the inner wall.
                  They orbit about cell center at omega_orbit=0.5 (frozen shape).

PASS criteria:
  - Post-swell phi_outer >= 0.78
  - At least 1 driven particle (inner-wall roughness elements)
  - All particle CMs remain inside annulus at every frame
  - GIF movie saved

Usage:
    python src/validation/test_C_couette_shear.py [--out PATH] [--ckpt PATH]
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.simulation.tf_sim import make_traj
from src.epd.particles import ParticleSpec
from src.epd.objects import CouetteCell
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out',  type=str, default='results/test_C_couette_shear.gif')
parser.add_argument('--ckpt', type=str, default='results/couette_phi08.npz')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test C (Part 2): Couette cell shear — rotating inner wall")
print("=" * 60)

# ── Geometry (identical to Part 1) ────────────────────────────────────────────
R_inner = 3.0
R_outer = 15.0
Lx = Ly = 2.0 * R_outer

print(f"\nGeometry: R_inner={R_inner}  R_outer={R_outer}  box={Lx}×{Ly}")

# ── System + spec (identical to Part 1) ───────────────────────────────────────
spec = ParticleSpec(count=50, nu=0.5, N_nodes=36, poly_dist=0.05)


def _make_system():
    s = System(Lx, Ly, periodic_x=False, periodic_y=False)
    cell = CouetteCell(inner_radius=R_inner, outer_radius=R_outer,
                       x0=R_outer, y0=R_outer)
    cell.set_render(color='#333333', linewidth=2.0, alpha=0.9)
    s.add_object(cell)
    s.add_particles(spec)
    return s, cell


# ── Initialize: load from checkpoint or re-swell ──────────────────────────────
if os.path.exists(args.ckpt):
    print(f"\nFast path: restoring post-swell state from {args.ckpt}")
    sys_c, cell = _make_system()
    # Quick initialize to build particles/params — positions will be overwritten
    sys_c.initialize(phi_target=0.80, seed=42, verbose=True,
                     relax_only=True, n_relax_init=0)
    sys_c.restore_state(args.ckpt)
    print(f"  Restored: Lx={sys_c.Lx:.3f}  phi_outer={sys_c.phi_outer:.4f}")
else:
    print(f"\nSlow path: no checkpoint found — running full swell (~30 min)")
    sys_c, cell = _make_system()
    sys_c.initialize(phi_target=0.80, seed=42, verbose=True,
                     n_relax=200, swell_alpha=10.0)
    sys_c.save_state(args.ckpt)
    print(f"  Checkpoint saved → {args.ckpt}")

phi_swell = sys_c.phi_outer
print(f"\nPost-swell phi_outer = {phi_swell:.4f}")
print(f"Post-swell box: Lx={sys_c.Lx:.3f}  Ly={sys_c.Ly:.3f}")

# ── Identify inner-wall particles ─────────────────────────────────────────────
# After swell, R_outer_final = Lx/2 (the box exactly fits the outer wall)
R0_mean       = float(np.mean([p.R0 for p in sys_c.particles]))
R_outer_final = sys_c.Lx / 2.0
R_inner_final = R_outer_final / 5.0
cx, cy        = R_outer_final, R_outer_final   # cell center in box coords

x_cm = sys_c.state['x_cm'].numpy()    # (P, 2)
r_cm = np.sqrt((x_cm[:, 0] - cx)**2 + (x_cm[:, 1] - cy)**2)

# Threshold: CM closer than 0.8 * mean_diameter from inner wall surface
thresh = R_inner_final + 0.8 * 2.0 * R0_mean
inner_idxs = [i for i in range(len(sys_c.particles)) if r_cm[i] < thresh]

print(f"\nInner-wall selection:")
print(f"  R_inner_final = {R_inner_final:.3f}  R0_mean = {R0_mean:.4f}")
print(f"  threshold r_cm < {thresh:.3f}  (= R_inner + 0.8×2×R0)")
print(f"  driven particles: {len(inner_idxs)}  indices={inner_idxs}")

assert len(inner_idxs) >= 1, "No driven particles found — check geometry or threshold"

# ── Assign orbital trajectories (frozen shape) ────────────────────────────────
Omega = 0.5   # rad / time-unit; CCW rotation of inner wall

traj_rows = [
    make_traj(omega_orbit_dc=Omega, r_ref=(cx, cy))
    for _ in inner_idxs
]
sys_c.set_driven_particles(inner_idxs, traj_rows, frozen=True)

# Color driven particles distinctly
for i in inner_idxs:
    sys_c._particle_colors[i] = '#e05252'   # red tint for driven

print(f"\nDriven particles assigned omega_orbit={Omega} about ({cx:.3f},{cy:.3f})")

# ── Run shear ─────────────────────────────────────────────────────────────────
n_steps      = 500
sample_every = max(1, n_steps // 40)
print(f"\nRunning {n_steps} steps (sample every {sample_every})…")
sys_c.run(n_steps, sample_every=sample_every, verbose=True)

# ── Physics check ─────────────────────────────────────────────────────────────
n_violations = 0
r_out_f = sys_c.Lx / 2.0
r_in_f  = r_out_f / 5.0
cx_f    = sys_c.Lx / 2.0
cy_f    = sys_c.Ly / 2.0
for fr in sys_c.frames:
    for cm in fr['x_cm']:
        r = np.sqrt((cm[0] - cx_f)**2 + (cm[1] - cy_f)**2)
        if r <= r_in_f or r >= r_out_f:
            n_violations += 1

phi_ok  = phi_swell >= 0.78
driven_ok = len(inner_idxs) >= 1
excl_ok = n_violations == 0
print(f"\n── Physics check ──")
print(f"  phi_outer          = {phi_swell:.4f}   {'PASS' if phi_ok else 'FAIL'} (≥0.78)")
print(f"  driven particles   = {len(inner_idxs)}   {'PASS' if driven_ok else 'FAIL'} (≥1)")
print(f"  CM outside annulus = {n_violations}   {'PASS' if excl_ok else 'FAIL'}")

# ── Movie ─────────────────────────────────────────────────────────────────────
print("\nRendering movie…")
sys_c.make_movie(
    args.out, fps=10,
    xlim=(-0.5, sys_c.Lx + 0.5),
    ylim=(-0.5, sys_c.Ly + 0.5),
    title=(f'Couette shear  Ω_inner={Omega}  ν=0.5  φ={phi_swell:.3f}  '
           f'N_driven={len(inner_idxs)}'),
)
print(f"Movie saved → {args.out}")

all_pass = phi_ok and driven_ok and excl_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test C Part 2 {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
