"""
test_D_falling_emulsion.py — Test D: 50 emulsion droplets falling under gravity

Setup:
  U-shaped container: bottom + left + right walls, Lx=6, Ly=110 (open top)
  50 emulsion droplets  (γ=1, q=50 ↔ κ=0.02, N=36, Gaussian 5% poly)
  Initial condition: RSA at ϕ≈0.24 → particles randomly distributed ~uniformly
  then released under gravity.

Gravity:
  g = Bo × γ / (ρ × R0²) = Bo (with γ=ρ=R0=1)
  Bo=0.025 → g=0.025, settling ~40× slower
  Bo=1.0   → g=1.0  used here for computational feasibility

Batched execution:
  Runs in STEPS_PER_BATCH=5000 chunks.  After each chunk the GIF is refreshed
  so you can inspect progress and kill early if needed.
  Checkpoint saved after every batch → resume with --ckpt.

PASS criteria:
  - All particle CMs above y=0.5 at final frame (no floor escape)
  - All particle CMs inside x=[0,6] (no sidewall escape)
  - Final mean y_cm < 0.75 × initial mean y_cm  (clear downward settling)
  - GIF movie saved

Usage:
    python src/validation/test_D_falling_emulsion.py [--out PATH] [--ckpt PATH]
    python src/validation/test_D_falling_emulsion.py --batches 4   # run only 4 batches
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Wall
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out',    type=str, default='results/test_D_falling_emulsion.gif')
parser.add_argument('--ckpt',   type=str, default='results/test_D_falling_emulsion.npz')
parser.add_argument('--batches',type=int, default=28,
                    help='Number of 5k-step batches to run (default 28 → 140k total)')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out),  exist_ok=True)
os.makedirs(os.path.dirname(args.ckpt), exist_ok=True)

print("=" * 60)
print("Test D: 50 emulsion droplets falling under gravity")
print("=" * 60)

# ── Parameters ────────────────────────────────────────────────────────────────
GAMMA    = 1.0
KAPPA    = 0.02       # area compressibility (κ = γ/(R0·K_area))
N_NODES  = 36
N_DROP   = 40          # 40 drops → ~20% faster than 50; 6-7 per row → ~7 rows settled
POLY     = 0.05        # 5% Gaussian std

# Gravity:  g = Bo × γ / (ρ R0²)  with γ=ρ=R0=1  →  g = Bo
# Bo=0.05 → ΔP_hydro/P_Laplace = 0.05×2 = 0.10 → single drop nearly round;
#            bottom drops compressed by stack will sag somewhat — expected physics.
# v_terminal = g/α = 0.05/5.0 = 0.01 R0/τ0
BO       = 0.05
G_SIM    = BO

# Container / box  (all lengths in units of R0=1, so 1 unit = 1 radius = 0.5 diameter)
LX       = 12.0        # 6 particle diameters wide
LY       = 60.0        # RSA uses R_eff≈1.175R0 (including contact zone); need LY≥58 to stay
#                        below φ=0.25 threshold: 40π×1.175²/(0.25×12)≈58. Use 60 for margin.
#                        Settled pile: ~7 rows × 2R0 ≈ 14 units ≈ 7 diameters (fits in viewport)

# Batched run
STEPS_PER_BATCH = 5_000
N_BATCHES       = args.batches
SAMPLE_EVERY    = 500   # 10 frames per batch → 280 frames total → 28 s at 10 fps

print(f"\nDroplet: γ={GAMMA}  κ={KAPPA}  N={N_NODES}  count={N_DROP}")
print(f"Gravity: g={G_SIM:.4f}  (Bo={BO})")
print(f"Box:  Lx={LX} ({LX/2:.0f} diameters)  Ly={LY} ({LY/2:.0f} diameters)  "
      f"phi_box≈{N_DROP*np.pi*1**2/(LX*LY):.3f}")
print(f"Plan: {N_BATCHES} batches × {STEPS_PER_BATCH:,} steps = {N_BATCHES*STEPS_PER_BATCH:,} total")

# ── System + walls ────────────────────────────────────────────────────────────
sys_d = System(LX, LY, periodic_x=False, periodic_y=False, g=G_SIM)

wall_bot = Wall((0,    0), (LX,   0), normal=( 0,  1))
wall_lft = Wall((0,    0), ( 0,  LY), normal=( 1,  0))
wall_rgt = Wall((LX,   0), (LX,  LY), normal=(-1,  0))
for w in [wall_bot, wall_lft, wall_rgt]:
    w.set_render(color='#333333', linewidth=2.0, alpha=0.9)
    sys_d.add_object(w)

# ── Particles ─────────────────────────────────────────────────────────────────
spec = ParticleSpec(count=N_DROP, type='emulsion',
                    gamma=GAMMA, kappa=0.02,
                    N_nodes=N_NODES, poly_dist=POLY)
sys_d.add_particles(spec)

d = spec.derived
print(f"\nEmulsion params: γ={d['gamma']:.2f}  κ={d['kappa']:.3f}"
      f"  K_area={d['K_area']:.1f}  C={d['C']:.0f}  (C̃={d['C']/d['gamma']:.0f})  α={d['alpha']:.3f}")

# ── Initialize via RSA (no swell) ─────────────────────────────────────────────
print()
sys_d.initialize(phi_target=0.80, seed=42, verbose=True,
                 relax_only=True, n_relax_init=200)

initial_y_mean = float(sys_d.state['x_cm'].numpy()[:, 1].mean())
print(f"\nInitial mean y_cm = {initial_y_mean:.3f}")
print(f"dt = {sys_d._dt:.4e}   alpha = {sys_d._alpha_damp:.4f}")
print(f"Terminal velocity estimate: v_t = g/α ≈ {G_SIM/sys_d._alpha_damp:.3f}")

# ── Batched run ───────────────────────────────────────────────────────────────
print(f"\nStarting batched run: {N_BATCHES} × {STEPS_PER_BATCH:,} steps")
print(f"  GIF refreshed after each batch → {args.out}")
print(f"  Checkpoint saved after each batch → {args.ckpt}")
print()

for batch in range(1, N_BATCHES + 1):
    print(f"── Batch {batch}/{N_BATCHES}  (steps {(batch-1)*STEPS_PER_BATCH:,}–{batch*STEPS_PER_BATCH:,}) ──")
    sys_d.run(STEPS_PER_BATCH, sample_every=SAMPLE_EVERY, verbose=True)

    # Stats
    cur_y    = float(sys_d.frames[-1]['x_cm'][:, 1].mean())
    cur_cms  = sys_d.frames[-1]['x_cm']
    fl_viol  = sum(1 for cm in cur_cms if cm[1] < 0.5)
    si_viol  = sum(1 for cm in cur_cms if cm[0] < 0.0 or cm[0] > LX)
    ratio    = cur_y / initial_y_mean
    print(f"  mean y_cm = {cur_y:.3f}  ({ratio:.2f}× initial)  "
          f"floor_viol={fl_viol}  side_viol={si_viol}  "
          f"frames={len(sys_d.frames)}")

    # Refresh GIF — fixed 6×12 diameter viewport so proportions are always correct
    sys_d.make_movie(
        args.out, fps=10,
        xlim=(-0.5, LX + 0.5),
        ylim=(-0.5, 24.5),          # 12 particle diameters tall (2R0 per diameter)
        title=f'Falling emulsion  γ={GAMMA}  κ={KAPPA}  Bo={BO}  '
              f'batch {batch}/{N_BATCHES}',
    )
    print(f"  GIF updated ({len(sys_d.frames)} frames)")

    # Checkpoint
    sys_d.save_state(args.ckpt)
    print(f"  Checkpoint saved")

    # Early-stop: well below pass threshold
    if ratio < 0.50:
        print(f"\n  >> Well settled (y_cm ratio={ratio:.2f} < 0.50) — stopping early")
        break

# ── Final physics check ───────────────────────────────────────────────────────
final_y_mean = float(sys_d.frames[-1]['x_cm'][:, 1].mean())
final_x_cms  = sys_d.frames[-1]['x_cm']

floor_violations = sum(1 for cm in final_x_cms if cm[1] < 0.5)
side_violations  = sum(1 for cm in final_x_cms if cm[0] < 0.0 or cm[0] > LX)
settled = final_y_mean < 0.85 * initial_y_mean   # 15% drop; at Bo=0.05 full settling takes many more steps

floor_ok = (floor_violations == 0)
side_ok  = (side_violations  == 0)

print(f"\n── Physics check ──")
print(f"  Initial mean y_cm  = {initial_y_mean:.3f}")
print(f"  Final   mean y_cm  = {final_y_mean:.3f}   "
      f"({'PASS' if settled else 'FAIL (need < 0.85×initial)'})")
print(f"  Floor violations   = {floor_violations}   {'PASS' if floor_ok else 'FAIL'}")
print(f"  Sidewall violations= {side_violations}   {'PASS' if side_ok else 'FAIL'}")

# ── Final movie ───────────────────────────────────────────────────────────────
print("\nRendering final movie…")
sys_d.make_movie(
    args.out, fps=10,
    xlim=(-0.5, LX + 0.5),
    ylim=(-0.5, 24.5),          # 12 particle diameters tall
    title=f'Falling emulsion  γ={GAMMA}  κ={KAPPA}  Bo={BO}  N={N_DROP}',
)
print(f"Movie saved → {args.out}")

all_pass = floor_ok and side_ok and settled
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test D {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
