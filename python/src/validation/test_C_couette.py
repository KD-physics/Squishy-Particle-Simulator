"""
test_C_couette.py — Test C Part 1: Couette cell seeding + swell

Geometry:  annular cell, R_inner:R_outer = 1:5
           50 elastic particles, N_nodes=36, nu=0.5
           phi_target=0.80 (fraction of annulus area)
           periodic_x=False, periodic_y=False

PASS criteria:
  - phi_outer (annulus-corrected) >= 0.78 after swell
  - All particle CMs inside annulus at every recorded frame
  - GIF movie saved

Usage:
    python src/validation/test_C_couette.py [--out PATH]
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import CouetteCell
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=str, default='results/test_C_couette.gif')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test C (Part 1): Couette cell — seeding + swell")
print("=" * 60)

# ── Geometry ──────────────────────────────────────────────────────────────────
# R_inner : R_outer = 1 : 5
# Size chosen so 50 particles (R0≈1) fit comfortably at RSA phi~0.25
R_inner = 3.0
R_outer = 15.0      # ratio 1:5
Lx = Ly = 2.0 * R_outer   # box exactly wraps the outer circle

print(f"\nGeometry: R_inner={R_inner}  R_outer={R_outer}  box={Lx}×{Ly}")
print(f"  annulus area = {np.pi*(R_outer**2-R_inner**2):.2f}  "
      f"(box = {Lx*Ly:.2f}  ratio = {np.pi*(R_outer**2-R_inner**2)/(Lx*Ly):.4f})")

# ── System ────────────────────────────────────────────────────────────────────
sys_c = System(Lx, Ly, periodic_x=False, periodic_y=False)

cell = CouetteCell(inner_radius=R_inner, outer_radius=R_outer,
                   x0=R_outer, y0=R_outer)   # centred in box
cell.set_render(color='#333333', linewidth=2.0, alpha=0.9)

sys_c.add_object(cell)

spec = ParticleSpec(count=50, nu=0.5, N_nodes=36, poly_dist=0.05)
sys_c.add_particles(spec)

# ── Initialize ────────────────────────────────────────────────────────────────
print()
sys_c.initialize(phi_target=0.80, seed=42, verbose=True,
                 n_relax=200, swell_alpha=10.0)

phi_final = sys_c.phi_outer
ckpt_path = 'results/couette_phi08.npz'
sys_c.save_state(ckpt_path)
print(f"\nCheckpoint saved → {ckpt_path}")
print(f"Post-swell (annulus phi): {phi_final:.4f}")
print(f"Post-swell box: Lx={sys_c.Lx:.3f}  Ly={sys_c.Ly:.3f}")
r_outer_final = sys_c.Lx / 2.0
r_inner_final = r_outer_final / 5.0
print(f"Post-swell cell: R_inner≈{r_inner_final:.3f}  R_outer≈{r_outer_final:.3f}")
print(f"  mean R0 = {np.mean([p.R0 for p in sys_c._particles]):.4f}")

# ── Short run ─────────────────────────────────────────────────────────────────
n_steps      = 500
sample_every = max(1, n_steps // 40)
print(f"\nRunning {n_steps} steps (sample every {sample_every})…")
sys_c.run(n_steps, sample_every=sample_every, verbose=True)

# ── Physics check: all CMs in annulus ─────────────────────────────────────────
n_violations = 0
cx = sys_c.Lx / 2.0
cy = sys_c.Ly / 2.0
# After swell, R_outer and R_inner scale with box
r_out = sys_c.Lx / 2.0
r_in  = r_out / 5.0
for fr in sys_c.frames:
    for cm in fr['x_cm']:
        r = np.sqrt((cm[0] - cx)**2 + (cm[1] - cy)**2)
        if r <= r_in or r >= r_out:
            n_violations += 1

phi_ok = phi_final >= 0.78
excl_ok = n_violations == 0
print(f"\n── Physics check ──")
print(f"  phi_outer (annulus) = {phi_final:.4f}   {'PASS (≥0.78)' if phi_ok else 'FAIL'}")
print(f"  CM outside annulus  = {n_violations}   {'PASS' if excl_ok else 'FAIL'}")

# ── Movie ──────────────────────────────────────────────────────────────────────
print("\nRendering movie…")
sys_c.make_movie(
    args.out, fps=10,
    xlim=(-0.5, sys_c.Lx + 0.5),
    ylim=(-0.5, sys_c.Ly + 0.5),
    title=f'Couette cell (R_in:R_out=1:5)  ν=0.5  φ={phi_final:.3f}',
)
print(f"Movie saved → {args.out}")

all_pass = phi_ok and excl_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test C Part 1 {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
