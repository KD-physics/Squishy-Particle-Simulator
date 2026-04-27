"""
test_B1_seed.py — Test B1: Channel + SquareObstacle seeding

Verifies:
  1. Channel (periodic x, hard y) + SquareObstacle (center, exclusion='interior')
     are correctly assembled and passed to System.
  2. RSA seeder places 12 particles (N=32, nu=0.5) with NONE inside the box.
  3. phi_outer correctly subtracts obstacle area from denominator.
  4. A snapshot PNG is saved for visual inspection.

PASS criteria:
  - All particle CMs outside obstacle region
  - phi_box_raw    = sum(R0²*pi) / (Lx*Ly)         (naive, ignores obstacle)
  - phi_outer      = system.phi_outer                (corrected, subtracts obstacle)
  - phi_outer > phi_box_raw * Lx*Ly / accessible_area  (corrected > naive by ~box_area/total)

Usage:
    python src/validation/test_B1_seed.py [--out PATH]
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Channel, SquareObstacle
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=str, default='results/test_B1_seed.png')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test B1: Channel + SquareObstacle seeding")
print("=" * 60)

# ── Geometry ──────────────────────────────────────────────────────────────────
Lx, Ly = 20.0, 20.0
box_side = Lx / 4.0          # 5.0 = 1/4 of Lx
box_cx   = Lx / 2.0          # centered at (10, 10)
box_cy   = Ly / 2.0

# ── Build system ──────────────────────────────────────────────────────────────
sys_b1 = System(Lx, Ly, periodic_x=True)

channel = Channel(width=Lx, height=Ly, x0=Lx/2, y0=Ly/2, exclusion='exterior')
channel.set_render(color='#333333', linewidth=2.5, alpha=0.9)

box = SquareObstacle(side=box_side, x0=box_cx, y0=box_cy, exclusion='interior')
box.set_render(color='#aa4444', linewidth=2.0, alpha=0.85, fill=True)

sys_b1.add_object(channel)
sys_b1.add_object(box)

# ── Particles ─────────────────────────────────────────────────────────────────
spec = ParticleSpec(count=12, nu=0.5, N_nodes=32)
sys_b1.add_particles(spec)

# ── Initialize (RSA seed + swell to phi_target) ───────────────────────────────
print(f"\nGeometry: Lx={Lx}, Ly={Ly}, box_side={box_side}")
print(f"  box area      = {box_side**2:.2f}")
print(f"  accessible    = {Lx*Ly - box_side**2:.2f}  (out of {Lx*Ly:.2f})")

# phi_target=0.12 is below the RSA-placed phi (~0.13), so swell is skipped.
# n_relax=200 keeps the test fast (B1 only verifies seeding + phi accounting).
sys_b1.initialize(phi_target=0.12, seed=42, verbose=True, n_relax=200)

# ── Verification ──────────────────────────────────────────────────────────────
print("\n--- Verification ---")

# 1. All particle CMs outside the obstacle region
from src.epd.objects import _point_in_polygon
snap = sys_b1.snapshot()
x_cm = snap['x_cm']   # (P, 2)

poly = box.region_polygon(t=0.0)
verts = poly['vertices']

n_inside = 0
for i, cm in enumerate(x_cm):
    if _point_in_polygon(cm, verts):
        n_inside += 1
        print(f"  FAIL: particle {i} CM at {cm} is inside the box!")

if n_inside == 0:
    print(f"  PASS: all {len(x_cm)} particle CMs are outside the box")
else:
    print(f"  FAIL: {n_inside}/{len(x_cm)} particles inside box")

# 2. Phi calculation correctness
accessible = Lx * Ly - box_side**2
phi_outer_corrected = sys_b1.phi_outer

# Naive phi (ignoring obstacle)
from src.epd.initializer import compute_phi_outer
r_c_arr = sys_b1._params['r_c_per_p'].numpy()
phi_naive = compute_phi_outer(sys_b1._state, Lx, Ly, r_c_arr)   # no accessible_area → uses Lx*Ly

print(f"\n  phi_naive     = {phi_naive:.4f}  (denominator = Lx*Ly = {Lx*Ly:.1f})")
print(f"  phi_corrected = {phi_outer_corrected:.4f}  (denominator = accessible = {accessible:.1f})")
print(f"  ratio corrected/naive = {phi_outer_corrected/phi_naive:.4f}  (expected ≈ {Lx*Ly/accessible:.4f})")

ratio_ok = abs(phi_outer_corrected / phi_naive - Lx*Ly/accessible) < 1e-10
if ratio_ok:
    print(f"  PASS: phi ratio matches Lx*Ly/accessible to machine precision")
else:
    print(f"  FAIL: ratio mismatch")

# 3. phi_outer is what system reports at title
print(f"\n  sys_b1.phi_outer = {sys_b1.phi_outer:.4f}")

# ── Snapshot ──────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(-0.5, Lx + 0.5)
ax.set_ylim(-0.5, Ly + 0.5)
ax.set_aspect('equal')
ax.set_title(f'Test B1 — seed φ={phi_outer_corrected:.3f}  '
             f'(accessible = {accessible:.0f} / {Lx*Ly:.0f})',
             fontsize=10)
ax.set_xlabel('x'); ax.set_ylabel('y')

# Draw obstacle box (filled)
from matplotlib.patches import Polygon as MplPolygon
box_patch = MplPolygon(verts, closed=True, facecolor='#ffdddd', edgecolor='#aa4444', lw=2)
ax.add_patch(box_patch)

# Draw channel walls
from src.epd.objects import _rotate_translate
for d in channel.resolved(0.0):
    prim = d['prim']
    ax.plot([prim.p0[0], prim.p1[0]], [prim.p0[1], prim.p1[1]],
            color='#333333', lw=2.5)

# Draw particles
colors = sys_b1._particle_colors
x_all_np = snap['x_all']   # (P, N, 2)
for i, (xy, cm) in enumerate(zip(x_all_np, x_cm)):
    c = colors[i] if i < len(colors) else 'royalblue'
    ax.fill(xy[:, 0], xy[:, 1], alpha=0.55, color=c)
    ax.plot(np.append(xy[:, 0], xy[0, 0]),
            np.append(xy[:, 1], xy[0, 1]), color=c, lw=0.8)

fig.tight_layout()
fig.savefig(args.out, dpi=120)
print(f"\nSnapshot saved → {args.out}")

# Summary
all_pass = (n_inside == 0) and ratio_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test B1 {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
