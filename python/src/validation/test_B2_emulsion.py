"""
test_B2_emulsion.py — Test B2 with emulsion droplets

Same geometry as test_B2_spin.py (spinning+translating SquareObstacle, periodic_x)
but with emulsion particles (surface tension, no bending stiffness).

PASS criteria:
  - phi_outer >= 0.77 after swell
  - 0 droplet CMs inside obstacle at every recorded frame
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Channel, SquareObstacle, _point_in_polygon
from src.epd.motion import MotionSpec
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=str, default='results/test_B2_emulsion.gif')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test B2 (emulsion): SquareObstacle spin + translation")
print("=" * 60)

Lx, Ly    = 20.0, 20.0
box_side  = 5.0
box_cx    = Lx / 2.0
box_cy    = Ly / 2.0
omega_box = 0.2
vx_box    = 0.5 / 3.0
vy_box    = 0.3 / 3.0

sys_em = System(Lx, Ly, periodic_x=True)

channel = Channel(width=Lx, height=Ly, x0=Lx/2, y0=Ly/2, exclusion='exterior')
channel.set_render(color='#333333', linewidth=2.5, alpha=0.9)

box = SquareObstacle(side=box_side, x0=box_cx, y0=box_cy, exclusion='interior')
box.set_motion(MotionSpec(omega=omega_box, vx=vx_box, vy=vy_box,
                          r_ref=(box_cx, box_cy)))
box.set_render(color='#aa4444', linewidth=2.0, alpha=0.85, fill=True)

sys_em.add_object(channel)
sys_em.add_object(box)

# Emulsion droplets: surface tension γ=1, κ=0.2 (q=5), no bending
spec = ParticleSpec(count=25, type='emulsion', gamma=1.0, kappa=0.2,
                    N_nodes=48, poly_dist=0.05)
sys_em.add_particles(spec)

d = spec.derived
print(f"\nEmulsion params: γ={d['gamma']:.2f}  κ={d['kappa']:.3f}"
      f"  K_area={d['K_area']:.3f}  C={d['C']:.1f}  α={d['alpha']:.1f}")

sys_em.initialize(phi_target=0.78, seed=7, verbose=True, n_relax=200,
                  swell_alpha=10.0)

print(f"\nPost-swell: Lx={sys_em.Lx:.3f}  phi_outer={sys_em.phi_outer:.4f}")

dt = sys_em._dt
n_steps      = 2000
sample_every = max(1, n_steps // 80)
print(f"\ndt={dt:.6f}, n_steps={n_steps}, sample_every={sample_every}")

sys_em.run(n_steps, sample_every=sample_every, verbose=True)

t_final = sys_em.t
print(f"Simulation ended at t={t_final:.4f} s")

# Check: all CMs outside obstacle at every frame
n_violations = 0
for fr in sys_em.frames:
    verts = box.region_polygon(fr['t'])['vertices']
    for cm in fr['x_cm']:
        if _point_in_polygon(cm, verts):
            n_violations += 1

phi_ok  = sys_em.phi_outer >= 0.77
excl_ok = n_violations == 0
print(f"\n── Physics check ──")
print(f"  phi_outer = {sys_em.phi_outer:.4f}   {'PASS' if phi_ok else 'FAIL'}")
print(f"  CM violations: {n_violations}   {'PASS' if excl_ok else 'FAIL'}")

# Movie
print("\nRendering movie…")
sys_em.make_movie(
    args.out, fps=12,
    xlim=(-0.5, sys_em.Lx + 0.5),
    ylim=(-0.5, sys_em.Ly + 0.5),
    title=f'Emulsion + spinning obstacle  (γ=1, κ=0.2, φ={sys_em.phi_outer:.3f})',
)
print(f"Movie saved → {args.out}")

all_pass = phi_ok and excl_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test B2 (emulsion) {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
