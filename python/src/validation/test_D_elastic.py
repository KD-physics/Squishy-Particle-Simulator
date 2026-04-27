"""
test_D_elastic.py — Test D (elastic): 40 elastic capsules falling under gravity.

Same geometry as test_D_falling_emulsion.py but with elastic particles (ν=0.5).
10k-step spot check only.

Usage:
    python src/validation/test_D_elastic.py [--out PATH]
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
parser.add_argument('--out', type=str, default='results/test_D_elastic.gif')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test D (elastic): 40 elastic capsules falling under gravity")
print("=" * 60)

NU       = 0.5
N_NODES  = 36
N_DROP   = 40
POLY     = 0.05
G_SIM    = 0.05

LX       = 12.0
LY       = 60.0

STEPS_PER_BATCH = 5_000
N_BATCHES       = 2
SAMPLE_EVERY    = 500

print(f"\nElastic: ν={NU}  N={N_NODES}  count={N_DROP}  g={G_SIM}")

sys_e = System(LX, LY, periodic_x=False, periodic_y=False, g=G_SIM)

for w_pts, normal in [
    (((0, 0), (LX, 0)),  (0,  1)),
    (((0, 0), (0,  LY)), (1,  0)),
    (((LX,0), (LX, LY)),(-1,  0)),
]:
    w = Wall(*w_pts, normal=normal)
    w.set_render(color='#333333', linewidth=2.0, alpha=0.9)
    sys_e.add_object(w)

spec = ParticleSpec(count=N_DROP, type='elastic',
                    nu=NU, N_nodes=N_NODES, poly_dist=POLY)
sys_e.add_particles(spec)

d = spec.derived
print(f"Elastic params: ν={NU}  q={d['q']:.3f}  El_t={d['El_t']:.3f}"
      f"  K_area={d['K_area']:.3f}  C={d['C']:.1f}  α={d['alpha']:.3f}")

sys_e.initialize(phi_target=0.80, seed=42, verbose=True,
                 relax_only=True, n_relax_init=200)

initial_y_mean = float(sys_e.state['x_cm'].numpy()[:, 1].mean())
print(f"\nInitial mean y_cm = {initial_y_mean:.3f}")
print(f"dt = {sys_e._dt:.4e}   alpha = {sys_e._alpha_damp:.4f}")
print(f"Terminal velocity estimate: v_t = g/α ≈ {G_SIM/sys_e._alpha_damp:.3f}")

print(f"\nStarting: {N_BATCHES} × {STEPS_PER_BATCH:,} steps")

for batch in range(1, N_BATCHES + 1):
    print(f"── Batch {batch}/{N_BATCHES}  (steps {(batch-1)*STEPS_PER_BATCH:,}–{batch*STEPS_PER_BATCH:,}) ──")
    sys_e.run(STEPS_PER_BATCH, sample_every=SAMPLE_EVERY, verbose=True)

    cur_y   = float(sys_e.frames[-1]['x_cm'][:, 1].mean())
    cur_cms = sys_e.frames[-1]['x_cm']
    fl_viol = sum(1 for cm in cur_cms if cm[1] < 0.5)
    si_viol = sum(1 for cm in cur_cms if cm[0] < 0.0 or cm[0] > LX)
    ratio   = cur_y / initial_y_mean
    print(f"  mean y_cm = {cur_y:.3f}  ({ratio:.2f}× initial)  "
          f"floor_viol={fl_viol}  side_viol={si_viol}  frames={len(sys_e.frames)}")

    sys_e.make_movie(
        args.out, fps=10,
        xlim=(-0.5, LX + 0.5),
        ylim=(-0.5, 24.5),
        title=f'Falling elastic capsules  ν={NU}  g={G_SIM}  batch {batch}/{N_BATCHES}',
    )
    print(f"  GIF updated ({len(sys_e.frames)} frames)")

# Physics check
final_y_mean = float(sys_e.frames[-1]['x_cm'][:, 1].mean())
floor_violations = sum(1 for cm in sys_e.frames[-1]['x_cm'] if cm[1] < 0.5)
side_violations  = sum(1 for cm in sys_e.frames[-1]['x_cm'] if cm[0] < 0.0 or cm[0] > LX)
settled  = final_y_mean < 0.95 * initial_y_mean  # just 5% drop for a 10k spot check

print(f"\n── Physics check ──")
print(f"  Initial mean y_cm  = {initial_y_mean:.3f}")
print(f"  Final   mean y_cm  = {final_y_mean:.3f}   ({'PASS' if settled else 'not yet settled (spot check only)'})")
print(f"  Floor violations   = {floor_violations}   {'PASS' if floor_violations==0 else 'FAIL'}")
print(f"  Sidewall violations= {side_violations}   {'PASS' if side_violations==0 else 'FAIL'}")

print(f"\nMovie saved → {args.out}")
