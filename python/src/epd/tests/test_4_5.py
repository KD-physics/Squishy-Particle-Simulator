"""
test_4_5.py — Phase 4.5: system.py + checkpoint.py verification.

Tests
-----
1. Bit-identical resume
   Initialize P=8, phi_target=0.32 (fast swell).
   Run 100 steps. Save checkpoint. Load. Run 100 more steps.
   Compare to reference that ran 200 steps uninterrupted.
   Gate: max|Δx_cm| < 1e-12

2. Sine-wave wall motion serialization
   System with sine-wave MotionSpec on a Box. Save/load. Verify wall velocity.
   Gate: |v_wall - expected| < 1e-10

3. system.t round-trip
   t before and after save/load must be identical float64 values.
   Gate: t_before == t_after exactly

Run: python src/epd/tests/test_4_5.py
Output: results/phase45_system/resume_test.png
"""

import sys, os, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)

os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
import src.simulation.tf_sim as tf_sim_mod
import tensorflow as tf
tf_sim_mod.set_dtype(tf.float64)

from src.epd.system import System
from src.epd.particles import ParticleSpec
from src.epd.motion import MotionSpec
from src.epd.objects import Box

OUTDIR   = os.path.join(ROOT, 'results', 'phase45_system')
CKPT_DIR = os.path.join(OUTDIR, 'checkpoint_test')
os.makedirs(OUTDIR, exist_ok=True)

results = []

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# Shared fast-swell params (small n_relax for test speed)
# phi from ~0.25 RSA → 0.32: only 7 compression steps at dphi=0.010
SWELL_KW = dict(dphi_init=0.010, dphi_max=0.015, n_relax=50,
                max_extra_relax=100)
# All tests use same P/N to reuse TF JIT graph after first compilation
P_ALL = 8
N_ALL = 32
NU    = 0.5

def _spec():
    return ParticleSpec(count=P_ALL, nu=NU, N_nodes=N_ALL)


# ══════════════════════════════════════════════════════════════════════════════
# Test 1 — bit-identical resume
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 1: Bit-identical resume P=8, phi=0.32 ──────────────────────────")

SEED      = 7
PHI_TGT   = 0.32
N_STEPS_A = 100
N_STEPS_B = 100

# Reference: 200 steps uninterrupted
print("  Building reference system (200 steps)...")
sys_ref = System(Lx=10.0, Ly=10.0)
sys_ref.add_particles(_spec())
sys_ref.initialize(phi_target=PHI_TGT, seed=SEED, verbose=False, **SWELL_KW)
sys_ref.step(N_STEPS_A + N_STEPS_B)
x_cm_ref = sys_ref.state['x_cm'].numpy().copy()
t_ref     = sys_ref.t

# Interrupted: 100, save, load, 100
print("  Building interrupted system (100 + save/load + 100)...")
sys_int = System(Lx=10.0, Ly=10.0)
sys_int.add_particles(_spec())
sys_int.initialize(phi_target=PHI_TGT, seed=SEED, verbose=False, **SWELL_KW)
sys_int.step(N_STEPS_A)

ckpt_path = os.path.join(CKPT_DIR, 'test1')
sys_int.save(ckpt_path)
print(f"  Saved to {ckpt_path}")

sys_loaded = System.from_file(ckpt_path)
print(f"  Loaded: t={sys_loaded.t:.6f}, step={sys_loaded.step_count}")
sys_loaded.step(N_STEPS_B)

x_cm_loaded = sys_loaded.state['x_cm'].numpy()
max_dx = float(np.max(np.abs(x_cm_loaded - x_cm_ref)))

check("1.1 Bit-identical x_cm after resume",
      max_dx < 1e-12,
      f"max|Δx_cm|={max_dx:.3e}")
check("1.2 Step count agrees",
      sys_loaded.step_count == N_STEPS_A + N_STEPS_B,
      f"got={sys_loaded.step_count} expected={N_STEPS_A + N_STEPS_B}")
check("1.3 Time agrees",
      abs(sys_loaded.t - t_ref) < 1e-12,
      f"loaded={sys_loaded.t:.10f} ref={t_ref:.10f}")
check("1.4 phi_outer in range",
      0.20 < sys_loaded.phi_outer < 1.0,
      f"phi_outer={sys_loaded.phi_outer:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 2 — MotionSpec round-trip (parametric sine wave wall)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 2: MotionSpec round-trip (sine-wave wall velocity) ─────────────")

freq          = 1.0
amp           = 0.5
t_quarter     = math.pi / (2.0 * freq)
expected_vx_0 = amp * math.cos(0.0)
expected_vx_q = amp * math.cos(freq * t_quarter)

ms_wall  = MotionSpec(vx_ac=amp, freq_x=freq, phase_x=0.0)

# Build a periodic system, initialize it, then add a Wall with motion before save.
# (Adding the wall after initialize avoids the RSA domain issue.)
sys2 = System(Lx=10.0, Ly=10.0)
sys2.add_particles(_spec())
sys2.initialize(phi_target=0.32, seed=42, verbose=False, **SWELL_KW)

# Add a bottom wall with sine-wave motion (after swell; won't affect RSA)
from src.epd.objects import Wall
w = Wall(p0=(0.0, 0.0), p1=(10.0, 0.0), normal=np.array([0.0, 1.0]))
w.set_motion(ms_wall)
sys2._objects = [w]

ckpt2 = os.path.join(CKPT_DIR, 'test2')
sys2.save(ckpt2)

from src.epd.checkpoint import load_checkpoint
sys2_loaded = load_checkpoint(ckpt2)

check("2.1 Object count preserved",
      len(sys2_loaded._objects) == 1,
      f"n_objects={len(sys2_loaded._objects)}")

ms_loaded = sys2_loaded._objects[0]._motion if sys2_loaded._objects else None
if ms_loaded is not None:
    vx0, _, _ = ms_loaded.velocity(0.0)
    vxq, _, _ = ms_loaded.velocity(t_quarter)
    check("2.2 Wall vx at t=0 correct",
          abs(vx0 - expected_vx_0) < 1e-10,
          f"vx={vx0:.6f} expected={expected_vx_0:.6f}")
    check("2.3 Wall vx at t=T/4 correct",
          abs(vxq - expected_vx_q) < 1e-10,
          f"vx={vxq:.6f} expected={expected_vx_q:.6f}")
else:
    check("2.2 Wall motion loaded (not None)", False, "ms_loaded is None")
    check("2.3 Wall vx at t=T/4 correct",     False, "ms_loaded is None")


# ══════════════════════════════════════════════════════════════════════════════
# Test 3 — system.t float64 round-trip
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 3: system.t float64 round-trip ─────────────────────────────────")

sys_t = System(Lx=8.0, Ly=8.0)
sys_t.add_particles(_spec())
sys_t.initialize(phi_target=0.28, seed=11, verbose=False, **SWELL_KW)
sys_t.step(37)
t_before    = sys_t.t
step_before = sys_t.step_count

ckpt3 = os.path.join(CKPT_DIR, 'test3')
sys_t.save(ckpt3)
sys_t3 = System.from_file(ckpt3)

check("3.1 system.t bit-identical after save/load",
      sys_t3.t == t_before,
      f"before={t_before:.15e} after={sys_t3.t:.15e}")
check("3.2 step_count preserved",
      sys_t3.step_count == step_before,
      f"before={step_before} after={sys_t3.step_count}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary figure
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

def _draw(ax, sys_obj, title):
    ax.set_aspect('equal')
    ax.set_facecolor('#f5f5f5')
    ax.set_xlim(-0.2, sys_obj.Lx + 0.2)
    ax.set_ylim(-0.2, sys_obj.Ly + 0.2)
    x_all = sys_obj.state['x_all'].numpy()
    P     = x_all.shape[0]
    colors = plt.cm.tab20(np.linspace(0, 1, P))
    for i in range(P):
        xy = x_all[i]
        ax.fill(xy[:, 0], xy[:, 1], color=colors[i], alpha=0.55)
        ax.plot(np.append(xy[:, 0], xy[0, 0]),
                np.append(xy[:, 1], xy[0, 1]), 'k-', lw=0.4)
    ax.plot([0, sys_obj.Lx, sys_obj.Lx, 0, 0],
            [0, 0, sys_obj.Ly, sys_obj.Ly, 0], 'b--', lw=1.0, alpha=0.5)
    ax.set_title(title)
    ax.grid(True, alpha=0.15)

_draw(axes[0], sys_ref,    f'Reference (200 steps)\nφ={sys_ref.phi_outer:.4f}')
_draw(axes[1], sys_loaded, f'Resumed (100+100 steps)\nmax|Δx|={max_dx:.1e}')

plt.tight_layout()
outpath = os.path.join(OUTDIR, 'resume_test.png')
plt.savefig(outpath, dpi=100)
plt.close(fig)
print(f"\nPlot saved: {outpath}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
n_pass = sum(r[1] for r in results)
n_fail = sum(not r[1] for r in results)
print(f"\n{'='*60}")
print(f"TOTAL: {n_pass}/{len(results)} PASS")
if n_fail:
    print("FAILED tests:")
    for name, passed, detail in results:
        if not passed:
            print(f"  - {name}  ({detail})")

sys.exit(0 if n_fail == 0 else 1)
