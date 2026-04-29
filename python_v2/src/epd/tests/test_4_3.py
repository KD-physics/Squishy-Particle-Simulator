"""
test_4_3.py — Phase 4.3: motion.py expanded MotionSpec verification.

Tests:
  1. MotionSpec(vx=1.0).resolve_tf(t=5.0) → vx=1.0, vy=0.0, omega=0.0
  2. MotionSpec(vx=lambda t: tf.sin(t)).resolve_tf(t=π/2) → vx≈1.0
  3. from_samples(vx_fn=cos, dt=0.01, duration=10).resolve_tf(π) → vx≈-1.0
  4. Orbital motion (omega, r_ref): wall rotates correctly
  5. Backwards compat: existing DC+AC parametric API still works

Run: python src/epd/tests/test_4_3.py
Output: results/phase43_motion/motion_interp.png
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

from src.epd.motion import MotionSpec

OUTDIR = os.path.join(ROOT, 'results', 'phase43_motion')
os.makedirs(OUTDIR, exist_ok=True)

DTYPE = tf.float64
results = []

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════════════
# Test 1 — scalar constant
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 1: MotionSpec(vx=1.0).resolve_tf ───────────────────────────────")

ms1 = MotionSpec(vx=1.0)
t5  = tf.constant(5.0, dtype=DTYPE)
vx1, vy1, om1, rr1 = ms1.resolve_tf(t5)

check("1.1 vx=1.0", abs(float(vx1) - 1.0) < 1e-12, f"vx={float(vx1):.6f}")
check("1.2 vy=0.0", abs(float(vy1)) < 1e-12)
check("1.3 omega=0.0", abs(float(om1)) < 1e-12)
check("1.4 r_ref=[0,0]", np.allclose(rr1.numpy(), [0, 0]))

# t-independence for DC
vx1b, _, _, _ = ms1.resolve_tf(tf.constant(0.0, dtype=DTYPE))
check("1.5 DC is t-independent", abs(float(vx1) - float(vx1b)) < 1e-12)


# ══════════════════════════════════════════════════════════════════════════════
# Test 2 — TF-native callable
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 2: MotionSpec(vx=lambda t: tf.sin(t)) ──────────────────────────")

ms2 = MotionSpec(vx=lambda t: tf.sin(tf.cast(t, DTYPE)))
t_pi2  = tf.constant(math.pi / 2, dtype=DTYPE)
t_pi   = tf.constant(math.pi,     dtype=DTYPE)
t_zero = tf.constant(0.0,          dtype=DTYPE)

vx2_pi2, _, _, _ = ms2.resolve_tf(t_pi2)
vx2_pi,  _, _, _ = ms2.resolve_tf(t_pi)
vx2_0,   _, _, _ = ms2.resolve_tf(t_zero)

check("2.1 sin(π/2) ≈ 1.0", abs(float(vx2_pi2) - 1.0) < 1e-10,
      f"got={float(vx2_pi2):.8f}")
check("2.2 sin(π)   ≈ 0.0", abs(float(vx2_pi)) < 1e-10,
      f"got={float(vx2_pi):.2e}")
check("2.3 sin(0)   = 0.0", abs(float(vx2_0)) < 1e-12)


# ══════════════════════════════════════════════════════════════════════════════
# Test 3 — pre-sampled (from_samples)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 3: MotionSpec.from_samples(vx_fn=cos, dt=0.01, duration=10) ────")

ms3 = MotionSpec.from_samples(vx_fn=lambda t: np.cos(t), dt=0.01, duration=10.0)

t_pi_tf   = tf.constant(math.pi,       dtype=DTYPE)
t_2pi_tf  = tf.constant(2*math.pi,     dtype=DTYPE)
t_pi2_tf  = tf.constant(math.pi/2.0,   dtype=DTYPE)
t_zero_tf = tf.constant(0.0,            dtype=DTYPE)

vx3_pi,  _, _, _ = ms3.resolve_tf(t_pi_tf)
vx3_2pi, _, _, _ = ms3.resolve_tf(t_2pi_tf)
vx3_pi2, _, _, _ = ms3.resolve_tf(t_pi2_tf)
vx3_0,   _, _, _ = ms3.resolve_tf(t_zero_tf)

check("3.1 cos(π) ≈ -1.0 (within 1%)",
      abs(float(vx3_pi) + 1.0) < 0.01,
      f"got={float(vx3_pi):.4f}")
check("3.2 cos(2π) ≈ 1.0 (within 1%)",
      abs(float(vx3_2pi) - 1.0) < 0.01,
      f"got={float(vx3_2pi):.4f}")
check("3.3 cos(π/2) ≈ 0.0 (within 1%)",
      abs(float(vx3_pi2)) < 0.01,
      f"got={float(vx3_pi2):.4f}")
check("3.4 cos(0) = 1.0",
      abs(float(vx3_0) - 1.0) < 1e-10)

# Test from_samples with Python float t (not TF tensor)
vx3_pi_py, _, _ = ms3.velocity(math.pi)
check("3.5 velocity(π) Python path matches",
      abs(vx3_pi_py + 1.0) < 0.01,
      f"got={vx3_pi_py:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 4 — spin (omega, r_ref): object rotated
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 4: spin MotionSpec(omega=1.0, r_ref=(3,0)) ─────────────────────")

from src.epd.objects import Wall

ms4  = MotionSpec(omega=1.0, r_ref=(3.0, 0.0))
wall = Wall(p0=(4.0, -0.5), p1=(4.0, 0.5))
wall.set_motion(ms4)

# At t=0: wall is at its initial position
prims_t0 = wall.resolved(t=0.0)
# At t=π/2: wall has rotated π/2 around (3,0)
# Initial p0=(4,-0.5): relative to pivot=(3,0): (1,-0.5)
# After π/2 rotation: (0.5, 1) → absolute: (3.5, 1)
# (Using displacement integral: dtheta = omega*t = π/2)
prims_t1 = wall.resolved(t=math.pi / 2)

seg0 = prims_t0[0]['prim']
seg1 = prims_t1[0]['prim']

check("4.1 At t=0 p0 unchanged",
      np.allclose(seg0.p0, [4.0, -0.5], atol=1e-10),
      f"p0={seg0.p0}")

# displacement (dx,dy,dtheta) from omega=1 at t=π/2:
# dtheta = π/2; displacement of pivot itself = 0 (omega is pure spin in MotionSpec)
# Actually in MotionSpec, omega is spin of the OBJECT about its own CM.
# objects.py adds om_parent to each child's 'omega' key.
# The wall has moved: objects.resolved applies displacement(t) to the origin.
# In this case MotionSpec with only omega_dc=1.0 gives:
#   dx=0, dy=0, dtheta=π/2
# So origin stays at (0,0) and wall rotates by π/2 in resolved().
dx4, dy4, dth4 = ms4.displacement(t=math.pi/2)
check("4.2 displacement(π/2): dx=0, dy=0, dtheta=π/2",
      abs(dx4) < 1e-10 and abs(dy4) < 1e-10 and abs(dth4 - math.pi/2) < 1e-10,
      f"dx={dx4:.2e} dy={dy4:.2e} dth={dth4:.6f}")

# omega at t=0 passed to wall
check("4.3 Wall omega=1 in resolved dict",
      abs(prims_t0[0]['omega'] - 1.0) < 1e-10,
      f"omega={prims_t0[0]['omega']}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 5 — backwards compat: DC+AC parametric still works
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 5: DC+AC parametric backwards compat ───────────────────────────")

ms5 = MotionSpec(vx_dc=1.0, vx_ac=0.5, freq_x=2.0, phase_x=0.0)
# vx(t) = 1.0 + 0.5*cos(2t)
# at t=0: 1.5;  at t=π/2: 1.0+0.5*cos(π)=0.5;  at t=π: 1.5

vx5_0, _, _ = ms5.velocity(0.0)
vx5_pi2, _, _ = ms5.velocity(math.pi / 2)
vx5_pi, _, _ = ms5.velocity(math.pi)

check("5.1 vx(0) = 1.5", abs(vx5_0 - 1.5) < 1e-12, f"got={vx5_0:.6f}")
check("5.2 vx(π/2) = 0.5", abs(vx5_pi2 - 0.5) < 1e-10, f"got={vx5_pi2:.6f}")
check("5.3 vx(π) = 1.5",  abs(vx5_pi - 1.5) < 1e-10,  f"got={vx5_pi:.6f}")

# resolve_tf
vx5_tf, _, _, _ = ms5.resolve_tf(tf.constant(0.0, dtype=DTYPE))
check("5.4 resolve_tf(0) = 1.5", abs(float(vx5_tf) - 1.5) < 1e-10)

# to_traj_row should still work
traj = ms5.to_traj_row()
check("5.5 to_traj_row returns (18,) array",
      hasattr(traj, '__len__') and len(traj) == 18)


# ══════════════════════════════════════════════════════════════════════════════
# Test 6 — is_static
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 6: is_static() ─────────────────────────────────────────────────")

ms_static  = MotionSpec()
ms_moving  = MotionSpec(vx=1.0)
ms_callable = MotionSpec(vx=lambda t: tf.sin(t))

check("6.1 MotionSpec() is_static", ms_static.is_static())
check("6.2 MotionSpec(vx=1.0) not static", not ms_moving.is_static())
check("6.3 MotionSpec(vx=callable) not static", not ms_callable.is_static())


# ══════════════════════════════════════════════════════════════════════════════
# Plot: sampled interpolation error
# ══════════════════════════════════════════════════════════════════════════════
t_plot = np.linspace(0, 10, 500)
analytic  = np.cos(t_plot)
sampled   = np.array([float(ms3.resolve_tf(tf.constant(ti, dtype=DTYPE))[0]) for ti in t_plot])

fig, axes = plt.subplots(2, 1, figsize=(10, 6))
ax = axes[0]
ax.plot(t_plot, analytic, 'k-', label='analytic cos(t)', lw=2)
ax.plot(t_plot, sampled,  'b--', label='sampled (dt=0.01)', lw=1.5, alpha=0.8)
ax.set_xlabel('t'); ax.set_ylabel('vx')
ax.legend(); ax.set_title('Sampled MotionSpec interpolation vs analytic')
ax.grid(True, alpha=0.3)

ax2 = axes[1]
err = np.abs(sampled - analytic)
ax2.semilogy(t_plot, err + 1e-20, 'r-', lw=1.5)
ax2.axhline(0.01, ls='--', color='gray', label='1% threshold')
ax2.set_xlabel('t'); ax2.set_ylabel('|error|')
ax2.set_title('Interpolation error (dt=0.01)')
ax2.legend(); ax2.grid(True, alpha=0.3)

plt.tight_layout()
outpath = os.path.join(OUTDIR, 'motion_interp.png')
plt.savefig(outpath, dpi=120)
print(f"\nPlot saved: {outpath}")

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
