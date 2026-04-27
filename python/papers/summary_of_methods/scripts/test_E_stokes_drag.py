"""
test_E_stokes_drag.py — Phase 5 benchmark: Stokes drag + parameter update

Tests
-----
E1 — Terminal velocity (emulsion + elastic)
     Single particle under gravity; v_cm plateaus at g/(2·Oh) for emulsion.

E2 — Constant background flow tracking
     Particle at rest, U_bg = constant; v_cm → U_bg exponentially.

E3 — Shear flow drift and rotation
     U_bg = (rate·y, 0); CM drifts, particle rotates at Ω ≈ rate/2.

E4 — Post-init set_param update
     Run to terminal velocity, call sys.set_param('Oh', new), confirm new v_t.

Usage:
    python src/validation/test_E_stokes_drag.py [--test E1 E2 E3 E4]
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
from src.epd.drag import oh_from_terminal_velocity

os.makedirs('results', exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('--test', nargs='+', default=['E1','E2','E3','E4'],
                    help='Tests to run (default: all)')
args = parser.parse_args()

TESTS = [t.upper() for t in args.test]
results = {}

# ── helpers ────────────────────────────────────────────────────────────────────

def _build_single_drop(Oh=None, g=0.05, N=36, particle_type='emulsion', nu=0.5,
                        LX=20.0, LY=40.0, with_floor=True):
    """Build a single-particle system in a tall box."""
    sys_ = System(LX, LY, periodic_x=False, periodic_y=False, g=g)
    if with_floor:
        w = Wall((0, 0), (LX, 0), normal=(0, 1))
        w.set_render(color='#333333', linewidth=2.0)
        sys_.add_object(w)
    if particle_type == 'emulsion':
        spec = ParticleSpec(count=1, type='emulsion', gamma=1.0, kappa=0.02,
                            N_nodes=N, Oh=Oh)
    else:
        spec = ParticleSpec(count=1, type='elastic', nu=nu,
                            N_nodes=N, Oh=Oh)
    sys_.add_particles(spec)
    sys_.initialize(phi_target=0.80, seed=0, verbose=False,
                    relax_only=True, n_relax_init=0)
    return sys_, spec


def _plateau_velocity(sys_, n_run=6000, sample=200):
    """Run until v_cm_y plateaus; return final mean over last 10 samples."""
    sys_.run(n_run, sample_every=sample, verbose=False)
    vy_series = np.array([fr['x_cm'][0, 1] for fr in sys_.frames])
    # finite-difference velocity from CM y positions
    dt_sample = sample * sys_._dt
    if len(vy_series) >= 4:
        vy_vals = np.diff(vy_series[-6:]) / dt_sample
        return float(np.mean(vy_vals))
    return 0.0


# ══════════════════════════════════════════════════════════════════════════════
# Test E1 — Terminal velocity
# ══════════════════════════════════════════════════════════════════════════════

if 'E1' in TESTS:
    print("=" * 60)
    print("Test E1: Terminal velocity")
    print("=" * 60)

    cases = [
        dict(label='emul Oh=0.25 g=0.05', type='emulsion', Oh=0.25, g=0.05),
        dict(label='emul Oh=0.50 g=0.05', type='emulsion', Oh=0.50, g=0.05),
        dict(label='emul Oh=0.50 g=0.10', type='emulsion', Oh=0.50, g=0.10),
        dict(label='elas Oh=0.50 g=0.05', type='elastic',  Oh=0.50, g=0.05, nu=0.5),
    ]

    e1_pass = True
    for c in cases:
        g   = c['g']
        Oh  = c['Oh']
        ptype = c['type']
        sys_, spec = _build_single_drop(Oh=Oh, g=g, particle_type=ptype,
                                         nu=c.get('nu', 0.5),
                                         LX=20.0, LY=60.0)
        # Place particle near top so it can fall
        x_cm_np = sys_.state['x_cm'].numpy()
        x_all_np = sys_.state['x_all'].numpy()
        dy = 30.0 - float(x_cm_np[0, 1])
        x_cm_np[0, 1] += dy
        x_all_np[0, :, 1] += dy
        sys_._state['x_cm'] = tf.constant(x_cm_np, dtype=tf.float64)
        sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

        v_meas = _plateau_velocity(sys_, n_run=8000, sample=200)
        v_meas = abs(v_meas)   # falling = negative y, take magnitude

        d = spec.derived
        xi = d['xi']
        v_theory = spec.terminal_velocity(g)
        err = abs(v_meas - v_theory) / max(v_theory, 1e-9)
        passed = err < 0.15   # 15% tolerance — initial transient included
        e1_pass = e1_pass and passed
        print(f"  {c['label']:30s}  v_theory={v_theory:.4f}  v_meas={v_meas:.4f}"
              f"  err={err:.1%}  {'PASS' if passed else 'FAIL'}")

    results['E1'] = e1_pass
    print(f"\nE1: {'PASS' if e1_pass else 'FAIL'}\n")


# ══════════════════════════════════════════════════════════════════════════════
# Test E2 — Constant background flow tracking
# ══════════════════════════════════════════════════════════════════════════════

if 'E2' in TESTS:
    print("=" * 60)
    print("Test E2: Constant background flow tracking")
    print("=" * 60)

    U0 = 0.10
    Oh = 0.5

    sys_, spec = _build_single_drop(Oh=Oh, g=0.0, with_floor=False,
                                     LX=40.0, LY=40.0)
    # Centre the particle
    x_cm_np  = sys_.state['x_cm'].numpy()
    x_all_np = sys_.state['x_all'].numpy()
    dx = 20.0 - float(x_cm_np[0, 0])
    dy = 20.0 - float(x_cm_np[0, 1])
    x_cm_np[0]  += [dx, dy]
    x_all_np[0] += [dx, dy]
    sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
    sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

    sys_.U_background = ('constant', {'U': (U0, 0.0)})

    # Time constant: tau = M / zeta = rho_d*pi*R0^2 / (xi*2*pi*R0) = R0/(2*xi)
    xi = spec.derived['xi']
    tau_drag = spec.R0_mean / (2.0 * xi)
    n_tau = int(5.0 * tau_drag / sys_._dt)
    n_run = max(n_tau, 4000)

    sys_.run(n_run, sample_every=max(1, n_run // 80), verbose=False)

    vx_final = float(np.diff([fr['x_cm'][0, 0] for fr in sys_.frames[-4:]])[-1]) / \
               (max(1, n_run // 80) * sys_._dt)
    err = abs(vx_final - U0) / U0
    e2_pass = err < 0.10

    print(f"  U_bg = {U0}   Oh = {Oh}   xi = {xi:.4f}   tau_drag = {tau_drag:.2f}")
    print(f"  v_cm_x final = {vx_final:.5f}   target = {U0}   err = {err:.1%}")
    print(f"  E2: {'PASS' if e2_pass else 'FAIL'}")
    results['E2'] = e2_pass
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Test E3 — Shear flow drift and rotation
# ══════════════════════════════════════════════════════════════════════════════

if 'E3' in TESTS:
    print("=" * 60)
    print("Test E3: Shear flow drift and rotation")
    print("=" * 60)

    shear_rate = 0.05
    Oh = 0.5
    y_cm0 = 20.0

    sys_, spec = _build_single_drop(Oh=Oh, g=0.0, with_floor=False,
                                     LX=200.0, LY=40.0)
    # Place particle at y = y_cm0, x = 100
    x_cm_np  = sys_.state['x_cm'].numpy()
    x_all_np = sys_.state['x_all'].numpy()
    dx = 100.0 - float(x_cm_np[0, 0])
    dy = y_cm0  - float(x_cm_np[0, 1])
    x_cm_np[0]  += [dx, dy]
    x_all_np[0] += [dx, dy]
    sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
    sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

    sys_.U_background = ('shear', {'rate': shear_rate})

    n_run = 4000
    sample = 100
    sys_.run(n_run, sample_every=sample, verbose=False)

    # Measure drift rate in x
    x_cms = np.array([fr['x_cm'][0, 0] for fr in sys_.frames])
    t_arr = np.array([fr['t'] for fr in sys_.frames])
    if len(x_cms) > 5:
        coeffs = np.polyfit(t_arr[3:], x_cms[3:], 1)
        dx_dt_meas = coeffs[0]
    else:
        dx_dt_meas = 0.0

    # Expected drift: U_x = shear_rate * y_cm (particle advects with local flow)
    y_cm_final = float(sys_.frames[-1]['x_cm'][0, 1])
    U_x_theory = shear_rate * y_cm0   # y_cm should not change much

    drift_err = abs(dx_dt_meas - U_x_theory) / max(abs(U_x_theory), 1e-9)

    # Measure rotation rate
    omega_vals = np.array([sys_.state['omega'].numpy()[0]
                           if fr is sys_.frames[-1] else 0.0
                           for fr in sys_.frames])
    omega_meas = float(sys_._state['omega'].numpy()[0])
    # Shear U_x = rate*y → vorticity = -rate → rigid-body rotation = -rate/2
    omega_theory = -shear_rate / 2.0
    rot_err = abs(omega_meas - omega_theory) / max(abs(omega_theory), 1e-9)

    drift_ok = drift_err < 0.15
    rot_ok   = rot_err   < 0.30   # looser: elastic modes contribute
    e3_pass  = drift_ok and rot_ok

    print(f"  shear_rate = {shear_rate}   Oh = {Oh}")
    print(f"  Drift: measured {dx_dt_meas:.5f}  theory {U_x_theory:.5f}  err {drift_err:.1%}  {'PASS' if drift_ok else 'FAIL'}")
    print(f"  Rotation: measured {omega_meas:.5f}  theory {omega_theory:.5f}  err {rot_err:.1%}  {'PASS' if rot_ok else 'FAIL'}")
    print(f"  E3: {'PASS' if e3_pass else 'FAIL'}")
    results['E3'] = e3_pass
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Test E4 — Post-init set_param update
# ══════════════════════════════════════════════════════════════════════════════

if 'E4' in TESTS:
    print("=" * 60)
    print("Test E4: Post-init set_param update")
    print("=" * 60)

    Oh1 = 0.25
    Oh2 = 1.00
    g   = 0.05
    v_t1_theory = g / (2.0 * Oh1)   # 0.10
    v_t2_theory = g / (2.0 * Oh2)   # 0.025

    sys_, spec = _build_single_drop(Oh=Oh1, g=g, LX=20.0, LY=80.0)
    # Place particle near top
    x_cm_np  = sys_.state['x_cm'].numpy()
    x_all_np = sys_.state['x_all'].numpy()
    dy = 50.0 - float(x_cm_np[0, 1])
    x_cm_np[0, 1]  += dy
    x_all_np[0, :, 1] += dy
    sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
    sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

    # Run to first terminal velocity
    sys_.run(6000, sample_every=200, verbose=False)
    vy_series1 = np.array([fr['x_cm'][0, 1] for fr in sys_.frames])
    dt_s = 200 * sys_._dt
    v_meas1 = abs(float(np.mean(np.diff(vy_series1[-5:]))) / dt_s)

    err1 = abs(v_meas1 - v_t1_theory) / v_t1_theory

    # Update Oh post-init — no re-initialization
    sys_.set_param('Oh', Oh2, particles='all')

    # Verify TF params updated
    xi_new = spec.derived['xi']
    xi_in_params = float(sys_._params['xi_drag_per_p'].numpy()[0])
    param_updated = abs(xi_new - xi_in_params) < 1e-10

    # Run to second terminal velocity from current state
    sys_.clear_recording()
    sys_.run(6000, sample_every=200, verbose=False)
    vy_series2 = np.array([fr['x_cm'][0, 1] for fr in sys_.frames])
    v_meas2 = abs(float(np.mean(np.diff(vy_series2[-5:]))) / dt_s)

    err2 = abs(v_meas2 - v_t2_theory) / v_t2_theory

    phase1_ok  = err1 < 0.15
    param_ok   = param_updated
    phase2_ok  = err2 < 0.15
    e4_pass    = phase1_ok and param_ok and phase2_ok

    print(f"  Phase 1 (Oh={Oh1}): v_theory={v_t1_theory:.4f}  v_meas={v_meas1:.4f}  err={err1:.1%}  {'PASS' if phase1_ok else 'FAIL'}")
    print(f"  set_param('Oh', {Oh2}) → xi_spec={xi_new:.5f}  xi_params={xi_in_params:.5f}  {'PASS' if param_ok else 'FAIL (params not updated)'}")
    print(f"  Phase 2 (Oh={Oh2}): v_theory={v_t2_theory:.4f}  v_meas={v_meas2:.4f}  err={err2:.1%}  {'PASS' if phase2_ok else 'FAIL'}")
    print(f"  E4: {'PASS' if e4_pass else 'FAIL'}")
    results['E4'] = e4_pass
    print()


# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("SUMMARY")
print("=" * 60)
all_pass = True
for name, passed in results.items():
    print(f"  {name}: {'PASS' if passed else 'FAIL'}")
    all_pass = all_pass and passed

print(f"\n{'PASS' if all_pass else 'FAIL'}: Test E Stokes drag")
if not all_pass:
    sys.exit(1)
