"""
test_A_elastic_api.py — Test A (elastic): Two-disk oscillatory squeeze via System API

Uses identical physics parameters and geometry to test_A_elastic_tf.py (ground truth).
Machine-precision agreement is guaranteed by construction: both call make_state() with
identical particles and step_full_tf() with identical wall positions at every step.

The explicit cross-check runs 10 steps both ways and verifies max|Δx_all| < 1e-12.

Usage:
    python src/validation/test_A_elastic_api.py [--cycles N] [--out PATH]
"""

import sys, os, argparse, time
import numpy as np
import tensorflow as tf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)
from src.simulation.tf_sim import (make_state, step_full_tf, make_prim_data,
                                    DTYPE, NP_DTYPE)
from src.simulation.candidacy_manager import CandidacyManager
from src.simulation.capsule_shell import CapsuleParticle, CapsuleSim
from src.simulation.contact_primitives import LineSegment

from src.epd.particles import ParticleSpec
from src.epd.objects import Wall, OscillatingWall
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--cycles', type=int, default=3)
parser.add_argument('--out',    type=str, default='results/test_A_elastic_api.gif')
args = parser.parse_args()

# ── Parameters (identical to TF test) ────────────────────────────────────────
R0, N, nu = 1.0, 32, 0.5
spec = ParticleSpec(count=2, nu=nu, N_nodes=N)
d    = spec.derived
tau, C, K_area, El_t = d['tau'], d['C'], d['K_area'], d['El_t']
S    = 1.0

p_ref = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area)
L0, r_c = p_ref.L0, p_ref.r_c

c_edge     = np.sqrt(El_t / (p_ref.rho_d * L0))
T_wave     = 2.0 * np.pi * R0 / c_edge
alpha_damp = 2.0 / T_wave

T_osc     = 5.0 * T_wave
omega_osc = 2.0 * np.pi / T_osc
A_wall    = 0.1 * R0

sim_probe = CapsuleSim([p_ref])
dt_max, _ = sim_probe.estimate_dt_max()
dt        = 0.4 * dt_max

n_cycles = args.cycles
t_total  = n_cycles * T_osc
n_steps  = int(np.ceil(t_total / dt))

print("=" * 60)
print(f"Test A (elastic API): Two-disk oscillatory squeeze  (nu={nu})")
print("=" * 60)
print(f"  q={d['q']:.4f}  TAU={d['TAU']:.4f}")
print(f"  tau={tau:.4f}  El_t={El_t:.4f}  K_area={K_area:.4f}")
print(f"  C={C:.2f}  alpha={alpha_damp:.4f}")

# ── Geometry (identical to TF test) ──────────────────────────────────────────
cy_disk  = R0 + L0
half_w   = R0 + 2.0 * L0
y_top0   =  2.0 * R0 + 3.0 * L0
y_bot0   = -(2.0 * R0 + 3.0 * L0)
r_c_wall = L0

print(f"\n  Geometry:  R0={R0}  N={N}  L0={L0:.4f}  r_c={r_c:.4f}")
print(f"  Dynamics:  dt={dt:.5f}  T_wave={T_wave:.4f}  alpha={alpha_damp:.4f}")
print(f"  Oscillation: A={A_wall}  T_osc={T_osc:.4f}  omega={omega_osc:.4f}")
print(f"  Total: {n_cycles} cycles  {n_steps} steps  t_total={t_total:.3f}")
print(f"  Box: half_w={half_w:.4f}  y_top0={y_top0:.4f}")

# ── Build System via API ──────────────────────────────────────────────────────
sys_api = System(Lx=20.0, Ly=20.0, periodic_x=False, periodic_y=False)

top_wall_obj = OscillatingWall(y0=y_top0, half_w=half_w, A=A_wall,
                                omega=omega_osc, is_top=True,  r_c_wall=r_c_wall)
bot_wall_obj = OscillatingWall(y0=y_bot0, half_w=half_w, A=A_wall,
                                omega=omega_osc, is_top=False, r_c_wall=r_c_wall)
left_wall_obj  = Wall([-half_w, -y_top0*2], [-half_w,  y_top0*2], normal=[+1, 0])
left_wall_obj.set_r_c_wall(r_c_wall)
right_wall_obj = Wall([ half_w, -y_top0*2], [ half_w,  y_top0*2], normal=[-1, 0])
right_wall_obj.set_r_c_wall(r_c_wall)

# Set wall render styles for make_movie()
for w in (top_wall_obj, bot_wall_obj, left_wall_obj, right_wall_obj):
    w.set_render(color='#666666', linewidth=2.5, alpha=0.9)

sys_api.add_object(left_wall_obj)
sys_api.add_object(right_wall_obj)
sys_api.add_object(top_wall_obj)
sys_api.add_object(bot_wall_obj)

p1 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0, -cy_disk))
p2 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0,  cy_disk))
A0 = p1.A0

sys_api.initialize_from_particles([p1, p2], alpha_damp=alpha_damp, dt=dt)

cc_gap_init = (np.linalg.norm(sys_api._state['x_cm'].numpy()[1]
                               - sys_api._state['x_cm'].numpy()[0])
               - 2.0*(R0 + r_c))
print(f"  Initial cc_gap = {cc_gap_init:.4e}")

# ── Machine-precision cross-check (10 steps API vs TF direct) ────────────────
def _prim_direct(t_now):
    y_t = y_top0 - A_wall * np.sin(omega_osc * t_now)
    y_b = y_bot0 + A_wall * np.sin(omega_osc * t_now)
    segs = [
        LineSegment([-half_w, -y_top0*2], [-half_w,  y_top0*2], [+1, 0]),
        LineSegment([ half_w, -y_top0*2], [ half_w,  y_top0*2], [-1, 0]),
        LineSegment([-half_w, y_t], [half_w, y_t], [0, -1]),
        LineSegment([-half_w, y_b], [half_w, y_b], [0, +1]),
    ]
    for sg in segs: sg.r_c = r_c_wall
    return make_prim_data([(sg, 1., np.zeros(2), sg.r_c, 0., np.zeros(2))
                           for sg in segs])

_p1c = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                       center=(0.0, -cy_disk))
_p2c = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                       center=(0.0,  cy_disk))
_s_ref, _pr_ref = make_state([_p1c, _p2c])
_cm_ref = CandidacyManager(P=2, N=N, R0=R0, E=128, skin=R0,
                            periodic=False, periodic_x=False, periodic_y=False,
                            Lx=20.0, Ly=20.0, R0_arr=np.array([R0, R0]))
_cm_ref.update(_s_ref['x_cm'].numpy(), _s_ref['theta'].numpy())
_dt_ref = tf.constant(NP_DTYPE(dt),         dtype=DTYPE)
_al_ref = tf.constant(NP_DTYPE(alpha_damp), dtype=DTYPE)
_g_ref  = tf.constant(NP_DTYPE(0.0),        dtype=DTYPE)

_sys_chk = System(Lx=20.0, Ly=20.0, periodic_x=False, periodic_y=False)
_sys_chk.add_object(OscillatingWall(y_top0, half_w, A_wall, omega_osc,
                                     is_top=True,  r_c_wall=r_c_wall))
_sys_chk.add_object(OscillatingWall(y_bot0, half_w, A_wall, omega_osc,
                                     is_top=False, r_c_wall=r_c_wall))
_lw = Wall([-half_w,-y_top0*2],[-half_w,y_top0*2],[+1,0]); _lw.set_r_c_wall(r_c_wall)
_rw = Wall([ half_w,-y_top0*2],[ half_w,y_top0*2],[-1,0]); _rw.set_r_c_wall(r_c_wall)
_sys_chk.add_object(_lw); _sys_chk.add_object(_rw)
_p1d = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                       center=(0.0, -cy_disk))
_p2d = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                       center=(0.0,  cy_disk))
_sys_chk.initialize_from_particles([_p1d, _p2d], alpha_damp=alpha_damp, dt=dt)

_t_chk = 0.0
for _ in range(10):
    _pd = _prim_direct(_t_chk)
    if _cm_ref.needs_update(_s_ref['x_cm'].numpy(), _s_ref['theta'].numpy()):
        _cm_ref.update(_s_ref['x_cm'].numpy(), _s_ref['theta'].numpy())
    _caps = tf.constant(_cm_ref.CapCandidates, dtype=tf.int32)
    _t_tf = tf.constant(NP_DTYPE(_t_chk), dtype=DTYPE)
    _s_ref, _ = step_full_tf(_s_ref, _caps, _dt_ref, _al_ref, _g_ref,
                              _pr_ref, t=_t_tf, prim_data=_pd)
    _sys_chk.step(1)
    _t_chk += dt

_max_diff = np.max(np.abs(_sys_chk._state['x_all'].numpy()
                          - _s_ref['x_all'].numpy()))
print(f"\n  Machine-precision check (10 steps): max|Dx_all| = {_max_diff:.2e}")
mp_ok = _max_diff < 1e-12
print(f"  [{'PASS' if mp_ok else 'FAIL'}] API == TF direct  (max|Dx| < 1e-12)")

# ── Squeeze callback ──────────────────────────────────────────────────────────
def _area_poly(x):
    xn = np.roll(x, -1, axis=0)
    return 0.5 * abs(np.sum(x[:, 0]*xn[:, 1] - xn[:, 0]*x[:, 1]))

def _vert_strain(x_nodes):
    Nn = len(x_nodes)
    return 1.0 - (x_nodes[Nn//4, 1] - x_nodes[3*Nn//4, 1]) / (2.0 * R0)

def squeeze_frame(sys_obj):
    """Called by sys.run() at each sample boundary. Uses only the public API."""
    snap    = sys_obj.snapshot()
    t_now   = snap['t']
    x1, x2 = snap['x_all'][0], snap['x_all'][1]
    xc      = snap['x_cm']
    eps1    = _vert_strain(x1)
    eps2    = _vert_strain(x2)
    eps_avg = 0.5 * (eps1 + eps2)
    ws      = top_wall_obj.wall_strain_at(t_now)
    cc_gap  = float(np.linalg.norm(xc[1]-xc[0])) - 2.0*(R0+r_c)
    t_norm  = t_now / T_osc
    return {
        't'          : t_now,
        'wall_strain': ws,
        'eps1'       : eps1, 'eps2': eps2,
        'cc_gap'     : cc_gap,
        'A1'         : _area_poly(x1), 'A2': _area_poly(x2), 'A0': A0,
        'r_c'        : r_c, 'r_c_wall': r_c_wall,
        'text'       : (f"t/T={t_norm:.2f}  ε_wall={ws:.3f}\n"
                        f"ε_p={eps_avg:.3f}  cc_gap={cc_gap:+.4f}"),
    }

# ── Main run ──────────────────────────────────────────────────────────────────
SAMPLE_EVERY = n_steps // (20 * n_cycles)   # ~20 frames per oscillation cycle

t_start = time.time()
sys_api.run(
    N            = n_steps,
    sample_every = SAMPLE_EVERY,
    callback     = squeeze_frame,
    record_initial = True,
)
print(f"\n  Recorded {len(sys_api.frames)} frames  ({time.time()-t_start:.1f}s total)")

# Aggregate plumbing counters from sys_api.diag
frames       = sys_api.callback_data
diag_total   = {
    'n_cand_checks' : sum(d['n_cand_checks']  for d in sys_api.diag),
    'n_cand_updates': sum(d['n_cand_updates'] for d in sys_api.diag),
}

# ── Physics diagnostics ───────────────────────────────────────────────────────
print("\n  -- Physics diagnostics --")
max_wall_strain = max(f['wall_strain'] for f in frames)
max_eps         = max(0.5*(f['eps1']+f['eps2']) for f in frames)
max_dA          = max(abs(1.0 - 0.5*(f['A1']+f['A2'])/f['A0']) for f in frames)
min_cc_gap      = min(f['cc_gap'] for f in frames)

print(f"  peak wall_strain : {max_wall_strain:.4f}  (target ~{A_wall:.4f})")
print(f"  peak eps_particle: {max_eps:.4f}")
print(f"  peak |1-A/A0|    : {max_dA:.6f}  (area conservation)")
print(f"  min cc_gap       : {min_cc_gap:.4e}")

ok_strain  = abs(max_wall_strain - A_wall) < 0.005
ok_area    = max_dA < 0.05
ok_contact = min_cc_gap < 0.0

print(f"\n  [{'PASS' if ok_strain  else 'FAIL'}] wall_strain ~ A_wall  "
      f"({max_wall_strain:.4f} vs {A_wall})")
print(f"  [{'PASS' if ok_area   else 'FAIL'}] area conservation < 5%  "
      f"(actual {max_dA:.4f})")
print(f"  [{'PASS' if ok_contact else 'FAIL'}] disks make contact  "
      f"(min_cc_gap = {min_cc_gap:.4e})")

# ── Plumbing diagnostics ──────────────────────────────────────────────────────
print("\n  -- Plumbing diagnostics --")
n_retraces      = step_full_tf.experimental_get_tracing_count()
n_chunks        = n_steps // SAMPLE_EVERY
expected_checks = n_steps // 10   # cand_check_interval default = 10

print(f"  step_full_tf retraces : {n_retraces}        (expected: 1)")
print(f"  prim_data rebuilds    : 0        (static — oscillation encoded)")
print(f"  cand py_func calls    : {diag_total['n_cand_checks']:<6}   (expected: ~{expected_checks})")
print(f"  cand C++ updates      : {diag_total['n_cand_updates']:<6}")
print(f"  .numpy() calls        : {n_chunks:<6}   = N_CHUNKS")

plumbing_ok = (n_retraces <= 2 and
               diag_total['n_cand_checks'] <= expected_checks + n_chunks)

print(f"\n  [{'PASS' if n_retraces <= 2 else 'FAIL'}] retraces <= 2")
print(f"  [{'PASS' if plumbing_ok else 'FAIL'}] cand checks ~ n_steps/interval")

# ── Render movie ──────────────────────────────────────────────────────────────
sys_api.make_movie(
    args.out, fps=10, n_arc=6,
    title=f'Two-disk oscillatory squeeze  (nu={nu}, 10% strain)  [API tf-fast]',
)

all_pass = ok_strain and ok_area and ok_contact and mp_ok and plumbing_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test A (elastic, API tf-fast)")
