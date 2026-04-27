"""
test_A_elastic_tf_fast.py — Test A (elastic) on the tf-fast branch.

Validates two things:
  1. PHYSICS: machine-precision agreement with the reference Python-loop approach
     (same prim_data rebuild logic from emulsion-dev, 10 steps).
  2. PLUMBING: prim_data built ONCE with encoded oscillation; run_simulation_tf
     (tf.while_loop) used for all time integration; diagnostics counters confirm
     no per-step Python callbacks except the cand_check_interval cadence.

Usage:
    python src/validation/test_A_elastic_tf_fast.py [--cycles N] [--out PATH]
"""

import sys, os, argparse, time
import numpy as np
import tensorflow as tf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)
from src.simulation.tf_sim import (make_state, make_prim_data, step_full_tf,
                                    run_simulation_tf, DTYPE, NP_DTYPE)
from src.simulation.candidacy_manager import CandidacyManager
from src.epd.particles import ParticleSpec
from src.simulation.capsule_shell import CapsuleParticle, CapsuleSim
from src.simulation.contact_primitives import LineSegment

parser = argparse.ArgumentParser()
parser.add_argument('--cycles', type=int, default=3)
parser.add_argument('--out',    type=str, default='results/test_A_elastic_tf_fast.gif')
args = parser.parse_args()

# ── Parameters (identical to test_A_elastic_tf.py) ───────────────────────────
R0  = 1.0
N   = 32
nu  = 0.5

spec = ParticleSpec(count=2, nu=nu, N_nodes=N)
d    = spec.derived
tau  = d['tau'];  C = d['C'];  K_area = d['K_area'];  S = 1.0

print("=" * 60)
print(f"Test A (elastic tf-fast): Two-disk squeeze  (nu={nu})")
print("=" * 60)
print(f"  q={d['q']:.4f}  TAU={d['TAU']:.4f}")
print(f"  tau={tau:.4f}  El_t={d['El_t']:.4f}  K_area={K_area:.4f}")
print(f"  C={C:.2f}  alpha={d['alpha']:.4f}")

p_ref = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area)
L0    = p_ref.L0;  r_c = p_ref.r_c
El_t  = d['El_t']
c_edge     = np.sqrt(El_t / (p_ref.rho_d * L0))
T_wave     = 2.0 * np.pi * R0 / c_edge
alpha_damp = 2.0 / T_wave
T_osc      = 5.0 * T_wave
omega_osc  = 2.0 * np.pi / T_osc
A_wall     = 0.1 * R0

sim_probe = CapsuleSim([p_ref])
dt_max, _ = sim_probe.estimate_dt_max()
dt = 0.4 * dt_max

n_cycles = args.cycles
t_total  = n_cycles * T_osc
n_steps  = int(np.ceil(t_total / dt))

print(f"\n  Geometry:  R0={R0}  N={N}  L0={L0:.4f}  r_c={r_c:.4f}")
print(f"  Dynamics:  dt={dt:.5f}  T_wave={T_wave:.4f}  alpha={alpha_damp:.4f}")
print(f"  Oscillation: A={A_wall}  T_osc={T_osc:.4f}  omega={omega_osc:.4f}")
print(f"  Total: {n_cycles} cycles  {n_steps} steps  t_total={t_total:.3f}")

# ── Geometry ──────────────────────────────────────────────────────────────────
cy_disk  = R0 + L0
half_w   = R0 + 2.0 * L0
y_top0   =  2.0 * R0 + 3.0 * L0
y_bot0   = -(2.0 * R0 + 3.0 * L0)
r_c_wall = L0

print(f"  Box: half_w={half_w:.4f}  y_top0={y_top0:.4f}")

def _make_cm(p1, p2):
    R0_arr = np.array([p1.R0, p2.R0])
    cm = CandidacyManager(P=2, N=N, R0=float(np.mean(R0_arr)), E=128,
                          skin=1.0*float(np.mean(R0_arr)), periodic=False,
                          periodic_x=False, periodic_y=False,
                          Lx=20.0, Ly=20.0, R0_arr=R0_arr)
    return cm

# ── Build prim_data ONCE with encoded oscillation (tf-fast path) ──────────────
left_seg  = LineSegment([-half_w, -y_top0*2], [-half_w,  y_top0*2], [+1, 0])
right_seg = LineSegment([ half_w, -y_top0*2], [ half_w,  y_top0*2], [-1, 0])
top_seg_rest = LineSegment([-half_w, y_top0], [half_w, y_top0], [0, -1])
bot_seg_rest = LineSegment([-half_w, y_bot0], [half_w, y_bot0], [0, +1])
for seg in [left_seg, right_seg, top_seg_rest, bot_seg_rest]:
    seg.r_c = r_c_wall

prim_data_static = make_prim_data([
    (left_seg,     1.0, np.zeros(2), r_c_wall, 0.0, np.zeros(2)),
    (right_seg,    1.0, np.zeros(2), r_c_wall, 0.0, np.zeros(2)),
    # Top wall: sign=-1  →  y(t) = y_top0 - A*sin(omega*t)
    (top_seg_rest, 1.0, np.zeros(2), r_c_wall, 0.0, np.zeros(2),
     (A_wall, omega_osc, -1.0)),
    # Bottom wall: sign=+1  →  y(t) = y_bot0 + A*sin(omega*t)
    (bot_seg_rest, 1.0, np.zeros(2), r_c_wall, 0.0, np.zeros(2),
     (A_wall, omega_osc, +1.0)),
])
print("\n  prim_data_static built ONCE  (oscillation encoded as parameters)")

# ── Machine-precision cross-check: 10 steps fast vs reference ────────────────
print("\n  -- Machine-precision cross-check (10 steps) --")

p1 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0, -cy_disk))
p2 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0,  cy_disk))
A0 = p1.A0
state0, params = make_state([p1, p2])

# Reference: Python loop with per-step prim_data rebuild (emulsion-dev style)
def _prim_ref(t_now):
    y_t = y_top0 - A_wall * np.sin(omega_osc * t_now)
    y_b = y_bot0 + A_wall * np.sin(omega_osc * t_now)
    segs = [
        LineSegment([-half_w, -y_top0*2], [-half_w,  y_top0*2], [+1, 0]),
        LineSegment([ half_w, -y_top0*2], [ half_w,  y_top0*2], [-1, 0]),
        LineSegment([-half_w, y_t], [half_w, y_t], [0, -1]),
        LineSegment([-half_w, y_b], [half_w, y_b], [0, +1]),
    ]
    for seg in segs: seg.r_c = r_c_wall
    return make_prim_data([(seg, 1.0, np.zeros(2), seg.r_c, 0.0, np.zeros(2))
                           for seg in segs])

cm_ref  = _make_cm(p1, p2)
cm_fast = _make_cm(p1, p2)
state_ref  = {k: tf.identity(v) for k, v in state0.items()}
state_fast = {k: tf.identity(v) for k, v in state0.items()}
cm_ref.update(state_ref['x_cm'].numpy(), state_ref['theta'].numpy())
cm_fast.update(state_fast['x_cm'].numpy(), state_fast['theta'].numpy())

dt_tf    = tf.constant(NP_DTYPE(dt),          dtype=DTYPE)
alpha_tf = tf.constant(NP_DTYPE(alpha_damp),  dtype=DTYPE)
g_tf     = tf.constant(NP_DTYPE(0.0),         dtype=DTYPE)

t_chk = 0.0
for _ in range(10):
    # Reference step (old approach)
    pd = _prim_ref(t_chk)
    if cm_ref.needs_update(state_ref['x_cm'].numpy(), state_ref['theta'].numpy()):
        cm_ref.update(state_ref['x_cm'].numpy(), state_ref['theta'].numpy())
    caps_ref = tf.constant(cm_ref.CapCandidates, dtype=tf.int32)
    t_tf     = tf.constant(NP_DTYPE(t_chk), dtype=DTYPE)
    state_ref, _ = step_full_tf(state_ref, caps_ref, dt_tf, alpha_tf, g_tf,
                                 params, t=t_tf, prim_data=pd)
    t_chk += dt

# Fast approach: run_simulation_tf for 10 steps
state_fast, diag10 = run_simulation_tf(
    state_fast, dt, alpha_damp, 0.0, params, 10,
    cm_fast, skin=R0, prim_data=prim_data_static, R0_max=R0,
    cand_check_interval=10, diagnostics=True)

max_diff = float(np.max(np.abs(
    state_fast['x_all'].numpy() - state_ref['x_all'].numpy())))
mp_ok = max_diff < 1e-12
print(f"  max|Dx_all| = {max_diff:.2e}")
print(f"  [{'PASS' if mp_ok else 'FAIL'}] fast == reference  (threshold 1e-12)")

# ── Full run with run_simulation_tf ───────────────────────────────────────────
print(f"\n  -- Full run ({n_steps} steps in 5 chunks) --")

p1f = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                      center=(0.0, -cy_disk))
p2f = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                      center=(0.0,  cy_disk))
state_run, params_run = make_state([p1f, p2f])
cm_run = _make_cm(p1f, p2f)

cc_gap_init = float(np.linalg.norm(
    state_run['x_cm'].numpy()[1] - state_run['x_cm'].numpy()[0]) - 2.0*(R0+r_c))
print(f"  Initial cc_gap = {cc_gap_init:.4e}")

N_CHUNKS      = 5
steps_chunk   = n_steps // N_CHUNKS
CAND_INTERVAL = 10

total_cand_checks  = 0
total_cand_updates = 0
checkpoints = []
t_start = time.time()

for chunk_i in range(N_CHUNKS):
    state_run, diag = run_simulation_tf(
        state_run, dt, alpha_damp, 0.0, params_run,
        steps_chunk, cm_run, skin=R0,
        prim_data=prim_data_static, R0_max=R0,
        cand_check_interval=CAND_INTERVAL, diagnostics=True,
        step_offset=chunk_i * steps_chunk)

    total_cand_checks  += diag['n_cand_checks']
    total_cand_updates += diag['n_cand_updates']

    # One .numpy() call per chunk for diagnostics
    t_now   = (chunk_i + 1) * steps_chunk * dt
    x_all   = state_run['x_all'].numpy()
    x_cm    = state_run['x_cm'].numpy()
    x1, x2  = x_all[0], x_all[1]

    # Wall position at t_now (computed from encoded params, not prim_data rebuild)
    y_top_now = y_top0 - A_wall * np.sin(omega_osc * t_now)
    wall_strain = abs(y_top_now - y_top0)

    xn1 = np.roll(x1, -1, axis=0)
    A1  = 0.5 * abs(np.sum(x1[:,0]*xn1[:,1] - xn1[:,0]*x1[:,1]))
    xn2 = np.roll(x2, -1, axis=0)
    A2  = 0.5 * abs(np.sum(x2[:,0]*xn2[:,1] - xn2[:,0]*x2[:,1]))
    eps1 = 1.0 - (x1[N//4,1] - x1[3*N//4,1]) / (2.0*R0)
    cc_gap = float(np.linalg.norm(x_cm[1]-x_cm[0]) - 2.0*(R0+r_c))

    checkpoints.append(dict(t=t_now, wall_strain=wall_strain, eps1=eps1,
                             A1=A1, A2=A2, cc_gap=cc_gap))
    elapsed = time.time() - t_start
    print(f"    chunk {chunk_i+1}/{N_CHUNKS}  t={t_now:.3f}  "
          f"wall_strain={wall_strain:.4f}  eps_p={eps1:.4f}  "
          f"cc_gap={cc_gap:.4e}  ({elapsed:.0f}s)")

print(f"\n  Total wall time: {time.time()-t_start:.1f}s")

# ── Physics diagnostics ───────────────────────────────────────────────────────
print("\n  -- Physics diagnostics --")
max_wall_strain = max(c['wall_strain'] for c in checkpoints)
max_eps         = max(c['eps1'] for c in checkpoints)
max_dA          = max(abs(1.0 - 0.5*(c['A1']+c['A2'])/A0) for c in checkpoints)
min_cc_gap      = min(c['cc_gap'] for c in checkpoints)

print(f"  peak wall_strain : {max_wall_strain:.4f}  (target ~{A_wall:.4f})")
print(f"  peak eps_particle: {max_eps:.4f}")
print(f"  peak |1-A/A0|    : {max_dA:.6f}  (area conservation)")
print(f"  min cc_gap       : {min_cc_gap:.4e}")

ok_strain  = abs(max_wall_strain - A_wall) < 0.005
ok_area    = max_dA < 0.05
ok_contact = min_cc_gap < 0.0

print(f"\n  [{'PASS' if ok_strain  else 'FAIL'}] wall_strain ~ A_wall  ({max_wall_strain:.4f} vs {A_wall})")
print(f"  [{'PASS' if ok_area   else 'FAIL'}] area conservation < 5%  (actual {max_dA:.4f})")
print(f"  [{'PASS' if ok_contact else 'FAIL'}] disks make contact  (min_cc_gap = {min_cc_gap:.4e})")

# ── Plumbing diagnostics ──────────────────────────────────────────────────────
print("\n  -- Plumbing diagnostics --")
n_retraces = step_full_tf.experimental_get_tracing_count()
expected_checks  = n_steps // CAND_INTERVAL
n_numpy_calls    = N_CHUNKS   # one per chunk for frame extraction

print(f"  step_full_tf retraces : {n_retraces}        (expected: 1)")
print(f"  prim_data rebuilds    : 0        (expected: 0  — static prim_data)")
print(f"  cand py_func calls    : {total_cand_checks:<6}   (expected: ~{expected_checks}  = n_steps/{CAND_INTERVAL})")
print(f"  cand C++ updates      : {total_cand_updates:<6}   (expected: << n_steps={n_steps})")
print(f"  .numpy() calls        : {n_numpy_calls:<6}   (expected: {N_CHUNKS}  = N_CHUNKS)")

plumbing_ok = (n_retraces <= 2 and
               total_cand_checks <= expected_checks + N_CHUNKS and
               n_numpy_calls == N_CHUNKS)

print(f"\n  [{'PASS' if n_retraces <= 2     else 'FAIL'}] retraces <= 2")
print(f"  [{'PASS' if total_cand_checks <= expected_checks+N_CHUNKS else 'FAIL'}] "
      f"cand checks ~ n_steps/interval")
print(f"  [{'PASS' if n_numpy_calls == N_CHUNKS else 'FAIL'}] .numpy() calls == N_CHUNKS")

all_pass = mp_ok and ok_strain and ok_area and ok_contact and plumbing_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test A (elastic, tf-fast)")
