"""
test_40_patches.py — Phase 4.0 pre-requisite patch verification.

Patches tested:
  A: CandidacyManager per-pair R0 threshold (polydisperse / mixed-type)
  B: step_rb_tf shape_frozen (elastic DOFs zeroed for frozen particles)
  C: absolute time t threaded through step_full_tf → step_rb_tf
  D: CandidacyManager per-particle skin_arr (mixed-type)

Run: python src/epd/tests/test_40_patches.py
Produces: results/phase40_patches/patch_tests.png
All tests print PASS/FAIL lines.  Exit code 0 = all pass, 1 = any fail.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── path setup ────────────────────────────────────────────────────────────────
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)

os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')
import tensorflow as tf

import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.simulation.candidacy_manager import CandidacyManager
from src.simulation.capsule_shell import CapsuleParticle
from src.simulation.tf_sim import make_state, step_full_tf, make_prim_data, DTYPE, NP_DTYPE

OUTDIR = os.path.join(ROOT, 'results', 'phase40_patches')
os.makedirs(OUTDIR, exist_ok=True)

results = []   # list of (name, passed, detail)

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════════════
# Patch A — per-pair R0 threshold
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Patch A: per-pair R0 threshold ──────────────────────────────────────")

N_nodes = 32
skin    = 0.3

# A1: backward compat — uniform R0, same threshold as old code
mgr_old = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=skin)
mgr_new = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=skin, R0_arr=None)
# Both should give threshold = 2.0 + 0.3 = 2.3
check("A1 threshold compat (uniform R0)",
      abs(mgr_old._pair_threshold(0, 1) - (2.0 + skin)) < 1e-12,
      f"thresh={mgr_old._pair_threshold(0, 1):.4f}")

# A2: heterogeneous R0 — small particle near large particle should be detected
# R0=0.5 + R0=1.5: threshold = 0.5 + 1.5 + skin = 2.3
# R0=1.0 + R0=1.0: threshold = 1.0 + 1.0 + skin = 2.3 (same here)
R0_small, R0_large = 0.5, 1.5
mgr_hetero = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=skin,
                               R0_arr=[R0_small, R0_large])
expected = R0_small + R0_large + skin
check("A2 heterogeneous R0 threshold",
      abs(mgr_hetero._pair_threshold(0, 1) - expected) < 1e-12,
      f"expected={expected:.3f} got={mgr_hetero._pair_threshold(0,1):.3f}")

# A3: R0=1.0 + R0=2.0 — new threshold 3.0+skin; old would be 2.0+skin = miss
R0_a, R0_b = 1.0, 2.0
mgr_mix = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=32, skin=skin,
                            R0_arr=[R0_a, R0_b])
# Place particles so CM separation = 3.1 (beyond old threshold 2.3, inside new 3.3)
x_cm = np.array([[0.0, 0.0], [3.1, 0.0]])
theta = np.array([0.0, 0.0])
mgr_mix.update(x_cm, theta)
# Should have found some candidates (pair is within range)
found_any = np.any(mgr_mix.CapCandidates != 0)
check("A3 extended range detects R0=(1,2) at sep=3.1", found_any,
      f"any_candidates={found_any}")

# Old uniform manager should miss this pair (threshold=2.3 < sep=3.1)
mgr_old2 = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=32, skin=skin)
mgr_old2.update(x_cm, theta)
old_miss = not np.any(mgr_old2.CapCandidates != 0)
check("A3 old code misses R0=(1,2) at sep=3.1 (expected miss)", old_miss,
      f"old_any={np.any(mgr_old2.CapCandidates != 0)}")


# ══════════════════════════════════════════════════════════════════════════════
# Patch D — per-particle skin_arr
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Patch D: per-particle skin_arr ──────────────────────────────────────")

# D1: skin_arr uniform — same threshold as scalar skin
skin1, skin2 = 0.3, 0.3
mgr_d1 = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=0.3,
                           skin_arr=[skin1, skin2])
check("D1 skin_arr uniform = scalar skin",
      abs(mgr_d1._pair_threshold(0, 1) - (2.0 + 0.3)) < 1e-12,
      f"thresh={mgr_d1._pair_threshold(0,1):.4f}")

# D2: different skin per type — emulsion (larger skin) vs elastic (smaller skin)
skin_elastic, skin_emulsion = 0.2, 0.6
mgr_d2 = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=0.3,
                           skin_arr=[skin_elastic, skin_emulsion])
expected_d2 = 1.0 + 1.0 + 0.5*(skin_elastic + skin_emulsion)
check("D2 mixed-type skin threshold",
      abs(mgr_d2._pair_threshold(0, 1) - expected_d2) < 1e-12,
      f"expected={expected_d2:.4f} got={mgr_d2._pair_threshold(0,1):.4f}")

# D3: backward compat — skin_arr=None → falls back to scalar skin
mgr_d3 = CandidacyManager(P=2, N=N_nodes, R0=1.0, E=16, skin=0.4)
check("D3 skin_arr=None backward compat",
      abs(mgr_d3._pair_threshold(0, 1) - (2.0 + 0.4)) < 1e-12,
      f"thresh={mgr_d3._pair_threshold(0,1):.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Patch B — shape_frozen freezes elastic DOFs
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Patch B: shape_frozen in step_rb_tf ─────────────────────────────────")

def _make_two_particles(N=16, sep=2.5, tau_b=0.2, q=1.0, S=1.0):
    TAU    = np.sqrt(12.0 * tau_b)
    K_area = q * (12.0 * S / TAU**2)
    C      = 3000.0 * S * (1.0 + q)
    p0 = CapsuleParticle(N=N, R0=1.0, tau=TAU/np.sqrt(12), S=S, C=C, K_area=K_area,
                         center=(0.0, 0.0))
    p1 = CapsuleParticle(N=N, R0=1.0, tau=TAU/np.sqrt(12), S=S, C=C, K_area=K_area,
                         center=(sep, 0.0))
    return [p0, p1]

particles_b = _make_two_particles(N=16, sep=5.0)  # far apart, no contact forces
state_b, params_b = make_state(particles_b)

# Give particle 0 a nonzero initial elastic displacement
u_init = np.zeros((2, 16, 2), dtype=NP_DTYPE)
u_init[0] = 0.05 * np.random.default_rng(42).standard_normal((16, 2))
u_dot_init = np.zeros_like(u_init)
u_dot_init[0] = 0.01 * np.random.default_rng(7).standard_normal((16, 2))
state_b = {**state_b,
           'u':     tf.constant(u_init,     dtype=DTYPE),
           'u_dot': tf.constant(u_dot_init, dtype=DTYPE)}

# Freeze particle 0 (index 0), free particle 1
frozen = np.array([1.0, 0.0], dtype=NP_DTYPE)   # 1 = frozen, 0 = free
import copy
params_frozen = dict(params_b)
params_frozen['shape_frozen'] = tf.constant(frozen, dtype=DTYPE)

dt         = tf.constant(1e-4, dtype=DTYPE)
alpha_damp = tf.constant(2.0,  dtype=DTYPE)
g          = tf.constant(0.0,  dtype=DTYPE)
prim_data  = make_prim_data([])

# Candidacy (no contacts)
mgr_b = CandidacyManager(P=2, N=16, R0=1.0, E=16, skin=0.3)
x_cm0 = np.array([[0.0, 0.0], [5.0, 0.0]])
theta0 = np.zeros(2)
mgr_b.update(x_cm0, theta0)
CapCand = tf.constant(mgr_b.CapCandidates, dtype=tf.int32)

t0 = tf.constant(0.0, dtype=DTYPE)
new_state_b, _ = step_full_tf(state_b, CapCand, dt, alpha_damp, g,
                               params_frozen, t=t0, prim_data=prim_data)

u_new   = new_state_b['u'].numpy()
u_dot_new = new_state_b['u_dot'].numpy()

# Particle 0 (frozen): u and u_dot must be exactly zero
b1 = check("B1 frozen particle u=0",
           np.allclose(u_new[0], 0.0, atol=1e-15),
           f"max|u[0]|={np.max(np.abs(u_new[0])):.2e}")
b2 = check("B2 frozen particle u_dot=0",
           np.allclose(u_dot_new[0], 0.0, atol=1e-15),
           f"max|u_dot[0]|={np.max(np.abs(u_dot_new[0])):.2e}")

# Particle 1 (free): u and u_dot should be nonzero (internal forces act)
b3 = check("B3 free particle u evolves",
           np.max(np.abs(u_new[1])) > 0.0,
           f"max|u[1]|={np.max(np.abs(u_new[1])):.2e}")

# Run without freeze — particles must differ
params_free = dict(params_b)
params_free['shape_frozen'] = tf.constant([0.0, 0.0], dtype=DTYPE)
new_state_free, _ = step_full_tf(state_b, CapCand, dt, alpha_damp, g,
                                  params_free, t=t0, prim_data=prim_data)
u_free = new_state_free['u'].numpy()
check("B4 frozen vs free gives different u[0]",
      not np.allclose(u_new[0], u_free[0], atol=1e-10),
      f"diff_u0={np.max(np.abs(u_new[0]-u_free[0])):.2e}")


# ══════════════════════════════════════════════════════════════════════════════
# Patch C — absolute time t threads through to driven-particle interpolation
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Patch C: absolute time t threading ──────────────────────────────────")

from src.simulation.tf_sim import set_driven, make_traj

particles_c = _make_two_particles(N=16, sep=5.0)
state_c, params_c = make_state(particles_c)

# Drive particle 0 at v_x = 1.0 (DC), no AC
traj_row = make_traj(v_dc=(1.0, 0.0))
set_driven(params_c, [0], [traj_row])

prim_data_c = make_prim_data([])
CapCand_c   = tf.constant(mgr_b.CapCandidates, dtype=tf.int32)  # same zero candidacy
dt_c        = tf.constant(1e-3, dtype=DTYPE)

# Step at t=0 and at t=10 with DC-only motion — CM velocity should be identical
t_early = tf.constant(0.0,  dtype=DTYPE)
t_late  = tf.constant(10.0, dtype=DTYPE)

new_e, _ = step_full_tf(state_c, CapCand_c, dt_c, alpha_damp, g,
                         params_c, t=t_early, prim_data=prim_data_c)
new_l, _ = step_full_tf(state_c, CapCand_c, dt_c, alpha_damp, g,
                         params_c, t=t_late,  prim_data=prim_data_c)

# Driven particle 0: v_cm_x should equal 1.0 after one step (DC, t-independent)
check("C1 DC driven particle v_x=1 at t=0",
      abs(float(new_e['v_cm'][0, 0]) - 1.0) < 1e-10,
      f"v_x[0]={float(new_e['v_cm'][0, 0]):.6f}")
check("C2 DC driven particle v_x=1 at t=10 (t-independent)",
      abs(float(new_l['v_cm'][0, 0]) - 1.0) < 1e-10,
      f"v_x[0]={float(new_l['v_cm'][0, 0]):.6f}")

# AC motion: omega=1, amplitude=1 → v_x(t) = cos(t)
# v_x at t=0 should be 1.0; at t=π/2 should be ~0
traj_ac = make_traj(v_dc=(0.0, 0.0), v_ac=(1.0, 0.0), freq=(1.0, 0.0))
set_driven(params_c, [0], [traj_ac])

t_zero   = tf.constant(0.0,            dtype=DTYPE)
t_half_pi = tf.constant(np.pi / 2.0,   dtype=DTYPE)

new_0, _ = step_full_tf(state_c, CapCand_c, dt_c, alpha_damp, g,
                         params_c, t=t_zero,    prim_data=prim_data_c)
new_pi2, _ = step_full_tf(state_c, CapCand_c, dt_c, alpha_damp, g,
                           params_c, t=t_half_pi, prim_data=prim_data_c)

v0   = float(new_0['v_cm'][0, 0])
v_pi2 = float(new_pi2['v_cm'][0, 0])
check("C3 AC driven particle v_x=cos(0)=1 at t=0",
      abs(v0 - 1.0) < 1e-10,
      f"v_x={v0:.6f}")
check("C4 AC driven particle v_x=cos(π/2)≈0 at t=π/2",
      abs(v_pi2) < 1e-8,
      f"v_x={v_pi2:.2e}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary figure
# ══════════════════════════════════════════════════════════════════════════════
n_pass = sum(r[1] for r in results)
n_fail = sum(not r[1] for r in results)

fig, ax = plt.subplots(figsize=(8, len(results) * 0.4 + 1.5))
for i, (name, passed, detail) in enumerate(results):
    color = '#2ecc71' if passed else '#e74c3c'
    ax.barh(i, 1, color=color, height=0.7)
    label = f"{'PASS' if passed else 'FAIL'} — {name}"
    if detail:
        label += f"  ({detail})"
    ax.text(0.02, i, label, va='center', fontsize=9)
ax.set_xlim(0, 1)
ax.axis('off')
ax.set_title(f"Phase 4.0 Patch Tests — {n_pass} PASS / {n_fail} FAIL",
             fontsize=12, fontweight='bold')
plt.tight_layout()
outpath = os.path.join(OUTDIR, 'patch_tests.png')
plt.savefig(outpath, dpi=120)
print(f"\nSummary plot: {outpath}")
print(f"\n{'='*60}")
print(f"TOTAL: {n_pass}/{len(results)} PASS")
if n_fail:
    print("FAILED tests:")
    for name, passed, detail in results:
        if not passed:
            print(f"  - {name}  ({detail})")

sys.exit(0 if n_fail == 0 else 1)
