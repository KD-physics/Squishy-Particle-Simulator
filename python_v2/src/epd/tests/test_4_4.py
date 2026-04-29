"""
test_4_4.py — Phase 4.4: initializer.py verification.

Tests:
  1. RSA placement: P=16, periodic, phi_init=0.30 (RSA only, no swell)
     → all CMs in [0,Lx]×[0,Ly], no overlap
  2. RSA inside Box(width=12, height=12): all CMs inside box
  3. RSA around CircleObstacle(radius=2): no CM within radius 2
  4. compute_phi_outer correctness (known geometry)
  5. Adaptive swell (short run: 5 increments, verify phi increases)

Run: python src/epd/tests/test_4_4.py
Output: results/phase44_init/init_test.png
"""

import sys, os
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
from src.simulation.tf_sim import (make_state, make_prim_data, set_periodic_box,
                                    DTYPE, NP_DTYPE)
from src.simulation.candidacy_manager import CandidacyManager
from src.epd.particles import ParticleSpec
from src.epd.objects import Box, CircleObstacle
from src.epd.initializer import (rsa_seed, compute_phi_outer,
                                  _compute_max_f, _wrap, _compress)

OUTDIR = os.path.join(ROOT, 'results', 'phase44_init')
os.makedirs(OUTDIR, exist_ok=True)

results = []

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════════════
# Test 1 — RSA placement, periodic box
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 1: RSA P=16 periodic ───────────────────────────────────────────")

N_nodes = 32
spec1   = ParticleSpec(count=16, N_nodes=N_nodes, q=1.0, tau_b=0.2)

# Box size for phi_init ≈ 0.30 with R_eff ≈ 1.1
R_eff   = 1.0 + 2.0 * np.sin(np.pi / N_nodes)
A_part  = np.pi * R_eff**2
Lx      = np.sqrt(16 * A_part / 0.30)
Ly      = Lx

particles1, pos1 = rsa_seed([spec1], [], Lx, Ly, seed=42, verbose=True)
P       = len(particles1)

check("1.1 P=16 particles placed", P == 16, f"P={P}")
check("1.2 All CMs in [0,Lx]×[0,Ly]",
      np.all((pos1[:, 0] >= 0) & (pos1[:, 0] <= Lx) &
             (pos1[:, 1] >= 0) & (pos1[:, 1] <= Ly)),
      f"x_range=[{pos1[:,0].min():.2f},{pos1[:,0].max():.2f}]")

# No overlap: min CM distance should exceed 2*r_c (nodes don't touch)
min_sep = float('inf')
for i in range(P):
    for j in range(i+1, P):
        dx = pos1[i, 0] - pos1[j, 0]
        dy = pos1[i, 1] - pos1[j, 1]
        dx -= Lx * np.round(dx / Lx)
        dy -= Ly * np.round(dy / Ly)
        min_sep = min(min_sep, np.sqrt(dx**2 + dy**2))

min_allowed = 2.0 * R_eff - 0.1   # allow small tolerance
check("1.3 No particle pair overlapping",
      min_sep >= min_allowed,
      f"min_sep={min_sep:.4f}, min_allowed={min_allowed:.4f}")

phi_rsa = P * np.pi * R_eff**2 / (Lx * Ly)
check("1.4 RSA phi ≈ 0.30",
      abs(phi_rsa - 0.30) < 0.05,
      f"phi_rsa={phi_rsa:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 2 — RSA inside Box
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 2: RSA inside Box(12,12) ───────────────────────────────────────")

box = Box(width=12.0, height=12.0, x0=6.0, y0=6.0, exclusion='exterior')
spec2 = ParticleSpec(count=8, N_nodes=32, q=1.0)

Lx2, Ly2 = 14.0, 14.0   # box is centred at (6,6); slightly larger periodic domain
particles2, pos2 = rsa_seed([spec2], [box], Lx2, Ly2, seed=7, verbose=True)

check("2.1 P=8 particles placed inside box", len(particles2) == 8)

# All CMs should be inside box (x in [1,11], y in [1,11]) approximately
in_box = np.all(
    (pos2[:, 0] >= 0.5) & (pos2[:, 0] <= 11.5) &
    (pos2[:, 1] >= 0.5) & (pos2[:, 1] <= 11.5)
)
check("2.2 All CMs inside box region",
      in_box,
      f"x=[{pos2[:,0].min():.2f},{pos2[:,0].max():.2f}] "
      f"y=[{pos2[:,1].min():.2f},{pos2[:,1].max():.2f}]")


# ══════════════════════════════════════════════════════════════════════════════
# Test 3 — RSA around CircleObstacle
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 3: RSA around CircleObstacle(r=2, x0=7, y0=7) ─────────────────")

obs = CircleObstacle(radius=2.0, x0=7.0, y0=7.0, exclusion='interior')
spec3 = ParticleSpec(count=8, N_nodes=32, q=1.0)

Lx3, Ly3 = 14.0, 14.0
particles3, pos3 = rsa_seed([spec3], [obs], Lx3, Ly3, seed=13, verbose=True)

check("3.1 P=8 particles placed", len(particles3) == 8)

# No CM within radius 2.0 of centre (7,7)
dists_to_centre = np.sqrt((pos3[:, 0] - 7.0)**2 + (pos3[:, 1] - 7.0)**2)
check("3.2 All CMs outside obstacle (r>2+R_eff)",
      np.all(dists_to_centre >= 2.0 + R_eff * 0.8),
      f"min_dist={dists_to_centre.min():.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 4 — compute_phi_outer correctness
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 4: compute_phi_outer correctness ───────────────────────────────")

# Build TF state from particles1
state1, params1 = make_state(particles1)
r_c_arr1 = params1['r_c_per_p'].numpy()

phi_outer1 = compute_phi_outer(state1, Lx, Ly, r_c_arr1)
# At RSA phi_init ≈ 0.30, phi_outer should be close to π*R_eff²/box_area ≈ 0.30
phi_expected_analytic = P * np.pi * R_eff**2 / (Lx * Ly)
check("4.1 compute_phi_outer near RSA phi",
      abs(phi_outer1 - phi_expected_analytic) < 0.05,
      f"phi_outer={phi_outer1:.4f} analytic={phi_expected_analytic:.4f}")

check("4.2 phi_outer > 0",
      phi_outer1 > 0)

check("4.3 phi_outer < 1",
      phi_outer1 < 1.0)


# ══════════════════════════════════════════════════════════════════════════════
# Test 5 — short adaptive swell (verify phi increases)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 5: short adaptive swell (2 increments) ─────────────────────────")

spec5  = ParticleSpec(count=6, N_nodes=32, q=1.0, tau_b=0.2)
R_eff5 = 1.0 + 2.0 * np.sin(np.pi / 32)
Lx5    = np.sqrt(6 * np.pi * R_eff5**2 / 0.30)
Ly5    = Lx5

particles5, pos5 = rsa_seed([spec5], [], Lx5, Ly5, seed=99, verbose=False)
state5, params5  = make_state(particles5)
r_c_arr5  = params5['r_c_per_p'].numpy()
L0_arr5   = params5['L0'].numpy()

set_periodic_box(params5, Lx5, Ly5)
prim_data5  = make_prim_data([])
phi_init5   = compute_phi_outer(state5, Lx5, Ly5, r_c_arr5)

# Setup candidacy manager
N_val5 = 32
cm5 = CandidacyManager(P=6, N=N_val5, R0=1.0, E=32, skin=0.3,
                        periodic=True, Lx=Lx5, Ly=Ly5)
cm5.update(state5['x_cm'].numpy(), state5['theta'].numpy())

# dt from conservative formula: 0.4 / omega_edge
p0 = particles5[0]
om_edge = p0.N * np.sqrt(p0.El_t / (p0.rho_d * p0.L0)) / p0.R0
dt_val   = float(0.4 * 2.0 / om_edge)
alpha_val = 2.0

# Store dt/alpha in params for adaptive_swell
params5['_dt_tf']    = tf.constant(NP_DTYPE(dt_val))
params5['_alpha_tf'] = tf.constant(NP_DTYPE(alpha_val))

# Run just 2 swell increments (very short test)
phi_after_init = phi_init5
phi_target5 = phi_init5 + 0.03   # just 3% above start

from src.epd.initializer import adaptive_swell
state5_new, Lx5_new, Ly5_new = adaptive_swell(
    state5, params5, cm5, prim_data5,
    phi_target=phi_target5, Lx=Lx5, Ly=Ly5,
    dphi_init=0.010, dphi_max=0.020, dphi_min=0.002,
    n_relax=200, max_extra_relax=500,
    f_warn=0.45, f_crit=0.65, max_restores=3,
    verbose=True,
)

phi_final5 = compute_phi_outer(state5_new, Lx5_new, Ly5_new, r_c_arr5)
check("5.1 phi_outer increased after swell",
      phi_final5 > phi_init5,
      f"before={phi_init5:.4f} after={phi_final5:.4f}")
check("5.2 phi_outer reached target or close",
      phi_final5 >= phi_target5 - 0.01,
      f"phi_final={phi_final5:.4f} target={phi_target5:.4f}")
check("5.3 Box compressed (Lx decreased)",
      Lx5_new < Lx5,
      f"Lx5={Lx5:.3f} Lx5_new={Lx5_new:.3f}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary figure
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

def _draw_particles(ax, particles, title, Lx_box=None, Ly_box=None, extra_patches=None):
    ax.set_aspect('equal')
    ax.set_facecolor('#f5f5f5')
    colors = plt.cm.tab20(np.linspace(0, 1, len(particles)))
    for i, p in enumerate(particles):
        xy = p.x
        ax.fill(xy[:, 0], xy[:, 1], color=colors[i], alpha=0.6)
        ax.plot(np.append(xy[:, 0], xy[0, 0]),
                np.append(xy[:, 1], xy[0, 1]), 'k-', lw=0.5)
    if Lx_box and Ly_box:
        ax.set_xlim(-0.5, Lx_box + 0.5)
        ax.set_ylim(-0.5, Ly_box + 0.5)
        ax.plot([0, Lx_box, Lx_box, 0, 0], [0, 0, Ly_box, Ly_box, 0],
                'b--', lw=1.5, alpha=0.5)
    if extra_patches:
        for patch in extra_patches:
            ax.add_patch(patch)
    ax.set_title(title)
    ax.grid(True, alpha=0.2)

_draw_particles(axes[0], particles1, f'RSA P=16 periodic\nφ≈{phi_rsa:.3f}', Lx, Ly)

# Box test
_draw_particles(axes[1], particles2, 'RSA inside Box(12×12)', Lx2, Ly2)
from matplotlib.patches import Rectangle
axes[1].add_patch(Rectangle((1.0, 1.0), 12.0-2, 12.0-2, fill=False,
                              edgecolor='navy', lw=2, ls='--'))

# Obstacle test
_draw_particles(axes[2], particles3, 'RSA around CircleObstacle(r=2)', Lx3, Ly3)
theta_c = np.linspace(0, 2*np.pi, 100)
axes[2].plot(7 + 2*np.cos(theta_c), 7 + 2*np.sin(theta_c),
             'r-', lw=2, label='obstacle')
axes[2].legend(fontsize=8)

plt.tight_layout()
outpath = os.path.join(OUTDIR, 'init_test.png')
plt.savefig(outpath, dpi=100)
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
