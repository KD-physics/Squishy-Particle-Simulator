"""
test_E_figures.py — Generate paper figures for the viscous drag section.

Figure E1: Oh sweep — measured vs. theoretical terminal velocity (emulsion + elastic)
Figure E2: Constant background flow tracking — v_cm_x(t) with tau_drag annotated

Outputs (relative to repo root):
  papers/summary_of_methods/figures/drag_terminal_velocity.png
  papers/summary_of_methods/figures/drag_flow_tracking.png

Usage:
    python src/validation/test_E_figures.py
or from papers/summary_of_methods/scripts/:
    python test_E_figures.py
"""

import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Wall
from src.epd.system import System

FIG_DIR = "papers/summary_of_methods/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# ── shared helpers ────────────────────────────────────────────────────────────

def _build_single(Oh, g, ptype='emulsion', nu=0.5, N=36, LX=20.0, LY=60.0,
                  with_floor=True):
    sys_ = System(LX, LY, periodic_x=False, periodic_y=False, g=g)
    if with_floor:
        w = Wall((0, 0), (LX, 0), normal=(0, 1))
        sys_.add_object(w)
    if ptype == 'emulsion':
        spec = ParticleSpec(count=1, type='emulsion', gamma=1.0, kappa=0.02,
                            N_nodes=N, Oh=Oh)
    else:
        spec = ParticleSpec(count=1, type='elastic', nu=nu, N_nodes=N, Oh=Oh)
    sys_.add_particles(spec)
    sys_.initialize(phi_target=0.80, seed=0, verbose=False,
                    relax_only=True, n_relax_init=0)
    return sys_, spec


def _place_high(sys_, y_target=30.0):
    x_cm_np  = sys_.state['x_cm'].numpy()
    x_all_np = sys_.state['x_all'].numpy()
    dy = y_target - float(x_cm_np[0, 1])
    x_cm_np[0, 1]    += dy
    x_all_np[0, :, 1] += dy
    sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
    sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)


def _plateau_v(sys_, n_run=8000, sample=200):
    sys_.run(n_run, sample_every=sample, verbose=False)
    vy = np.array([fr['x_cm'][0, 1] for fr in sys_.frames])
    dt_s = sample * sys_._dt
    if len(vy) >= 6:
        return abs(float(np.mean(np.diff(vy[-6:]))) / dt_s)
    return 0.0


# ══════════════════════════════════════════════════════════════════════════════
# Figure E1 — Terminal velocity: Oh sweep
# ══════════════════════════════════════════════════════════════════════════════

print("── Figure E1: terminal velocity Oh sweep ────────────────────────────")

g = 0.05

# Emulsion: sweep over 5 Oh values
Oh_emul  = [0.10, 0.25, 0.50, 1.00, 2.00]
vt_emul_theory = [g / (2.0 * Oh) for Oh in Oh_emul]
vt_emul_meas   = []

for Oh in Oh_emul:
    print(f"  emulsion Oh={Oh}", flush=True)
    sys_, spec = _build_single(Oh=Oh, g=g, ptype='emulsion', LX=20.0, LY=60.0)
    _place_high(sys_, 30.0)
    vt_emul_meas.append(_plateau_v(sys_))

# Elastic: two spot-checks
Oh_elas  = [0.50, 1.00]
vt_elas_meas   = []
vt_elas_theory = []
for Oh in Oh_elas:
    print(f"  elastic  Oh={Oh}", flush=True)
    sys_, spec = _build_single(Oh=Oh, g=g, ptype='elastic', nu=0.5,
                                LX=20.0, LY=60.0)
    _place_high(sys_, 30.0)
    v_meas = _plateau_v(sys_)
    v_th   = spec.terminal_velocity(g)
    vt_elas_meas.append(v_meas)
    vt_elas_theory.append(v_th)
    print(f"    v_theory={v_th:.5f}  v_meas={v_meas:.5f}")

# ── plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5.5, 5.0))

# y = x reference line
v_all = vt_emul_theory + vt_elas_theory
v_max = max(v_all) * 1.15
v_min = 0.0
ax.plot([v_min, v_max], [v_min, v_max], 'k--', lw=1.0, label=r'$v_{\rm meas} = v_{\rm theory}$')

# Emulsion points
ax.scatter(vt_emul_theory, vt_emul_meas,
           s=70, color='#1f77b4', marker='o', zorder=5, label='Emulsion')
for th, me, Oh in zip(vt_emul_theory, vt_emul_meas, Oh_emul):
    ax.annotate(f'$Oh={Oh}$', xy=(th, me), xytext=(5, 3),
                textcoords='offset points', fontsize=7.5, color='#1f77b4')

# Elastic points
ax.scatter(vt_elas_theory, vt_elas_meas,
           s=70, color='#d62728', marker='s', zorder=5, label='Elastic')
for th, me, Oh in zip(vt_elas_theory, vt_elas_meas, Oh_elas):
    ax.annotate(f'$Oh={Oh}$', xy=(th, me), xytext=(5, -10),
                textcoords='offset points', fontsize=7.5, color='#d62728')

ax.set_xlabel(r'$v_{\rm theory} = g\,/\,(2\,Oh\,v_{\rm ref}/R_0)$', fontsize=11)
ax.set_ylabel(r'$v_{\rm meas}$', fontsize=11)
ax.set_title(r'Terminal velocity: measured vs.\ theory', fontsize=11)
ax.set_xlim(v_min, v_max)
ax.set_ylim(v_min, v_max)
ax.legend(fontsize=9, loc='upper left')
ax.set_aspect('equal')
fig.tight_layout()
path_e1 = os.path.join(FIG_DIR, "drag_terminal_velocity.png")
fig.savefig(path_e1, dpi=150)
plt.close(fig)
print(f"  Saved {path_e1}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure E2 — Constant background flow tracking
# ══════════════════════════════════════════════════════════════════════════════

print("── Figure E2: background flow tracking ──────────────────────────────")

U0  = 0.10
Oh  = 0.50

sys_, spec = _build_single(Oh=Oh, g=0.0, ptype='emulsion',
                            LX=40.0, LY=40.0, with_floor=False)

# Centre particle
x_cm_np  = sys_.state['x_cm'].numpy()
x_all_np = sys_.state['x_all'].numpy()
dx = 20.0 - float(x_cm_np[0, 0])
dy = 20.0 - float(x_cm_np[0, 1])
x_cm_np[0]  += [dx, dy]
x_all_np[0] += [dx, dy]
sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

sys_.U_background = ('constant', {'U': (U0, 0.0)})

xi       = spec.derived['xi']
tau_drag = spec.R0_mean / (2.0 * xi)
n_tau    = int(5.0 * tau_drag / sys_._dt)
n_run    = max(n_tau, 4000)
sample   = max(1, n_run // 120)

sys_.run(n_run, sample_every=sample, verbose=False)

t_arr  = np.array([fr['t'] for fr in sys_.frames])
vx_arr = np.array([
    float(np.diff([sys_.frames[max(0, i-1)]['x_cm'][0, 0],
                   fr['x_cm'][0, 0]])[0]) / (sample * sys_._dt)
    for i, fr in enumerate(sys_.frames)
])
# Smooth with a 5-point running mean to remove noise
from numpy.lib.stride_tricks import sliding_window_view
pad = 2
vx_smooth = np.convolve(vx_arr, np.ones(5)/5, mode='same')

# Theoretical exponential: v_x(t) = U0 * (1 - exp(-t/tau_drag))
t_theory = np.linspace(0, t_arr[-1], 300)
vx_theory = U0 * (1.0 - np.exp(-t_theory / tau_drag))

# ── plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6.5, 4.0))
ax.plot(t_arr, vx_smooth, color='#1f77b4', lw=1.5, label=r'$v_{{\rm cm},x}$ (measured)')
ax.plot(t_theory, vx_theory, 'k--', lw=1.2, label=r'$U_0(1-e^{-t/\tau_{\rm drag}})$')
ax.axhline(U0, color='gray', ls=':', lw=1.0, label=f'$U_0 = {U0}$')
ax.axvline(tau_drag, color='#ff7f0e', ls='--', lw=1.0,
           label=fr'$\tau_{{\rm drag}} = {tau_drag:.2f}$')
ax.axvline(3*tau_drag, color='#ff7f0e', ls=':', lw=0.8, alpha=0.5,
           label=r'$3\tau_{\rm drag}$')
ax.set_xlabel(r'$t$', fontsize=11)
ax.set_ylabel(r'$v_{{\rm cm},x}$', fontsize=11)
ax.set_title(r'Exponential approach to background flow ($Oh = 0.5$, $U_0 = 0.1$)', fontsize=10)
ax.legend(fontsize=8.5, loc='lower right')
ax.set_xlim(0, t_arr[-1])
ax.set_ylim(0, U0 * 1.12)
fig.tight_layout()
path_e2 = os.path.join(FIG_DIR, "drag_flow_tracking.png")
fig.savefig(path_e2, dpi=150)
plt.close(fig)
print(f"  Saved {path_e2}")

print("\nFigure generation complete.")
print(f"  {path_e1}")
print(f"  {path_e2}")
