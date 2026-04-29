"""
test_E_movies.py — Render supplemental movies for the viscous drag section.

Movie S12: Single emulsion droplet falling under gravity with Stokes drag,
           reaching terminal velocity (Oh=0.25, g=0.05).

Movie S13: Single emulsion particle in simple shear flow (Γ=0.05, Oh=0.5),
           showing CM drift and clockwise rotation.

Outputs:
  papers/summary_of_methods/movies/S12_terminal_velocity_drag.gif
  papers/summary_of_methods/movies/S13_shear_rotation.gif

Usage:
    python src/validation/test_E_movies.py
"""

import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
import imageio.v2 as imageio
import io

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Wall
from src.epd.system import System

MOV_DIR = "papers/summary_of_methods/movies"
os.makedirs(MOV_DIR, exist_ok=True)


def _render_frame_terminal(nodes, x_cm, t, v_cm_y, v_theory, LX, LY,
                            floor_y=0.0, title_extra=""):
    """Render one frame of the terminal-velocity movie."""
    fig, axes = plt.subplots(1, 2, figsize=(8, 5),
                              gridspec_kw={'width_ratios': [2, 1]})
    ax, ax2 = axes

    # Left: particle in box
    ax.set_xlim(0, LX); ax.set_ylim(-0.5, LY)
    ax.set_aspect('equal')
    ax.axhline(floor_y, color='#555555', lw=2.0)
    poly = plt.Polygon(nodes, closed=True, fill=True, facecolor='#aec6e8',
                       edgecolor='#1f77b4', lw=1.5, alpha=0.85)
    ax.add_patch(poly)
    ax.plot(*x_cm, 'r+', ms=8, lw=1.5)
    ax.set_xlabel(r'$x$', fontsize=10)
    ax.set_ylabel(r'$y$', fontsize=10)
    ax.set_title(fr'$t = {t:.1f}$    $v_{{cm,y}} = {abs(v_cm_y):.4f}$'
                 + (f'\n{title_extra}' if title_extra else ''), fontsize=9)

    # Right: velocity gauge
    v_max = max(v_theory * 1.3, 0.01)
    ax2.barh(0, abs(v_cm_y), color='#1f77b4', height=0.4, label='measured')
    ax2.barh(-0.5, v_theory, color='k', height=0.05, alpha=0.5)
    ax2.axvline(v_theory, color='k', ls='--', lw=1.0, label=f'$v_t={v_theory:.4f}$')
    ax2.set_xlim(0, v_max)
    ax2.set_yticks([])
    ax2.set_xlabel(r'$|v_{cm,y}|$', fontsize=9)
    ax2.set_title('Speed', fontsize=9)
    ax2.legend(fontsize=7, loc='lower right')

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=90)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


def _render_frame_shear(nodes, x_cm, theta, t, omega, omega_theory,
                         LX, LY, shear_rate):
    """Render one frame of the shear-flow movie."""
    fig, axes = plt.subplots(1, 2, figsize=(8, 4),
                              gridspec_kw={'width_ratios': [3, 1]})
    ax, ax2 = axes

    # Left: particle in wide box with shear-flow arrows
    ax.set_xlim(x_cm[0] - 4, x_cm[0] + 4)
    ax.set_ylim(0, LY)
    ax.set_aspect('equal')
    # Shear-flow background arrows
    for y_arr in np.linspace(1, LY - 1, 8):
        u_arr = shear_rate * y_arr
        ax.annotate('', xy=(x_cm[0] + u_arr * 0.4, y_arr),
                    xytext=(x_cm[0], y_arr),
                    arrowprops=dict(arrowstyle='->', color='lightgray', lw=1.0))
    poly = plt.Polygon(nodes, closed=True, fill=True, facecolor='#ffbb78',
                       edgecolor='#d62728', lw=1.5, alpha=0.85)
    ax.add_patch(poly)
    ax.plot(*x_cm, 'r+', ms=8, lw=1.5)
    ax.set_xlabel(r'$x$', fontsize=10)
    ax.set_ylabel(r'$y$', fontsize=10)
    ax.set_title(fr'Shear flow $\Gamma={shear_rate}$    $t={t:.1f}$', fontsize=9)

    # Right: rotation rate gauge
    ax2.barh(0, omega, color='#d62728', height=0.4, label=fr'$\omega={omega:.4f}$')
    ax2.axvline(omega_theory, color='k', ls='--', lw=1.0,
                label=fr'$\omega_{{th}}={omega_theory:.4f}$')
    omega_range = max(abs(omega_theory) * 1.5, 0.01)
    ax2.set_xlim(-omega_range, omega_range)
    ax2.set_yticks([])
    ax2.set_xlabel(r'$\omega$', fontsize=9)
    ax2.set_title('Rotation', fontsize=9)
    ax2.legend(fontsize=7, loc='upper right')
    ax2.axvline(0, color='gray', lw=0.5)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=90)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


# ══════════════════════════════════════════════════════════════════════════════
# Movie S12 — Terminal velocity with drag (Oh=0.25, g=0.05)
# ══════════════════════════════════════════════════════════════════════════════

print("── Movie S12: terminal velocity ─────────────────────────────────────")

Oh12 = 0.25
g12  = 0.05
LX12, LY12 = 20.0, 60.0

sys12 = System(LX12, LY12, periodic_x=False, periodic_y=False, g=g12)
w12   = Wall((0, 0), (LX12, 0), normal=(0, 1))
sys12.add_object(w12)
spec12 = ParticleSpec(count=1, type='emulsion', gamma=1.0, kappa=0.02,
                      N_nodes=36, Oh=Oh12)
sys12.add_particles(spec12)
sys12.initialize(phi_target=0.80, seed=0, verbose=False,
                 relax_only=True, n_relax_init=0)

# Place near top
x_cm_np  = sys12.state['x_cm'].numpy()
x_all_np = sys12.state['x_all'].numpy()
dy = 45.0 - float(x_cm_np[0, 1])
x_cm_np[0, 1]    += dy
x_all_np[0, :, 1] += dy
sys12._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
sys12._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

v_theory12 = spec12.terminal_velocity(g12)
print(f"  v_theory = {v_theory12:.4f}")

# Run and collect frames
n_total = 10000
chunk   = 200
sample  = 40
frames12 = []

for i in range(n_total // chunk):
    sys12.run(chunk, sample_every=sample, verbose=False)
    fr = sys12.frames[-1]
    nodes   = fr['x_all'][0]
    x_cm    = fr['x_cm'][0]
    t_now   = fr['t']
    v_cm_y  = float(sys12._state['v_cm'].numpy()[0, 1])
    img = _render_frame_terminal(nodes, x_cm, t_now, v_cm_y, v_theory12,
                                  LX12, LY12, floor_y=0.0,
                                  title_extra=f'Oh={Oh12}  g={g12}  v_theory={v_theory12:.4f}')
    frames12.append(img)
    sys12.clear_recording()
    if x_cm[1] < 2.5:   # particle reached floor region
        print(f"  Reached floor at t={t_now:.1f}, stopping early.")
        break

path_s12 = os.path.join(MOV_DIR, "S12_terminal_velocity_drag.gif")
imageio.mimwrite(path_s12, frames12, duration=0.12, loop=0)
print(f"  Saved {path_s12} ({len(frames12)} frames)")


# ══════════════════════════════════════════════════════════════════════════════
# Movie S13 — Shear flow: drift and rotation (Γ=0.05, Oh=0.5)
# ══════════════════════════════════════════════════════════════════════════════

print("── Movie S13: shear flow rotation ───────────────────────────────────")

shear13 = 0.05
Oh13    = 0.50
LX13, LY13 = 200.0, 40.0

sys13 = System(LX13, LY13, periodic_x=False, periodic_y=False, g=0.0)
spec13 = ParticleSpec(count=1, type='emulsion', gamma=1.0, kappa=0.02,
                      N_nodes=36, Oh=Oh13)
sys13.add_particles(spec13)
sys13.initialize(phi_target=0.80, seed=0, verbose=False,
                 relax_only=True, n_relax_init=0)

# Centre at y=20, x=20 (leaves room to drift right)
x_cm_np  = sys13.state['x_cm'].numpy()
x_all_np = sys13.state['x_all'].numpy()
dx13 = 20.0 - float(x_cm_np[0, 0])
dy13 = 20.0 - float(x_cm_np[0, 1])
x_cm_np[0]  += [dx13, dy13]
x_all_np[0] += [dx13, dy13]
sys13._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
sys13._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)

sys13.U_background = ('shear', {'rate': shear13})

omega_theory13 = -shear13 / 2.0
print(f"  omega_theory = {omega_theory13:.4f}")

n_total13 = 6000
chunk13   = 200
sample13  = 50
frames13  = []

for i in range(n_total13 // chunk13):
    sys13.run(chunk13, sample_every=sample13, verbose=False)
    fr = sys13.frames[-1]
    nodes   = fr['x_all'][0]
    x_cm    = fr['x_cm'][0]
    t_now   = fr['t']
    omega   = float(sys13._state['omega'].numpy()[0])
    img = _render_frame_shear(nodes, x_cm, fr.get('theta', np.array([0.0]))[0],
                               t_now, omega, omega_theory13, LX13, LY13, shear13)
    frames13.append(img)
    sys13.clear_recording()

path_s13 = os.path.join(MOV_DIR, "S13_shear_rotation.gif")
imageio.mimwrite(path_s13, frames13, duration=0.14, loop=0)
print(f"  Saved {path_s13} ({len(frames13)} frames)")

print("\nMovie generation complete.")
print(f"  {path_s12}")
print(f"  {path_s13}")
