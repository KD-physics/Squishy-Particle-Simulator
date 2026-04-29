"""
make_swelling_gif.py — Swelling initialization: two disks grow inside a fixed-wall box

Physical setup
--------------
Two equal elastic disks inside a rigid box of FIXED width.
Disk radius R grows from R_start (no contact) to R_target.
As R grows:
  - Disks first approach walls and each other
  - Contact forms at R ≈ R_target - δ_target/4
  - Both wall-disk and disk-disk forces build simultaneously
  - Bisection (run_two_disk_wall_contact) finds equilibrium cx_shift at each step

Scene layout
------------
Left (60%): disks + walls.
  - Grey walls FIXED at ±W
  - Dashed ring: current undeformed circle (grows each frame)
  - Dotted ring: final target circle
  - Perimeter coloured by |u| (hot colormap)
  - Red dots: disk-disk contact nodes
  - Blue dots: disk-wall contact nodes
Right top (40%): F_disk-disk vs R
Right bottom   : F_wall vs R

Run:
  source .venv/bin/activate
  python src/visualization/make_swelling_gif.py
"""

import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os, sys, time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from src.simulation.two_disk_contact import run_two_disk_wall_contact
from src.simulation.disk_mesh import make_disk_mesh

os.makedirs('results', exist_ok=True)

# ── Palette ────────────────────────────────────────────────────────────────────
BG          = '#1a1a2e'
PANEL       = '#16213e'
WALL_COLOR  = '#7f8c8d'
COL_D1      = '#2ecc71'    # disk 1 (left)
COL_D2      = '#e67e22'    # disk 2 (right)
COL_DD      = '#e74c3c'    # disk-disk contact
COL_DW      = '#3498db'    # disk-wall contact


# ── Simulation ─────────────────────────────────────────────────────────────────

def _make_no_contact_frame(R_step, E, nu, W, N_perimeter):
    """Dummy result for frames where disks haven't yet reached the walls."""
    mesh1 = make_disk_mesh(R=R_step, N_perimeter=N_perimeter)
    mesh2 = make_disk_mesh(R=R_step, N_perimeter=N_perimeter)
    mesh1['center'] = np.array([-R_step, 0.0])
    mesh2['center'] = np.array([+R_step, 0.0])
    N1  = mesh1['vertices'].shape[0]
    N2  = mesh2['vertices'].shape[0]
    Np  = len(mesh1['perimeter_ids'])
    return {
        'converged': True,
        'F_contact': 0.0, 'F_wall1': 0.0, 'F_wall2': 0.0,
        'cx_shift_eq': 0.0, 'delta_over_R': 0.0,
        'u1': np.zeros((N1, 2)), 'u2': np.zeros((N2, 2)),
        'mesh1': mesh1, 'mesh2': mesh2,
        'node_forces1':      np.zeros((Np, 2)),
        'node_forces2':      np.zeros((Np, 2)),
        'node_forces_wall1': np.zeros((Np, 2)),
        'node_forces_wall2': np.zeros((Np, 2)),
        'R': R_step, 'R2': R_step, 'E': E, 'nu': nu,
        'N_perimeter': N_perimeter,
    }


def collect_frames(R_start, R_target, delta_target, E, nu, N_steps, N_perimeter):
    """
    Collect one result per frame.

    Fixed-box half-width:  W = 2·R_target − δ_target/2
    At radius R_step, effective compression in the box:
      δ_step = max(0,  4·R_step − 4·R_target + δ_target)
    Contact onset: R_step > R_target − δ_target/4

    run_two_disk_wall_contact handles bisection for the equilibrium
    cx_shift at each step.
    """
    W       = 2.0 * R_target - delta_target / 2.0   # fixed wall half-width
    R_steps = np.linspace(R_start, R_target, N_steps + 1)[1:]
    frames  = []

    for R_step in R_steps:
        delta_step = max(0.0, 4.0 * R_step - 4.0 * R_target + delta_target)
        d_over_R   = delta_step / R_step
        in_contact = d_over_R > 1e-6

        if in_contact:
            res = run_two_disk_wall_contact(
                R=R_step, E=E, nu=nu,
                delta_over_R=d_over_R,
                N_perimeter=N_perimeter,
                verbose=False,
            )
        else:
            res = _make_no_contact_frame(R_step, E, nu, W, N_perimeter)

        res['_R_step']     = R_step
        res['_in_contact'] = in_contact
        res['_W']          = W

        if in_contact:
            print(f'    R={R_step:.4f}  d/R={d_over_R:.4f}  '
                  f'F_dd={res["F_contact"]:.1f}  F_wall={res["F_wall1"]:.1f}  '
                  f'cx_shift={res["cx_shift_eq"]:.4f}')
        else:
            print(f'    R={R_step:.4f}  [no contact]')

        frames.append(res)

    return frames, W


# ── Drawing ────────────────────────────────────────────────────────────────────

def _draw_disk(ax, mesh, u, nf_dd, nf_wall, col, R_target, cx_target):
    """Draw one disk: reference ring, deformed perimeter, contact highlights."""
    cx  = mesh['center']
    pid = mesh['perimeter_ids']
    pp  = mesh['perimeter_pos']
    R   = np.max(np.linalg.norm(pp, axis=1))   # current disk radius

    theta = np.linspace(0, 2 * np.pi, 200)

    # Dashed ring: current undeformed radius
    ax.plot(cx[0] + R * np.cos(theta),
            cx[1] + R * np.sin(theta),
            '--', color=col, alpha=0.28, lw=1.0, zorder=3)

    # Dotted ring: final target radius
    ax.plot(cx_target + R_target * np.cos(theta),
            R_target * np.sin(theta),
            ':', color=col, alpha=0.12, lw=0.8, zorder=2)

    # Deformed fill
    xy_def    = cx + pp + u[pid]
    xy_closed = np.vstack([xy_def, xy_def[0]])
    ax.fill(xy_closed[:, 0], xy_closed[:, 1], alpha=0.18, color=col, zorder=4)

    # Perimeter segments coloured by |u|
    umag = np.linalg.norm(u[pid], axis=1)
    umax = umag.max() + 1e-30
    N_p  = len(pid)
    for j in range(N_p):
        jn = (j + 1) % N_p
        c  = plt.cm.YlOrRd(umag[j] / umax * 0.88 + 0.05)
        ax.plot([cx[0] + pp[j,  0] + u[pid[j],  0],
                 cx[0] + pp[jn, 0] + u[pid[jn], 0]],
                [cx[1] + pp[j,  1] + u[pid[j],  1],
                 cx[1] + pp[jn, 1] + u[pid[jn], 1]],
                color=c, lw=2.2, solid_capstyle='round', alpha=0.92, zorder=5)

    # Disk-disk contact nodes (red)
    Fmag_dd = np.linalg.norm(nf_dd, axis=1)
    if Fmag_dd.max() > 1e-12:
        in_c = Fmag_dd > 0.03 * Fmag_dd.max()
        if in_c.sum() >= 1:
            xs = cx[0] + pp[in_c, 0] + u[pid[in_c], 0]
            ys = cx[1] + pp[in_c, 1] + u[pid[in_c], 1]
            ax.scatter(xs, ys, s=22, c=COL_DD, zorder=7, alpha=0.90,
                       edgecolors='none', label='disk–disk')

    # Disk-wall contact nodes (blue)
    Fmag_w = np.linalg.norm(nf_wall, axis=1)
    if Fmag_w.max() > 1e-12:
        in_w = Fmag_w > 0.03 * Fmag_w.max()
        if in_w.sum() >= 1:
            xs = cx[0] + pp[in_w, 0] + u[pid[in_w], 0]
            ys = cx[1] + pp[in_w, 1] + u[pid[in_w], 1]
            ax.scatter(xs, ys, s=22, c=COL_DW, zorder=7, alpha=0.90,
                       edgecolors='none', label='disk–wall')

    # Centre marker
    ax.plot(cx[0], cx[1], '+', color=col, ms=5, mew=1.4, zorder=8, alpha=0.6)


def _draw_scene(ax, res, R_target, delta_target):
    ax.set_facecolor(PANEL)

    W     = res['_W']
    R_step = res['_R_step']
    pad   = R_target * 0.18
    xlim  = (-W - pad, W + pad)
    ylim  = (-R_target - pad, R_target + pad)

    # Fixed walls (grey bands)
    ax.axvspan(xlim[0], -W, color=WALL_COLOR, alpha=0.80, zorder=2)
    ax.axvspan(+W, xlim[1],  color=WALL_COLOR, alpha=0.80, zorder=2)
    ax.axhline(+R_target, color=WALL_COLOR, lw=2.0, alpha=0.55, zorder=2)
    ax.axhline(-R_target, color=WALL_COLOR, lw=2.0, alpha=0.55, zorder=2)

    # Wall labels
    ax.text(-W - pad*0.5, 0, 'wall', color='white', fontsize=6,
            ha='center', va='center', rotation=90, alpha=0.6, zorder=3)
    ax.text(+W + pad*0.5, 0, 'wall', color='white', fontsize=6,
            ha='center', va='center', rotation=90, alpha=0.6, zorder=3)

    # Contact-plane divider (only when in contact)
    if res['_in_contact']:
        ax.axvline(0, color='white', lw=0.5, ls=':', alpha=0.20, zorder=2)

    # Disk 1 (left, green)
    cx1_target = -(R_target - res.get('cx_shift_eq', 0.0))
    _draw_disk(ax, res['mesh1'], res['u1'],
               res['node_forces1'], res['node_forces_wall1'],
               COL_D1, R_target, cx1_target)

    # Disk 2 (right, orange)
    cx2_target = +(R_target - res.get('cx_shift_eq', 0.0))
    _draw_disk(ax, res['mesh2'], res['u2'],
               res['node_forces2'], res['node_forces_wall2'],
               COL_D2, R_target, cx2_target)

    # Annotation box
    status = 'IN CONTACT' if res['_in_contact'] else 'no contact'
    fdd_str  = f'F_dd   = {res["F_contact"]:.1f} N/m' if res['_in_contact'] else 'F_dd   = 0'
    fwall_str = f'F_wall = {res["F_wall1"]:.1f} N/m' if res['_in_contact'] else 'F_wall = 0'
    d_str    = (f'δ/R = {res["delta_over_R"]:.4f}'
                if res['_in_contact'] else 'δ/R = —')
    label = (f'R = {R_step:.4f} m  [{status}]\n'
             f'{d_str}\n'
             f'{fdd_str}\n'
             f'{fwall_str}')
    ax.text(0.02, 0.97, label, transform=ax.transAxes,
            color='white', fontsize=7.5, va='top', ha='left', family='monospace',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='#0d0d1a', alpha=0.70),
            zorder=9)

    # Legend (contact markers)
    if res['_in_contact']:
        from matplotlib.lines import Line2D
        handles = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_DD,
                   markersize=6, label='disk–disk contact', linestyle='None'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_DW,
                   markersize=6, label='disk–wall contact', linestyle='None'),
        ]
        ax.legend(handles=handles, loc='lower right', fontsize=6,
                  facecolor='#111', edgecolor='#333', labelcolor='white',
                  framealpha=0.7)

    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x (m)', color='white', fontsize=7)
    ax.set_ylabel('y (m)', color='white', fontsize=7)
    ax.tick_params(colors='white', labelsize=6)
    for sp in ax.spines.values():
        sp.set_edgecolor('#444')


def _draw_F_curve(ax, frames, current_idx, key, ylabel, col, label):
    ax.set_facecolor(PANEL)

    R_arr = np.array([f['_R_step']   for f in frames])
    F_arr = np.array([f[key]         for f in frames])

    # Full curve (faint)
    ax.plot(R_arr, F_arr, '-', color=col, alpha=0.22, lw=1.5, zorder=2)

    # Path so far (bright)
    ax.plot(R_arr[:current_idx + 1], F_arr[:current_idx + 1],
            '-', color=col, alpha=0.85, lw=2.2, zorder=3)

    # Current point
    ax.plot(R_arr[current_idx], F_arr[current_idx],
            'o', ms=8, markerfacecolor=col, markeredgecolor='white',
            mew=1.2, zorder=5)

    # Contact onset line
    F_max = F_arr.max()
    contact_frames = [f for f in frames if f['_in_contact']]
    if contact_frames:
        R_onset = contact_frames[0]['_R_step']
        ax.axvline(R_onset, color='white', ls=':', lw=0.9, alpha=0.35)
        ax.text(R_onset, F_max * 0.02, 'onset', color='white',
                fontsize=5.5, ha='center', va='bottom', alpha=0.50)

    ax.set_xlabel('R (m)', color='white', fontsize=7)
    ax.set_ylabel(ylabel, color='white', fontsize=7)
    ax.set_title(label, color=col, fontsize=8, pad=2)
    ax.set_xlim(R_arr[0] - 0.002, R_arr[-1] + 0.002)
    ax.set_ylim(-F_max * 0.06, F_max * 1.15)
    ax.tick_params(colors='white', labelsize=6)
    for sp in ax.spines.values():
        sp.set_edgecolor('#444')


# ══════════════════════════════════════════════════════════════════════════════

def make_swelling_gif():
    print('\n── Swelling GIF: two disks in a fixed-wall box ──')
    print('   polystyrene E=5×10⁵, ν=0.33,  δ_target/R = 0.10')

    R_target     = 1.0
    delta_target = 0.10     # larger than before → more wall contact visible
    E, nu        = 5e5, 0.33
    N_perimeter  = 120
    R_start      = 0.94     # contact onset at R ≈ 0.975
    N_steps      = 14

    frames, W = collect_frames(
        R_start, R_target, delta_target, E, nu, N_steps, N_perimeter)

    n_frames = len(frames)
    print(f'\n  Rendering {n_frames} frames  (box half-width W={W:.3f} m)...')

    fig = plt.figure(figsize=(12, 6.0))
    fig.patch.set_facecolor(BG)

    gs = fig.add_gridspec(2, 2, width_ratios=[2.4, 1.0],
                          hspace=0.38, wspace=0.28,
                          left=0.04, right=0.97, top=0.90, bottom=0.09)
    ax_scene = fig.add_subplot(gs[:, 0])    # full-height scene
    ax_Fdd   = fig.add_subplot(gs[0, 1])   # F disk-disk vs R
    ax_Fw    = fig.add_subplot(gs[1, 1])   # F wall vs R

    fig.suptitle(
        'Particle Swelling in Fixed Box — Disk-Disk & Disk-Wall Contact\n'
        'polystyrene  E=5×10⁵ Pa  ν=0.33   δ_target/R=0.10',
        color='white', fontsize=10, y=0.97)

    def draw_frame(i):
        ax_scene.cla()
        ax_Fdd.cla()
        ax_Fw.cla()

        _draw_scene(ax_scene, frames[i], R_target, delta_target)
        ax_scene.set_title(
            f'Walls fixed at x = ±{W:.3f} m  |  '
            f'{"-- dashed: current R  ·· dotted: final R"}',
            color='white', fontsize=7.5, pad=4)

        _draw_F_curve(ax_Fdd, frames, i, 'F_contact',
                      'F (N/m)', COL_DD, 'Disk–disk contact force')
        _draw_F_curve(ax_Fw,  frames, i, 'F_wall1',
                      'F (N/m)', COL_DW, 'Disk–wall contact force')

    anim = FuncAnimation(fig, draw_frame, frames=n_frames, interval=500)
    path = 'results/swelling_two_disk.gif'
    anim.save(path, writer=PillowWriter(fps=2.0), dpi=120)
    plt.close(fig)
    print(f'  Saved {path}')


if __name__ == '__main__':
    t0 = time.time()
    make_swelling_gif()
    print(f'\nDone in {time.time()-t0:.0f}s')
