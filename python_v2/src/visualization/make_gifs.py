"""
make_gifs.py — Animated GIFs of two-disk contact deformation

Physical setup: two disks in a box. Left and right walls translate inward,
squeezing the disks together. Walls are rigid; disks deform to comply.

GIF 1: delta_ramp_rubber_jello.gif  — rubber/jello, walls compress 0→15%
GIF 2: material_comparison.gif      — 5 materials side-by-side, same compression
GIF 3: size_ratio_4to1.gif          — 4:1 (glass + rubber/jello), walls compress

Run:
  source .venv/bin/activate
  python src/visualization/make_gifs.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch
import os, sys, time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from src.simulation.two_disk_contact import run_two_disk_contact_robust

os.makedirs('results', exist_ok=True)

BG   = '#1a1a2e'
PANEL= '#16213e'
WALL_COLOR  = '#95a5a6'
WALL_ALPHA  = 0.85
ARROW_COLOR = '#e74c3c'


# ── helpers ───────────────────────────────────────────────────────────────────

def run_batch(configs, verbose=True):
    out = []
    for i, (name, R1, E1, nu1, R2, E2, nu2, d) in enumerate(configs):
        res = run_two_disk_contact_robust(
            R=R1, E=E1, nu=nu1, R2=R2, E2=E2, nu2=nu2,
            delta_over_R=d, N_perimeter=120, verbose=False)
        res['_name'] = name
        out.append(res)
        if verbose:
            print(f'  [{i+1}/{len(configs)}] {name:20s} d={d:.3f}  '
                  f'conv={"Y" if res["converged"] else "N"}  F={res["F_contact"]:.3f}')
    return out


def draw_scene(ax, res, col1='#2980b9', col2='#e74c3c',
               show_pressure_overlay=False, label=''):
    """
    Draw the physical squeezing scene:
      - grey walls on left/right that have moved inward by delta/2
      - top/bottom frame lines
      - reference (undeformed) disk outlines (dashed)
      - deformed disk outlines (solid, filled)
      - perimeter coloured by |u| magnitude
      - contact pressure arrows (optional overlay)
    """
    ax.set_facecolor(PANEL)

    R1 = res['R']
    R2 = res['R2']
    delta = res['delta_over_R'] * R1

    # Wall positions in absolute coords
    # At delta=0: disk1 centre at -R1, disk2 at +R2; left wall at -2R1, right at +2R2
    # At delta>0: walls moved inward by delta/2 each
    x_left_wall  = -(2*R1 - delta/2)   # = cx1 - R1
    x_right_wall = +(2*R2 - delta/2)   # = cx2 + R2

    # Scene extents (fixed across frames so walls visibly move)
    x_wall0_left  = -2 * R1            # initial wall positions
    x_wall0_right = +2 * R2
    pad = max(R1, R2) * 0.18
    xlim = (x_wall0_left - pad, x_wall0_right + pad)
    ylim = (-max(R1, R2) - pad, max(R1, R2) + pad)

    # ── Walls ────────────────────────────────────────────────────────────────
    # Left wall: filled band from scene edge to x_left_wall
    ax.axvspan(xlim[0], x_left_wall, color=WALL_COLOR, alpha=WALL_ALPHA, zorder=2)
    # Right wall
    ax.axvspan(x_right_wall, xlim[1], color=WALL_COLOR, alpha=WALL_ALPHA, zorder=2)
    # Top/bottom frame
    ax.axhline( max(R1,R2), color=WALL_COLOR, lw=2, alpha=0.6, zorder=2)
    ax.axhline(-max(R1,R2), color=WALL_COLOR, lw=2, alpha=0.6, zorder=2)

    # Wall-motion arrows (show direction of squeeze)
    arr_y = max(R1,R2) * 0.75
    if delta > 0.001 * R1:
        ax.annotate('', xy=(x_left_wall + 0.02, arr_y),
                    xytext=(x_left_wall - max(R1,R2)*0.22, arr_y),
                    arrowprops=dict(arrowstyle='->', color=ARROW_COLOR,
                                   lw=1.8, mutation_scale=14), zorder=5)
        ax.annotate('', xy=(x_right_wall - 0.02, arr_y),
                    xytext=(x_right_wall + max(R1,R2)*0.22, arr_y),
                    arrowprops=dict(arrowstyle='->', color=ARROW_COLOR,
                                   lw=1.8, mutation_scale=14), zorder=5)

    # Contact plane
    ax.axvline(0, color='white', lw=0.6, alpha=0.25, ls=':', zorder=3)

    for (mesh, u, col) in [(res['mesh1'], res['u1'], col1),
                           (res['mesh2'], res['u2'], col2)]:
        cx  = mesh['center']
        pid = mesh['perimeter_ids']
        pp  = mesh['perimeter_pos']
        R   = mesh.get('R', np.max(np.linalg.norm(pp, axis=1)))

        # Reference circle (dashed, thin)
        theta = np.linspace(0, 2*np.pi, 200)
        ax.plot(cx[0] + R*np.cos(theta), cx[1] + R*np.sin(theta),
                '--', color=col, alpha=0.22, lw=1, zorder=3)

        # Deformed shape
        xy_def = cx + pp + u[pid]
        xy_def = np.vstack([xy_def, xy_def[0]])
        ax.fill(xy_def[:,0], xy_def[:,1], alpha=0.28, color=col, zorder=4)

        # Perimeter coloured by |u| (hot colormap)
        umag = np.linalg.norm(u[pid], axis=1)
        umax = umag.max() + 1e-30
        for j in range(len(pid)):
            jn = (j+1) % len(pid)
            c = plt.cm.YlOrRd(umag[j]/umax * 0.88 + 0.05)
            ax.plot([cx[0]+pp[j,0]+u[pid[j],0], cx[0]+pp[jn,0]+u[pid[jn],0]],
                    [cx[1]+pp[j,1]+u[pid[j],1], cx[1]+pp[jn,1]+u[pid[jn],1]],
                    color=c, lw=2.5, solid_capstyle='round', alpha=0.9, zorder=5)

    ax.set_xlim(xlim); ax.set_ylim(ylim)
    ax.set_aspect('equal')
    ax.tick_params(colors='white', labelsize=6)
    for sp in ax.spines.values():
        sp.set_edgecolor('#333')

    if label:
        ax.text(0.02, 0.97, label, transform=ax.transAxes,
                color='white', fontsize=7, va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#111', alpha=0.6))


def draw_pressure(ax, res, col, ylim_override=None):
    """Draw contact pressure profile on disk 1 surface."""
    ax.set_facecolor(PANEL)
    R  = res['R']
    nf = res['node_forces1']
    pp = res['mesh1']['perimeter_pos']
    h  = 2 * np.pi * R / 120
    P  = np.abs(nf[:,0]) / h
    y  = pp[:,1]
    mask = P > 0.005 * (P.max() + 1e-30)
    if mask.sum() > 1:
        ys = y[mask]; Ps = P[mask]; si = np.argsort(ys)
        ax.fill_betweenx(ys[si], 0, Ps[si], alpha=0.70, color=col)
        ax.plot(Ps[si], ys[si], color='white', lw=1.2)
    pad = R * 0.18
    ax.set_ylim(ylim_override or (-R-pad, R+pad))
    ax.tick_params(colors='white', labelsize=5)
    for sp in ax.spines.values():
        sp.set_edgecolor('#333')


# ══════════════════════════════════════════════════════════════════════════════
# GIF 1 — rubber/jello delta ramp
# ══════════════════════════════════════════════════════════════════════════════

def make_gif1():
    print('\n── GIF 1: rubber/jello walls-squeezing ramp ──')
    DELTAS = [0.000,                                          # zero-contact reference
              0.002, 0.004, 0.006, 0.010, 0.015, 0.020,
              0.030, 0.040, 0.050, 0.070, 0.10, 0.12, 0.15]
    E, nu  = 1e3, 0.49
    col    = '#e74c3c'

    # delta=0 is a dummy (no contact, zero displacement)
    configs = [('rj', 1.0, E, nu, None, None, None, d)
               for d in DELTAS if d > 0]
    print(f'  Running {len(configs)} simulations...')
    results_nonzero = run_batch(configs, verbose=True)

    # Build full list: prepend a fake zero-displacement result
    zero_res = dict(results_nonzero[0])   # copy mesh/material structure
    zero_res['u1'] = np.zeros_like(results_nonzero[0]['u1'])
    zero_res['u2'] = np.zeros_like(results_nonzero[0]['u2'])
    zero_res['F_contact']  = 0.0
    zero_res['delta_over_R'] = 0.0
    zero_res['node_forces1'] = np.zeros_like(results_nonzero[0]['node_forces1'])
    # Reset centres to zero-overlap positions
    import copy
    zero_res = copy.deepcopy(results_nonzero[0])
    zero_res['u1'][:] = 0; zero_res['u2'][:] = 0
    zero_res['F_contact'] = 0.0; zero_res['delta_over_R'] = 0.0
    zero_res['node_forces1'][:] = 0; zero_res['node_forces2'][:] = 0
    R = 1.0
    zero_res['mesh1']['center'] = np.array([-R, 0.0])
    zero_res['mesh2']['center'] = np.array([+R, 0.0])

    all_results = [zero_res] + results_nonzero

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.8),
                             gridspec_kw={'width_ratios': [2.2, 1]})
    fig.patch.set_facecolor(BG)

    def draw_frame(i):
        res = all_results[i]
        axes[0].cla(); axes[1].cla()
        draw_scene(axes[0], res, col1=col, col2=col,
                   label=(f'E = 1×10³ Pa,  ν = 0.49\n'
                          f'δ/R = {res["delta_over_R"]:.3f}\n'
                          f'F = {res["F_contact"]:.2f} N/m'))
        axes[0].set_xlabel('x (m)', color='white', fontsize=8)
        axes[0].set_ylabel('y (m)', color='white', fontsize=8)
        axes[0].set_title('Rubber / Jello — wall squeezing',
                          color='white', fontsize=9, pad=4)

        draw_pressure(axes[1], res, col)
        axes[1].set_xlabel('Contact pressure (Pa)', color='white', fontsize=7)
        axes[1].set_ylabel('y (m)', color='white', fontsize=7)
        axes[1].set_title('Pressure on\ndisk surface', color='white', fontsize=8, pad=3)

        fig.tight_layout()

    anim = FuncAnimation(fig, draw_frame, frames=len(all_results), interval=350)
    path = 'results/delta_ramp_rubber_jello.gif'
    anim.save(path, writer=PillowWriter(fps=2.5), dpi=110)
    plt.close(fig)
    print(f'  Saved {path}')


# ══════════════════════════════════════════════════════════════════════════════
# GIF 2 — material comparison (5 materials, same compression ramp)
# ══════════════════════════════════════════════════════════════════════════════

def make_gif2():
    print('\n── GIF 2: material comparison ──')

    MATERIALS = [
        ('glass',        1e8, 0.20, '#3498db'),
        ('polystyrene',  1e6, 0.33, '#2ecc71'),
        ('PDMS',         3e4, 0.45, '#f39c12'),
        ('rubber',       3e3, 0.48, '#e67e22'),
        ('rubber/jello', 1e3, 0.49, '#e74c3c'),
    ]
    DELTAS = [0.002, 0.010, 0.020, 0.040, 0.070, 0.10, 0.15]

    configs = [(name, 1.0, E, nu, None, None, None, d)
               for (name, E, nu, col) in MATERIALS for d in DELTAS]
    print(f'  Running {len(configs)} simulations...')
    all_res = run_batch(configs, verbose=True)

    nd = len(DELTAS); nm = len(MATERIALS)
    res_grid = [[all_res[i*nd + j] for j in range(nd)] for i in range(nm)]

    fig, axes = plt.subplots(2, nm, figsize=(15, 6.5),
                             gridspec_kw={'height_ratios': [2.5, 1]})
    fig.patch.set_facecolor(BG)

    def draw_frame(di):
        for mi, (name, E, nu, col) in enumerate(MATERIALS):
            res = res_grid[mi][di]
            axes[0, mi].cla(); axes[1, mi].cla()

            Fnorm = res['F_contact'] / (E * 1.0)
            draw_scene(axes[0, mi], res, col1=col, col2=col,
                       label=f'{name}\nE={E:.0e}, ν={nu}\nF/(ER)={Fnorm:.3e}')
            if mi == 0:
                axes[0, mi].set_ylabel('y (m)', color='white', fontsize=7)

            draw_pressure(axes[1, mi], res, col)
            axes[1, mi].set_xlabel('P (Pa)', color='white', fontsize=6)
            if mi == 0:
                axes[1, mi].set_ylabel('y (m)', color='white', fontsize=6)

        fig.suptitle(
            f'Wall compression  —  δ/R = {DELTAS[di]:.3f}  '
            f'({DELTAS[di]*100:.1f}% strain)',
            color='white', fontsize=11, y=1.00)
        fig.tight_layout(rect=[0, 0, 1, 0.97])

    anim = FuncAnimation(fig, draw_frame, frames=nd, interval=500)
    path = 'results/material_comparison.gif'
    anim.save(path, writer=PillowWriter(fps=2), dpi=100)
    plt.close(fig)
    print(f'  Saved {path}')


# ══════════════════════════════════════════════════════════════════════════════
# GIF 3 — 4:1 size ratio: large glass + small rubber/jello
# ══════════════════════════════════════════════════════════════════════════════

def make_gif3():
    print('\n── GIF 3: 4:1 size ratio (glass + rubber/jello) ──')

    R1, E1, nu1 = 1.0, 1e8, 0.20   # big glass
    R2, E2, nu2 = 0.25, 1e3, 0.49  # small rubber/jello
    DELTAS = [0.000, 0.005, 0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.075]

    configs = [('het', R1, E1, nu1, R2, E2, nu2, d) for d in DELTAS if d > 0]
    print(f'  Running {len(configs)} simulations...')
    results_nz = run_batch(configs, verbose=True)

    import copy
    zero_res = copy.deepcopy(results_nz[0])
    zero_res['u1'][:] = 0; zero_res['u2'][:] = 0
    zero_res['F_contact'] = 0.0; zero_res['delta_over_R'] = 0.0
    zero_res['node_forces1'][:] = 0; zero_res['node_forces2'][:] = 0
    zero_res['mesh1']['center'] = np.array([-R1, 0.0])
    zero_res['mesh2']['center'] = np.array([+R2, 0.0])
    all_results = [zero_res] + results_nz

    fig, axes = plt.subplots(1, 2, figsize=(10, 5.2),
                             gridspec_kw={'width_ratios': [2.8, 1]})
    fig.patch.set_facecolor(BG)

    def draw_frame(i):
        res = all_results[i]
        axes[0].cla(); axes[1].cla()
        draw_scene(axes[0], res, col1='#3498db', col2='#e74c3c',
                   label=(f'Blue:  glass  R=1.0  E=1×10⁸  ν=0.20\n'
                          f'Red:   rubber/jello  R=0.25  E=1×10³  ν=0.49\n'
                          f'δ/R₁={res["delta_over_R"]:.3f}   '
                          f'δ/R₂={res["delta_over_R"]*R1/R2:.3f}   '
                          f'F={res["F_contact"]:.3f} N/m'))
        axes[0].set_xlabel('x (m)', color='white', fontsize=8)
        axes[0].set_ylabel('y (m)', color='white', fontsize=8)
        axes[0].set_title('4:1 size ratio — wall squeezing',
                          color='white', fontsize=9, pad=4)

        # Pressure on SMALL disk (more interesting)
        ax_p = axes[1]; ax_p.set_facecolor(PANEL)
        R2r = res['R2']
        nf2 = res['node_forces2']
        h2  = 2*np.pi*R2r/120
        P2  = np.abs(nf2[:,0])/h2
        y2  = res['mesh2']['perimeter_pos'][:,1]
        mask2 = P2 > 0.01*(P2.max()+1e-30)
        if mask2.sum() > 1:
            ys = y2[mask2]; Ps = P2[mask2]; si = np.argsort(ys)
            ax_p.fill_betweenx(ys[si], 0, Ps[si], alpha=0.75, color='#e74c3c')
            ax_p.plot(Ps[si], ys[si], color='white', lw=1.5)
        # Large disk (faint)
        nf1 = res['node_forces1']; h1 = 2*np.pi*res['R']/120
        P1  = np.abs(nf1[:,0])/h1; y1 = res['mesh1']['perimeter_pos'][:,1]
        mask1 = P1 > 0.01*(P1.max()+1e-30)
        if mask1.sum() > 1:
            ys1 = y1[mask1]; Ps1 = P1[mask1]; si1 = np.argsort(ys1)
            ax_p.fill_betweenx(ys1[si1], 0, Ps1[si1], alpha=0.2, color='#3498db')
            ax_p.plot(Ps1[si1], ys1[si1], color='#3498db', lw=1, alpha=0.5, ls='--')
        pad = R1*0.18
        ax_p.set_ylim(-R1-pad, R1+pad)
        ax_p.set_xlabel('Pressure (Pa)', color='white', fontsize=7)
        ax_p.set_ylabel('y (m)', color='white', fontsize=7)
        ax_p.set_title('Contact pressure\nred=small  blue=large',
                       color='white', fontsize=8, pad=3)
        ax_p.tick_params(colors='white', labelsize=6)
        for sp in ax_p.spines.values(): sp.set_edgecolor('#333')

        fig.tight_layout()

    anim = FuncAnimation(fig, draw_frame, frames=len(all_results), interval=420)
    path = 'results/size_ratio_4to1.gif'
    anim.save(path, writer=PillowWriter(fps=2.5), dpi=110)
    plt.close(fig)
    print(f'  Saved {path}')


# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    t0 = time.time()
    make_gif1()
    make_gif2()
    make_gif3()
    print(f'\nAll GIFs done in {time.time()-t0:.0f}s')
    print('Files in results/')
