"""
Render the hopper test bed for the paper.

Loads the 50k checkpoint (mid-flow steady state) and produces a static PNG
showing the actual hopper wall geometry (funnel + outlet guides + corner
arcs, read directly from the HopperRegion objects in the simulation) with
particles drawn at their true capsule outer contour (Minkowski sum with
disk of radius r_c — same as the rest of the codebase via render_utils).

Output: papers/summary_of_methods/figures/hopper_testbed.png
"""
import os, sys
sys.path.insert(0, os.getcwd())
import argparse
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Arc as MplArc
from matplotlib.path import Path

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.validation.profile_hopper import build_test_bed
from src.simulation.render_utils import capsule_outer_contour, PARTICLE_COLORS

OUT_DIR = pathlib.Path("results/profiling")
FIG_DIR = pathlib.Path("papers/summary_of_methods/figures")
FIG_DIR.mkdir(parents=True, exist_ok=True)


def _draw_hopper_walls(ax, sys_h, lw=1.6, color='#222222'):
    """Draw line-segment + arc walls from every object in sys_h._objects.

    Reads the geometry directly from each object's ``_walls`` (list of Wall)
    and ``_arc_walls`` (list of ArcWall) attributes — the same primitives
    the simulator sees, so what is drawn is exactly what is colliding.
    """
    for obj in sys_h._objects:
        for w in getattr(obj, '_walls', []):
            p0, p1 = w._p0, w._p1
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
                    color=color, lw=lw, solid_capstyle='round', zorder=2)
        for aw in getattr(obj, '_arc_walls', []):
            c    = aw._center
            r    = aw._radius
            arng = aw._angle_range
            if arng is None:
                a0, a1 = 0.0, 2 * np.pi
            else:
                a0, a1 = arng
            ts = np.linspace(a0, a1, 32)
            xs = c[0] + r * np.cos(ts)
            ys = c[1] + r * np.sin(ts)
            ax.plot(xs, ys, color=color, lw=lw, zorder=2)


def render(P=300, tag='emulsion_Bo005', Bo=0.05, kappa=0.02, Oh=0.15,
           ptype='emulsion', nu=0.5, output='hopper_testbed.png',
           crop_y_lo=None, crop_y_hi=None, dpi=150):
    sys_h, info = build_test_bed(P=P, N_nodes=60, scale=1.0, verbose=False,
                                 Bo=Bo, kappa=kappa, Oh=Oh, ptype=ptype, nu=nu,
                                 E_candidates=64)
    state_path = OUT_DIR / f"hopper_N{P}_{tag}_state_50k.npz"
    sys_h.restore_state(state_path)

    x_all     = sys_h._state['x_all'].numpy()       # (P, N, 2)
    x_cm      = sys_h._state['x_cm'].numpy()        # (P, 2)
    r_c_per_p = np.array([p.r_c for p in sys_h._particles])

    LX = float(sys_h.Lx)

    # Default crop: focus on funnel + outlet + lower reservoir
    Y_BOT        = info['Y_BOT']
    Y_FUNNEL_TOP = Y_BOT + info['funnel_h']
    if crop_y_lo is None:
        crop_y_lo = -2.0
    if crop_y_hi is None:
        crop_y_hi = Y_FUNNEL_TOP + min(40.0, info['H_RES'])

    fig, ax = plt.subplots(figsize=(6, 9))
    ax.set_aspect('equal')
    ax.set_xlim(-1.0, LX + 1.0)
    ax.set_ylim(crop_y_lo, crop_y_hi)
    ax.set_xlabel(r'$x / R_0$', fontsize=11)
    ax.set_ylabel(r'$y / R_0$', fontsize=11)

    # ── Walls (read straight from HopperRegion._walls/_arc_walls) ────────────
    _draw_hopper_walls(ax, sys_h)
    ax.axhline(0, color='#888888', lw=0.8, ls=':', zorder=1)

    # ── Particles: true capsule outer contour, palette-cycled colors ─────────
    for i in range(P):
        ymin, ymax = x_all[i, :, 1].min(), x_all[i, :, 1].max()
        if ymax < crop_y_lo - 2 or ymin > crop_y_hi + 2:
            continue
        contour = capsule_outer_contour(x_all[i], r_c_per_p[i])
        verts   = np.vstack([contour, contour[:1]])
        codes   = ([Path.MOVETO]
                   + [Path.LINETO] * (len(contour) - 1)
                   + [Path.CLOSEPOLY])
        color   = PARTICLE_COLORS[i % len(PARTICLE_COLORS)]
        ax.add_patch(PathPatch(
            Path(verts, codes),
            facecolor=color, edgecolor='k',
            linewidth=0.4, alpha=0.88, zorder=3,
        ))

    title = (rf'Hopper test bed — {ptype}'
             + (rf' ($\kappa$={kappa}, Oh={Oh}, Bo={Bo})' if ptype == 'emulsion'
                else rf' ($\nu$={nu}, Bo={Bo})')
             + f', P={P}, mid-flow @ 50k steps')
    ax.set_title(title, fontsize=10)
    ax.text(0.02, 0.98, f'crop y∈[{crop_y_lo:.1f}, {crop_y_hi:.1f}]   LX={LX:.1f}',
            transform=ax.transAxes, va='top', fontsize=8, alpha=0.6)

    fig.tight_layout()
    out_path = FIG_DIR / output
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"saved {out_path}")
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p',     type=int, default=300)
    ap.add_argument('--tag',   type=str, default='emulsion_Bo005')
    ap.add_argument('--Bo',    type=float, default=0.05)
    ap.add_argument('--kappa', type=float, default=0.02)
    ap.add_argument('--Oh',    type=float, default=0.15)
    ap.add_argument('--ptype', type=str,   default='emulsion')
    ap.add_argument('--nu',    type=float, default=0.5)
    ap.add_argument('--output', type=str, default='hopper_testbed.png')
    args = ap.parse_args()

    render(P=args.p, tag=args.tag, Bo=args.Bo, kappa=args.kappa, Oh=args.Oh,
           ptype=args.ptype, nu=args.nu, output=args.output)


if __name__ == "__main__":
    main()
