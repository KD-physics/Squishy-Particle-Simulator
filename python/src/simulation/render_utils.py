"""
render_utils.py — Universal capsule rendering helpers.

Usage in any script:
    from src.simulation.render_utils import render_capsule_frame, save_gif

The rendered boundary is the true capsule outer contour:
    Minkowski sum of the N-node polygon with a disk of radius r_c.
For N >= 32 this is visually indistinguishable from the analytic offset
curve because the per-corner arc spans < 12° and is already smooth.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# ── default colours ──────────────────────────────────────────────────────────

PARTICLE_COLORS = [
    '#e74c3c', '#3498db', '#2ecc71', '#f39c12', '#9b59b6',
    '#1abc9c', '#e67e22', '#34495e', '#e91e63', '#00bcd4',
    '#8e44ad', '#27ae60', '#d35400', '#2980b9', '#c0392b',
]


# ── contour geometry ─────────────────────────────────────────────────────────

def capsule_outer_contour(nodes, r_c, n_arc=3):
    """
    True outer contour of a capsule = Minkowski sum of the N-gon polygon
    with a disk of radius r_c.

    At each node B the outward bisector of the two adjacent edge normals
    is computed; the corner arc (spanning the turn angle) is approximated
    with n_arc segments.  For N >= 32 the arc spans < 12° so n_arc=3 gives
    sub-pixel accuracy.

    Parameters
    ----------
    nodes : (N, 2) array   perimeter nodes in order (CW or CCW)
    r_c   : float          capsule contact radius
    n_arc : int            arc sample points per corner (default 3)

    Returns
    -------
    pts : (M, 2) array     closed contour  (M = N * n_arc approximately)
    """
    nodes = np.asarray(nodes, dtype=float)
    N = len(nodes)
    centroid = nodes.mean(axis=0)

    pts = []
    for i in range(N):
        A = nodes[(i - 1) % N]
        B = nodes[i]
        C = nodes[(i + 1) % N]

        e1 = B - A;  l1 = np.linalg.norm(e1)
        e2 = C - B;  l2 = np.linalg.norm(e2)
        if l1 < 1e-14 or l2 < 1e-14:
            continue
        e1 /= l1;  e2 /= l2

        # Candidate outward normals (rotate each edge +90°)
        n1 = np.array([-e1[1], e1[0]])
        n2 = np.array([-e2[1], e2[0]])

        # Ensure they point outward from centroid
        if np.dot(n1, B - centroid) < 0:
            n1 = -n1
        if np.dot(n2, B - centroid) < 0:
            n2 = -n2

        # Arc endpoints on the offset circle at B
        p1 = B + r_c * n1   # end of offset edge (A→B)
        p2 = B + r_c * n2   # start of offset edge (B→C)

        a1 = np.arctan2(p1[1] - B[1], p1[0] - B[0])
        a2 = np.arctan2(p2[1] - B[1], p2[0] - B[0])

        # Always sweep the short way (convex corner → |da| < π)
        da = a2 - a1
        if da >  np.pi: da -= 2 * np.pi
        if da < -np.pi: da += 2 * np.pi

        for j in range(n_arc + 1):
            a = a1 + da * j / n_arc
            pts.append(B + r_c * np.array([np.cos(a), np.sin(a)]))

    return np.array(pts)


# ── frame renderer ────────────────────────────────────────────────────────────

def render_capsule_frame(
    ax, fig,
    x_all,          # (P, N, 2) node positions
    x_cms,          # (P, 2)   centre-of-mass positions
    r_c_per_p,      # (P,)     per-particle capsule radius
    box,            # (x_lo, x_hi, y_bot, y_top) current box extents
    title="",
    colors=None,
    margin=0.3,
    label_particles=True,
    label_fontsize=5.5,
    wall_color='#bdc3c7',
    wall_lw=1.5,
    particle_alpha=0.88,
    particle_edge_color='k',
    particle_edge_lw=0.6,
):
    """
    Render one simulation frame onto *ax*.

    Parameters
    ----------
    ax, fig   : matplotlib axes and figure
    x_all     : (P, N, 2) float  –  node positions
    x_cms     : (P, 2) float     –  particle centroids
    r_c_per_p : (P,) float       –  per-particle capsule radius r_c
    box       : (x_lo, x_hi, y_bot, y_top) –  current wall positions
    title     : str  –  axes title
    colors    : list of colour strings; cycles if shorter than P
    margin    : float  –  whitespace outside the box walls

    Returns
    -------
    img : (H, W, 3) uint8 –  RGB frame
    """
    x_lo, x_hi, y_bot, y_top = box
    P = x_all.shape[0]
    colors = colors or PARTICLE_COLORS

    ax.clear()
    ax.set_xlim(x_lo - margin, x_hi + margin)
    ax.set_ylim(y_bot - margin, y_top + margin)
    ax.set_aspect('equal')
    ax.axis('off')

    # ── wall fill (grey regions outside box) ─────────────────────────────────
    # top wall
    ax.fill_between([x_lo - margin, x_hi + margin],
                    [y_top,          y_top         ],
                    [y_top + margin, y_top + margin],
                    color=wall_color, zorder=1)
    # bottom wall
    ax.fill_between([x_lo - margin, x_hi + margin],
                    [y_bot - margin, y_bot - margin],
                    [y_bot,          y_bot         ],
                    color=wall_color, zorder=1)
    # left wall
    ax.fill_betweenx([y_bot, y_top],
                     [x_lo - margin, x_lo - margin],
                     [x_lo,          x_lo         ],
                     color=wall_color, zorder=1)
    # right wall
    ax.fill_betweenx([y_bot, y_top],
                     [x_hi,          x_hi         ],
                     [x_hi + margin, x_hi + margin],
                     color=wall_color, zorder=1)

    # ── wall outlines ─────────────────────────────────────────────────────────
    for (xs, ys) in [
        ([x_lo, x_hi], [y_bot, y_bot]),
        ([x_lo, x_hi], [y_top, y_top]),
        ([x_lo, x_lo], [y_bot, y_top]),
        ([x_hi, x_hi], [y_bot, y_top]),
    ]:
        ax.plot(xs, ys, 'k-', lw=wall_lw, zorder=2)

    # ── particles ─────────────────────────────────────────────────────────────
    for i in range(P):
        contour = capsule_outer_contour(x_all[i], r_c_per_p[i])
        color   = colors[i % len(colors)]

        # Build a closed matplotlib Path from the contour
        verts = np.vstack([contour, contour[:1]])          # close the loop
        codes = ([Path.MOVETO]
                 + [Path.LINETO] * (len(contour) - 1)
                 + [Path.CLOSEPOLY])
        patch = PathPatch(
            Path(verts, codes),
            facecolor=color,
            edgecolor=particle_edge_color,
            linewidth=particle_edge_lw,
            alpha=particle_alpha,
            zorder=3,
        )
        ax.add_patch(patch)

        if label_particles:
            ax.text(
                x_cms[i, 0], x_cms[i, 1], str(i),
                ha='center', va='center',
                fontsize=label_fontsize, fontweight='bold',
                color='white', zorder=4,
            )

    ax.set_title(title, fontsize=9, pad=3)
    fig.tight_layout(pad=0.3)
    fig.canvas.draw()
    buf = fig.canvas.buffer_rgba()
    img = np.frombuffer(buf, dtype=np.uint8).reshape(
        fig.canvas.get_width_height()[::-1] + (4,))[:, :, :3]
    return img


# ── GIF writer ────────────────────────────────────────────────────────────────

def save_gif(images, path, fps=12, loop=0):
    """Save a list of (H, W, 3) uint8 arrays as an animated GIF."""
    import imageio
    imageio.mimsave(path, images, fps=fps, loop=loop)
    print(f"  Saved {len(images)} frames → {path}  ({len(images)/fps:.1f} s @ {fps} fps)")
