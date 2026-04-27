"""
rb_collision_movie.py — Animated GIFs of two-particle EPD collisions.

Generates movies for five representative conditions:

  M1: Equal mass, quasi-elastic  (q=1, α=0,  v0=0.05)
  M2: Equal mass, damped         (q=1, α=2,  v0=0.05)
  M3: Stiff shell (glass-like)   (q=0.1, α=0, v0=0.05)
  M4: Soft shell (rubber-like)   (q=10,  α=0, v0=0.05)
  M5: Mass asymmetry 1:4         (R0_2=2×R0_1, α=2, v0=0.2)

Output: results/rb_benchmarks/movies/collision_M{1..5}.gif
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from io import BytesIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.simulation.capsule_shell import CapsuleParticle, CapsuleSim


# ── Frame renderer ────────────────────────────────────────────────────────────

def render_frame(p1, p2, t, KE, PE, W_diss, label, xlim, ylim):
    """Return a PNG byte buffer for one animation frame."""
    # Portrait figure sized to match the data aspect ratio cleanly
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    fig_w  = 4.0
    fig_h  = fig_w * y_span / x_span + 0.9   # extra for title
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for p, color in [(p1, 'steelblue'), (p2, 'crimson')]:
        # Inner polygon: node positions at radius R0
        xs = np.append(p.x[:, 0], p.x[0, 0])
        ys = np.append(p.x[:, 1], p.x[0, 1])
        ax.fill(xs, ys, alpha=0.30, color=color)
        ax.plot(xs, ys, color=color, lw=1.0)
        # Outer contact shell: nodes expanded radially by r_c
        # This is the true surface where contact forces first activate
        dirs = p.x - p.x_cm          # (N, 2) vectors from CM to each node
        norms = np.linalg.norm(dirs, axis=1, keepdims=True)
        unit = dirs / np.where(norms > 0, norms, 1.0)
        outer = p.x + unit * p.r_c   # push each node outward by r_c
        oxs = np.append(outer[:, 0], outer[0, 0])
        oys = np.append(outer[:, 1], outer[0, 1])
        ax.fill(oxs, oys, alpha=0.10, color=color)
        ax.plot(oxs, oys, color=color, lw=1.5, linestyle='--', alpha=0.7)
        ax.plot(*p.x_cm, 'x', color=color, ms=6, mew=1.8)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=8)
    ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(
        f'{label}\n'
        f't={t:.3f}   KE/KE₀={KE:.3f}   W_d/KE₀={W_diss:.3f}',
        fontsize=7.5, pad=3)
    ax.grid(True, alpha=0.25, lw=0.5)

    buf = BytesIO()
    fig.tight_layout(pad=0.4)
    plt.savefig(buf, format='png', dpi=100)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


# ── Collision runner with frame capture ───────────────────────────────────────

def run_movie(R0_1, R0_2, q, alpha_damp, v0, tau, N, label, out_path,
              every=30, n_post=2000, fps=20):
    """
    every   : capture one frame every this many simulation steps (controls speed)
    n_post  : simulation steps to run after separation (longer = see more fly-apart)
    fps     : playback frame rate
    """
    S = 1.0
    C = 3000.0 * S * (1.0 + q)
    p1 = CapsuleParticle(N=N, R0=R0_1, tau=tau, S=S, C=C, rho_d=1.0)
    p2 = CapsuleParticle(N=N, R0=R0_2, tau=tau, S=S, C=C, rho_d=1.0)

    contact_thr = R0_1 + R0_2 + p1.r_c + p2.r_c
    # Tiny gap so approach is short (~0.5s) and the collision dominates
    start_dist  = contact_thr + 0.04 * max(R0_1, R0_2)

    p2.set_center([0.0,  0.0])
    p1.set_center([0.0, -start_dist])
    p1.set_velocity([0.0,  v0])
    p2.set_velocity([0.0, -v0])

    sim = CapsuleSim([p1, p2], primitives=[])
    dt_max, _ = sim.estimate_dt_max()
    dt = 0.35 * dt_max

    KE_0 = (0.5 * p1.M_disk * np.dot(p1.v_cm, p1.v_cm) +
            0.5 * p2.M_disk * np.dot(p2.v_cm, p2.v_cm))

    # Fixed view that fits both disks throughout entire trajectory.
    # After elastic collision, particles fly ~start_dist apart from the
    # midpoint, so y-range must cover ±(start_dist/2 + max_R + margin).
    r_max   = max(R0_1, R0_2)
    rc_max  = max(p1.r_c, p2.r_c)
    margin  = 0.25 * r_max
    # Horizontal: widest disk + margin
    x_half  = r_max + rc_max + margin
    # Vertical: p1 starts at -start_dist (bottom), p2 at 0; after bounce
    # the faster particle goes up to ~+start_dist.  Centre the view at the
    # mid-point between the two initial positions = -start_dist/2.
    y_centre = -start_dist / 2.0
    y_half   = start_dist / 2.0 + r_max + rc_max + margin
    xlim = (-x_half, x_half)
    ylim = (y_centre - y_half, y_centre + y_half)

    T_max_steps = int((start_dist / v0 * 3.5) / dt) + n_post + 200

    frames  = []
    W_diss  = 0.0
    phase   = 'approach'
    n_sep   = 0
    t       = 0.0

    for step in range(T_max_steps):
        if step % every == 0:
            KE = sim.kinetic_energy_rb() + sim.elastic_kinetic_energy()
            PE = sim.potential_energy()
            frames.append(render_frame(
                p1, p2, t, KE / KE_0, W_diss / KE_0, W_diss / KE_0,
                label, xlim, ylim))

        sim.step_rb(dt, alpha_damp=alpha_damp)
        W_diss += sim.dissipation_rate_rb(alpha_damp) * dt
        t += dt

        dist = np.linalg.norm(p1.x_cm - p2.x_cm)
        gap  = dist - contact_thr

        if phase == 'approach' and gap < 0.05 * max(R0_1, R0_2):
            phase = 'contact'
        elif phase == 'contact' and gap > 0.05 * max(R0_1, R0_2):
            phase = 'separated'
        elif phase == 'separated':
            n_sep += 1
            if n_sep >= n_post:
                break

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    imageio.mimsave(out_path, frames, fps=fps)
    duration = len(frames) / fps
    print(f"  → {out_path}  ({len(frames)} frames, {duration:.1f}s @ {fps}fps)")
    return len(frames)


# ── Movie specs ───────────────────────────────────────────────────────────────

MOVIES = [
    # v0=0.2 for all: 5× faster approach than v0=0.05, deformation still visible,
    # and particles separate noticeably in the post-collision phase.
    dict(label='M1: Equal mass, elastic (α=0)',
         R0_1=1.0, R0_2=1.0, q=1.0, alpha_damp=0.0, v0=0.2, tau=0.2, N=32,
         every=30, out='results/rb_benchmarks/movies/collision_M1_elastic.gif'),

    dict(label='M2: Equal mass, damped (α=2)',
         R0_1=1.0, R0_2=1.0, q=1.0, alpha_damp=2.0, v0=0.2, tau=0.2, N=32,
         every=30, out='results/rb_benchmarks/movies/collision_M2_damped.gif'),

    dict(label='M3: Glass-like shell (q=0.1, ν≈0.21, α=0)',
         R0_1=1.0, R0_2=1.0, q=0.1, alpha_damp=0.0, v0=0.2, tau=0.2, N=32,
         every=30, out='results/rb_benchmarks/movies/collision_M3_glass.gif'),

    dict(label='M4: Rubber-like shell (q=10, ν≈0.88, α=0)',
         R0_1=1.0, R0_2=1.0, q=10.0, alpha_damp=0.0, v0=0.2, tau=0.2, N=32,
         every=30, out='results/rb_benchmarks/movies/collision_M4_rubber.gif'),

    dict(label='M5: Mass ratio 1:4 (R0_2=2, α=2)',
         R0_1=1.0, R0_2=2.0, q=1.0, alpha_damp=2.0, v0=0.2, tau=0.2, N=32,
         every=50, out='results/rb_benchmarks/movies/collision_M5_massratio.gif'),
]


# ── Main ──────────────────────────────────────────────────────────────────────

print("Generating collision movies (EPD rigid-body integrator)")
print("=" * 60)

for spec in MOVIES:
    print(f"\n{spec['label']}")
    run_movie(
        R0_1=spec['R0_1'], R0_2=spec['R0_2'],
        q=spec['q'], alpha_damp=spec['alpha_damp'],
        v0=spec['v0'], tau=spec['tau'], N=spec['N'],
        label=spec['label'], out_path=spec['out'],
        every=spec['every'], n_post=2000, fps=20,
    )

print("\nDone. All movies in results/rb_benchmarks/movies/")
