"""
rb_extension_movies.py — Animated GIFs for Phase 1G extension benchmarks.

Movies generated:
  3-body collision:
    3B1: q=1,   α=0  (elastic,  equal mass)
    3B2: q=1,   α=2  (damped,   equal mass)
    3B3: q=10,  α=0  (rubber,   equal mass)

  Rearrangement squeeze (T1-event analogue):
    SQ1: q=1,   α=2, squeeze_gap=0.20 (nominal, N=32)
    SQ2: q=0.1, α=2, squeeze_gap=0.20 (glass-like, N=32)
    SQ3: q=10,  α=2, squeeze_gap=0.20 (rubber-like, N=32)

Output: results/rb_benchmarks/movies/extension_{3B1..3B3,SQ1..SQ3}.gif
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import imageio.v2 as imageio
from io import BytesIO

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.simulation.capsule_shell import CapsuleParticle, CapsuleSim
from src.simulation.contact_primitives import LineSegment


# ── Shared particle renderer ──────────────────────────────────────────────────

COLORS = ['steelblue', 'crimson', 'seagreen']

def _draw_particle(ax, p, color):
    """Draw inner polygon + outer contact shell of one particle."""
    xs = np.append(p.x[:, 0], p.x[0, 0])
    ys = np.append(p.x[:, 1], p.x[0, 1])
    ax.fill(xs, ys, alpha=0.30, color=color)
    ax.plot(xs, ys, color=color, lw=1.0)

    dirs  = p.x - p.x_cm
    norms = np.linalg.norm(dirs, axis=1, keepdims=True)
    unit  = dirs / np.where(norms > 0, norms, 1.0)
    outer = p.x + unit * p.r_c
    oxs   = np.append(outer[:, 0], outer[0, 0])
    oys   = np.append(outer[:, 1], outer[0, 1])
    ax.fill(oxs, oys, alpha=0.08, color=color)
    ax.plot(oxs, oys, color=color, lw=1.2, linestyle='--', alpha=0.7)
    ax.plot(*p.x_cm, 'x', color=color, ms=5, mew=1.5)


# ── 3-body movie ──────────────────────────────────────────────────────────────

def make_3body_frame(pA, pB, pC, t, KE_ratio, Wd_ratio, dpx, dpy, label, xlim, ylim):
    """
    dpx, dpy: momentum residuals Σpx - Σpx₀ and Σpy - Σpy₀, normalised by p_scale.
    Displayed as text so the viewer can see machine-precision conservation.
    """
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    fig_w  = 4.2
    fig_h  = fig_w * y_span / x_span + 1.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for p, col in zip([pA, pB, pC], COLORS):
        _draw_particle(ax, p, col)

    ax.set_xlim(*xlim);  ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=8);  ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f'{label}\nt={t:.3f}  KE/KE₀={KE_ratio:.3f}  W_d/KE₀={Wd_ratio:.3f}',
                 fontsize=7.5, pad=3)
    ax.grid(True, alpha=0.25, lw=0.5)

    handles = [mpatches.Patch(color=c, label=lbl, alpha=0.6)
               for c, lbl in zip(COLORS, ['A', 'B', 'C'])]
    ax.legend(handles=handles, fontsize=7, loc='upper right')

    # Momentum conservation display (machine precision)
    mom_txt = (f'Δpx/p₀ = {dpx:+.2e}\n'
               f'Δpy/p₀ = {dpy:+.2e}')
    ax.text(0.02, 0.02, mom_txt, transform=ax.transAxes,
            fontsize=6.5, va='bottom', family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.75, lw=0.5))

    buf = BytesIO()
    fig.tight_layout(pad=0.4)
    plt.savefig(buf, format='png', dpi=90)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


def run_3body_movie(q, alpha_damp, v_down, v_up, b, x_off,
                    label, out_path, every=20, n_post=100, fps=16):
    """
    Runs a 3-body collision movie with generous approach and post-collision coast.

    Changes vs original:
    - gap_init = 1.5 * contact_thr  (longer approach, ~3× more pre-collision frames)
    - n_post increased to 100 steps  (longer post-collision coast)
    - every = 20  (smoother animation)
    - fps = 16    (slightly slower, better for talks)
    - Adds per-frame Δpx/p₀ and Δpy/p₀ momentum residuals (machine precision)
    """
    S  = 1.0
    C  = 3000.0 * S * (1.0 + q)
    N  = 32
    R0 = 1.0

    pA = CapsuleParticle(N=N, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)
    pB = CapsuleParticle(N=N, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)
    pC = CapsuleParticle(N=N, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)

    contact_thr = 2.0 * R0 + pA.r_c + pC.r_c
    V = v_down + v_up

    # Start half as far away; post-collision coast uses doubled n_post for balance
    d_BC_init = contact_thr + 0.75 * contact_thr
    dx_BC     = b - x_off
    y_B0      = np.sqrt(max(d_BC_init**2 - dx_BC**2, 1e-10))

    y_contact_B = np.sqrt(max(contact_thr**2 - dx_BC**2, 0.0))
    t_star      = (y_B0 - y_contact_B) / V

    dx_AC       = b + x_off
    y_contact_A = np.sqrt(max(contact_thr**2 - dx_AC**2, 0.0))
    y_A0        = V * t_star + y_contact_A

    pA.set_center([-b,    y_A0])
    pB.set_center([+b,    y_B0])
    pC.set_center([x_off, 0.0 ])
    pA.set_velocity([0.0, -v_down])
    pB.set_velocity([0.0, -v_down])
    pC.set_velocity([0.0, +v_up  ])

    sim   = CapsuleSim([pA, pB, pC], primitives=[])
    dt_max, _ = sim.estimate_dt_max()
    dt    = 0.35 * dt_max

    KE0   = sim.kinetic_energy_rb()
    M_tot = pA.M_disk + pB.M_disk + pC.M_disk
    p0    = (pA.M_disk * pA.v_cm + pB.M_disk * pB.v_cm + pC.M_disk * pC.v_cm).copy()
    p_scale = float(np.linalg.norm(p0)) if np.linalg.norm(p0) > 1e-30 else (
        M_tot * max(v_down, v_up))

    # View
    r_max  = R0 + pA.r_c
    margin = 0.6 * R0
    x_half = b + r_max + margin
    y_top  = max(y_A0, y_B0) + r_max + margin
    y_bot  = -(y_top * 0.6 + margin)
    xlim   = (-x_half, x_half)
    ylim   = (y_bot, y_top)

    T_safety = int((y_top * 8 / min(v_down, v_up)) / dt) + n_post + 500

    frames      = []
    W_diss      = 0.0
    in_contact  = False
    n_sep       = 0

    for step in range(T_safety):
        if step % every == 0:
            KE = sim.kinetic_energy_rb() + sim.elastic_kinetic_energy()
            p_now = (pA.M_disk * pA.v_cm + pB.M_disk * pB.v_cm
                     + pC.M_disk * pC.v_cm)
            dp    = p_now - p0
            dpx   = float(dp[0]) / p_scale
            dpy   = float(dp[1]) / p_scale
            frames.append(make_3body_frame(
                pA, pB, pC, sim.t, KE / KE0, W_diss / KE0,
                dpx, dpy, label, xlim, ylim))

        sim.step_rb(dt, alpha_damp=alpha_damp)
        W_diss += sim.dissipation_rate_rb(alpha_damp) * dt

        dAC = np.linalg.norm(pA.x_cm - pC.x_cm) - contact_thr
        dBC = np.linalg.norm(pB.x_cm - pC.x_cm) - contact_thr
        dAB = np.linalg.norm(pA.x_cm - pB.x_cm) - contact_thr
        any_contact = (dAC < 0) or (dBC < 0) or (dAB < 0)
        all_sep     = (dAC > 0.05) and (dBC > 0.05) and (dAB > 0.05)

        if any_contact:
            in_contact = True
        if in_contact and all_sep:
            n_sep += 1
            if n_sep >= n_post:
                break

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    imageio.mimsave(out_path, frames, fps=fps)
    print(f"  → {out_path}  ({len(frames)} frames, {len(frames)/fps:.1f}s @ {fps}fps)")


# ── Squeeze movie ──────────────────────────────────────────────────────────────

def make_squeeze_frame(pC, pL, pR, wall_y, wall_hw,
                       outer_x, t, KE, PE_elastic, W_diss, KE0, label, xlim, ylim):
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    fig_w  = 4.5
    fig_h  = fig_w * y_span / x_span + 1.0
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # ── Outer cradling walls (vertical gray lines) ────────────────────────────
    wall_ht = pL.R0 + pL.r_c + 0.15
    for sign, label_w in [(-1, 'cradle L'), (+1, 'cradle R')]:
        wx = sign * outer_x
        ax.plot([wx, wx], [-wall_ht, +wall_ht], color='dimgray', lw=3.0, alpha=0.7,
                solid_capstyle='butt')
        # Hatch marks to indicate fixed wall
        for y_hatch in np.linspace(-wall_ht + 0.1, wall_ht - 0.1, 6):
            ax.plot([wx, wx + sign * 0.15], [y_hatch, y_hatch + 0.12],
                    color='dimgray', lw=1.0, alpha=0.5)

    # ── Side particles (cradled, deformable) ─────────────────────────────────
    _draw_particle(ax, pL, 'slategray')
    _draw_particle(ax, pR, 'slategray')

    # ── Moving wall (brown, below pC) ─────────────────────────────────────────
    ax.plot([-wall_hw, +wall_hw], [wall_y, wall_y],
            color='saddlebrown', lw=3.0, alpha=0.9, solid_capstyle='butt')
    ax.fill_between([-wall_hw, +wall_hw], wall_y - 0.10, wall_y,
                    color='saddlebrown', alpha=0.20)
    # Hatch marks on wall
    for x_hatch in np.linspace(-wall_hw + 0.1, wall_hw - 0.1, 8):
        ax.plot([x_hatch, x_hatch - 0.10], [wall_y, wall_y - 0.12],
                color='saddlebrown', lw=0.8, alpha=0.5)

    # ── Center particle ───────────────────────────────────────────────────────
    _draw_particle(ax, pC, 'steelblue')

    # ── Dashed gap guides at x = ±d_side ────────────────────────────────────
    d_side = float(pL.x_cm[0]) * -1.0   # pL is at -d_side
    ax.axvline(-d_side, color='slategray', lw=0.7, ls=':', alpha=0.35)
    ax.axvline(+d_side, color='slategray', lw=0.7, ls=':', alpha=0.35)

    handles = [
        mpatches.Patch(color='steelblue', alpha=0.6, label='pC (center)'),
        mpatches.Patch(color='slategray', alpha=0.6, label='pL / pR (cradled)'),
    ]
    ax.legend(handles=handles, fontsize=7, loc='upper right')

    ax.set_xlim(*xlim);  ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=8);  ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)
    KE_rat = KE / max(KE0, 1e-12)
    PE_rat = PE_elastic / max(KE0, 1e-12)
    ax.set_title(f'{label}\nt={t:.3f}  KE/KE₀={KE_rat:.2f}  ΔPE/KE₀={PE_rat:.1f}',
                 fontsize=7.5, pad=3)
    ax.grid(True, alpha=0.20, lw=0.5)

    buf = BytesIO()
    fig.tight_layout(pad=0.4)
    plt.savefig(buf, format='png', dpi=90)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


def run_squeeze_movie(q, alpha_damp, v0, A_osc, omega_osc, squeeze_gap,
                      label, out_path, every=15, fps=18):
    S    = 1.0
    R0   = 1.0
    C    = 3000.0 * S * (1.0 + q)
    N_SQ = 32

    # Three EPD particles
    pC = CapsuleParticle(N=N_SQ, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)
    pL = CapsuleParticle(N=N_SQ, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)
    pR = CapsuleParticle(N=N_SQ, R0=R0, tau=0.2, S=S, C=C, rho_d=1.0)

    contact_thr_pp = 2.0 * R0 + 2.0 * pC.r_c
    d_side         = contact_thr_pp - squeeze_gap

    pL.set_center([-d_side, 0.0]);  pL.set_velocity([0.0, 0.0])
    pR.set_center([+d_side, 0.0]);  pR.set_velocity([0.0, 0.0])

    y_start = -1.0
    pC.set_center([0.0, y_start])
    pC.set_velocity([0.0, v0])

    # Outer cradling walls
    outer_x = d_side + R0 + pC.r_c + 0.05
    wall_ht = R0 + pC.r_c + 0.1
    wall_outer_L = LineSegment(
        p0=[-outer_x, -wall_ht], p1=[-outer_x, +wall_ht],
        inward_normal=[+1.0, 0.0])
    wall_outer_R = LineSegment(
        p0=[+outer_x, -wall_ht], p1=[+outer_x, +wall_ht],
        inward_normal=[-1.0, 0.0])

    # Narrow moving wall below pC only
    wall_y  = y_start - R0 - (pC.r_c + 0.01)
    wall_hw = d_side - R0 - pC.r_c - 0.10
    wall    = LineSegment(p0=[-wall_hw, wall_y], p1=[+wall_hw, wall_y],
                          inward_normal=[0.0, 1.0])

    sim   = CapsuleSim([pC, pL, pR], primitives=[wall_outer_L, wall_outer_R, wall])
    dt_max, _ = sim.estimate_dt_max()
    dt    = 0.35 * dt_max

    KE0 = sim.kinetic_energy_rb()
    PE0 = sim.potential_energy()

    # View: covers center action + side particles + outer walls
    x_pad = outer_x + 0.3
    y_bot = wall_y - 0.3
    y_top = 3.0 * R0
    xlim  = (-x_pad, x_pad)
    ylim  = (y_bot, y_top)

    T_pass  = 3.0 / v0 * 1.5
    T_osc   = 3 * 2.0 * np.pi / omega_osc
    T_total = max(T_pass, T_osc) + 1.0
    n_total = int(T_total / dt) + 200

    frames = []
    W_diss = 0.0
    passed = False
    n_post_steps = 0
    # Run ~2 R0/v0 worth of steps after pC clears the gap as post-pass coast
    n_post_target = int(2.0 * R0 / v0 / dt)
    wall_frozen = False

    for step in range(n_total):
        t = sim.t

        # Push wall only while pC is still below the gap; freeze once it clears
        if not wall_frozen:
            v_wall = v0 + A_osc * np.sin(omega_osc * t)
            wall.p0[1] += v_wall * dt
            wall.p1[1] += v_wall * dt
            if pC.x_cm[1] > 0.0:
                wall_frozen = True

        if step % every == 0:
            KE = sim.kinetic_energy_rb()
            PE = sim.potential_energy() - PE0
            frames.append(make_squeeze_frame(
                pC, pL, pR, wall.p0[1], wall_hw, outer_x,
                t, KE, abs(PE), W_diss, KE0, label, xlim, ylim))

        sim.step_rb(dt, alpha_damp=alpha_damp)

        # Pin side particles (rigid-body DOFs only; elastic shape can deform)
        pL.v_cm[:] = 0.0;  pL.omega = 0.0
        pR.v_cm[:] = 0.0;  pR.omega = 0.0

        W_diss += sim.dissipation_rate_rb(alpha_damp) * dt

        if pC.x_cm[1] > 2.0 * R0 and not passed:
            passed = True
        if passed:
            n_post_steps += 1
            if n_post_steps >= n_post_target:
                break

    # Freeze frames after exit
    for _ in range(12):
        KE = sim.kinetic_energy_rb()
        PE = sim.potential_energy() - PE0
        frames.append(make_squeeze_frame(
            pC, pL, pR, wall.p0[1], wall_hw, outer_x,
            sim.t, KE, abs(PE), W_diss, KE0, label, xlim, ylim))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    imageio.mimsave(out_path, frames, fps=fps)
    print(f"  → {out_path}  ({len(frames)} frames, {len(frames)/fps:.1f}s @ {fps}fps)")


# ── Movie specs ────────────────────────────────────────────────────────────────

MOVIES_3B = [
    dict(label='3B1: q=1, α=0 (elastic)',
         q=1.0, alpha_damp=0.0, v_down=0.1, v_up=0.1,
         b=1.5, x_off=0.5, every=600,
         out='results/rb_benchmarks/movies/extension_3B1_elastic.gif'),

    dict(label='3B2: q=1, α=2 (damped)',
         q=1.0, alpha_damp=2.0, v_down=0.1, v_up=0.1,
         b=1.5, x_off=0.5, every=600,
         out='results/rb_benchmarks/movies/extension_3B2_damped.gif'),

    dict(label='3B3: q=10, α=0 (rubber, elastic)',
         q=10.0, alpha_damp=0.0, v_down=0.1, v_up=0.1,
         b=1.5, x_off=0.5, every=600,
         out='results/rb_benchmarks/movies/extension_3B3_rubber.gif'),
]

MOVIES_SQ = [
    dict(label='SQ1: q=1, α=2 (nominal gap=0.20)',
         q=1.0, alpha_damp=2.0, v0=0.20, A_osc=0.04, omega_osc=5.0,
         squeeze_gap=0.20, every=240,
         out='results/rb_benchmarks/movies/extension_SQ1_q1.gif'),

    dict(label='SQ2: q=0.1, α=2 (glass, gap=0.20)',
         q=0.1, alpha_damp=2.0, v0=0.20, A_osc=0.04, omega_osc=5.0,
         squeeze_gap=0.20, every=240,
         out='results/rb_benchmarks/movies/extension_SQ2_glass.gif'),

    dict(label='SQ3: q=10, α=2 (rubber, gap=0.20)',
         q=10.0, alpha_damp=2.0, v0=0.20, A_osc=0.04, omega_osc=5.0,
         squeeze_gap=0.20, every=240,
         out='results/rb_benchmarks/movies/extension_SQ3_rubber.gif'),
]


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    t_start = time.time()
    os.makedirs('results/rb_benchmarks/movies', exist_ok=True)

    print("Phase 1G — Extension Benchmark Movies")
    print("=" * 60)

    print("\n── 3-body collision movies ──────────────────────────────────")
    for spec in MOVIES_3B:
        print(f"\n  {spec['label']}")
        run_3body_movie(
            q=spec['q'], alpha_damp=spec['alpha_damp'],
            v_down=spec['v_down'], v_up=spec['v_up'],
            b=spec['b'], x_off=spec['x_off'],
            label=spec['label'], out_path=spec['out'],
            every=spec['every'], n_post=spec['every']*40, fps=25,
        )

    print("\n── Squeeze movies ────────────────────────────────────────────")
    for spec in MOVIES_SQ:
        print(f"\n  {spec['label']}")
        run_squeeze_movie(
            q=spec['q'], alpha_damp=spec['alpha_damp'],
            v0=spec['v0'], A_osc=spec['A_osc'],
            omega_osc=spec['omega_osc'],
            squeeze_gap=spec['squeeze_gap'],
            label=spec['label'], out_path=spec['out'],
            every=spec['every'], fps=25,
        )

    elapsed = time.time() - t_start
    print(f"\nDone. All movies in results/rb_benchmarks/movies/  ({elapsed:.0f}s)")
