"""
ellipse_collision_movies.py — Animated GIFs of rotating ellipse collisions.

Design philosophy (asymmetric L=0 validation):
  The collision is deliberately ASYMMETRIC — different tilts, different spin speeds,
  off-center approach — so neither p=0 nor L=0 holds trivially by symmetry.
  Any conservation violation is immediately visible.

  Initial conditions constructed so BOTH total linear AND angular momentum = 0:

  Linear momentum = 0:
    Equal mass:   v_B = −v_x  (equal speed, opposite direction)
    Unequal mass: v_B = −(M_A/M_B) * v_x  (slower/faster B to balance)

  Angular momentum = 0:
    L = I_A*ωA + I_B*ωB  +  M_A*(−offset_y * v_x)  +  0  = 0
    → offset_y = (I_A*ωA + I_B*ωB) / (M_A * v_x)
    This gives a nonzero y-offset that partially de-centers the collision.

  Spin direction: ωA < 0 (CW) + theta_A = +45°  →  upper-right tip sweeps toward B
                  ωB > 0 (CCW) + theta_B = −20°  →  upper-left tip sweeps toward A
                  Both leading faces drive into each other. ✓

  Monitoring: Δp/p₀ and ΔL/L_s displayed each frame.
  Since both start at zero, any deviation — however small — is clearly a violation.
  The asymmetry guarantees the test is non-trivial.

Movies:
  EC1: Equal ellipses, moderate asymmetric spin, offset_y ≈ 1.4 R0
  EC2: Equal ellipses, stronger spin, different tilt
  EC3: Different-size ellipses (I ratio 3.3×) — unequal v to zero p
  EC4: Same I, different mass — unequal v to zero p, offset cancels spin L

Output: results/ellipse_benchmarks/movies/ellipse_{EC1..EC4}.gif
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
from src.simulation.capsule_shell import CapsuleSim
from src.simulation.ellipse_particle import EllipseParticle

# ── Parameters ────────────────────────────────────────────────────────────────

TAU = 0.05
S   = 1.0
N   = 32
C   = 3.0 * 12.0 / TAU ** 2 * S   # ≈ 14400

COLORS = ['steelblue', 'crimson']


# ── Helpers ───────────────────────────────────────────────────────────────────

def _set_orientation(p, theta):
    """Rotate particle nodes to initial tilt theta (radians, CCW from +x)."""
    if theta == 0.0:
        return
    p.theta = float(theta)
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    p.x = p.x_cm + (R @ p.X_ref.T).T


def _system_momentum(particles):
    p = np.zeros(2)
    for pi in particles:
        p += pi.M_disk * pi.v_cm
    return p


def _system_angular_momentum(particles):
    """Total L about world origin: I*omega + M*(x_cm × v_cm)."""
    L = 0.0
    for pi in particles:
        L += pi.I_disk * pi.omega
        L += pi.M_disk * float(pi.x_cm[0] * pi.v_cm[1] - pi.x_cm[1] * pi.v_cm[0])
    return L


# ── Frame renderer ────────────────────────────────────────────────────────────

def _draw_ellipse_particle(ax, p, color, label=None):
    xs = np.append(p.x[:, 0], p.x[0, 0])
    ys = np.append(p.x[:, 1], p.x[0, 1])
    ax.fill(xs, ys, alpha=0.28, color=color)
    ax.plot(xs, ys, color=color, lw=1.2)

    dirs  = p.x - p.x_cm
    norms = np.linalg.norm(dirs, axis=1, keepdims=True)
    unit  = dirs / np.where(norms > 0, norms, 1.0)
    outer = p.x + unit * p.r_c
    oxs   = np.append(outer[:, 0], outer[0, 0])
    oys   = np.append(outer[:, 1], outer[0, 1])
    ax.fill(oxs, oys, alpha=0.07, color=color)
    ax.plot(oxs, oys, color=color, lw=1.0, linestyle='--', alpha=0.6)
    ax.plot(*p.x_cm, 'x', color=color, ms=5, mew=1.5)

    # Arrow pointing from CM toward first node (= semi-major axis tip in body frame)
    tip = p.x[0]
    ax.annotate('', xy=tip, xytext=p.x_cm,
                arrowprops=dict(arrowstyle='->', color=color, lw=1.0, alpha=0.75))

    if label:
        ax.text(p.x_cm[0], p.x_cm[1] + p.a + p.r_c + 0.12,
                label, ha='center', va='bottom', fontsize=8,
                color=color, fontweight='bold')


def make_ellipse_frame(pA, pB, t, KE_ratio, dpx, dpy, dL, label, xlim, ylim):
    x_span = xlim[1] - xlim[0]
    y_span = ylim[1] - ylim[0]
    fig_w  = 4.5
    fig_h  = fig_w * y_span / x_span + 1.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    _draw_ellipse_particle(ax, pA, COLORS[0], label='A')
    _draw_ellipse_particle(ax, pB, COLORS[1], label='B')

    ax.set_xlim(*xlim);  ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.set_xlabel('x', fontsize=8);  ax.set_ylabel('y', fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f'{label}\nt={t:.3f}   KE/KE₀={KE_ratio:.3f}', fontsize=7.5, pad=3)
    ax.grid(True, alpha=0.25, lw=0.5)

    handles = [mpatches.Patch(color=c, alpha=0.6, label=lbl)
               for c, lbl in zip(COLORS, ['A', 'B'])]
    ax.legend(handles=handles, fontsize=7, loc='upper right')

    # p₀=0 and L₀=0 by construction — deviations are conservation violations
    cons_txt = (f'Δpx/p_s = {dpx:+.2e}\n'
                f'Δpy/p_s = {dpy:+.2e}\n'
                f'ΔL/L_s  = {dL:+.2e}')
    ax.text(0.02, 0.02, cons_txt, transform=ax.transAxes,
            fontsize=6.2, va='bottom', family='monospace',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.82, lw=0.5))

    buf = BytesIO()
    fig.tight_layout(pad=0.4)
    plt.savefig(buf, format='png', dpi=80)
    plt.close(fig)
    buf.seek(0)
    return imageio.imread(buf)


# ── Movie runner ──────────────────────────────────────────────────────────────

def run_ellipse_movie(a_A, b_A, rho_A,
                      a_B, b_B, rho_B,
                      v_x, omega_A, omega_B,
                      theta_A, theta_B,
                      label='', out_path='out.gif',
                      alpha_damp=0.0,
                      n_post=3000, every=400, fps=25):
    """
    Asymmetric two-ellipse collision with p=0 AND L=0 at t=0.

    Momentum balance:
      v_A = (+v_x, 0),  v_B = (−v_x * M_A/M_B, 0)  →  total p = 0

    Angular momentum balance:
      offset_y = (I_A*omega_A + I_B*omega_B) / (M_A * v_x)
      A placed at (−d/2, offset_y),  B at (+d/2, 0)
      → L_spin + L_orb_A = 0  exactly
    """
    pA = EllipseParticle(a=a_A, b=b_A, N=N, tau=TAU, S=S, C=C, rho_d=rho_A)
    pB = EllipseParticle(a=a_B, b=b_B, N=N, tau=TAU, S=S, C=C, rho_d=rho_B)

    # Speed of B that zeros linear momentum
    v_Bx = v_x * (pA.M_disk / pB.M_disk)

    # y-offset that zeros total angular momentum
    L_spin   = pA.I_disk * omega_A + pB.I_disk * omega_B
    offset_y = L_spin / (pA.M_disk * v_x)

    gap_init = pA.r_c + pB.r_c + 0.50
    d_init   = a_A + a_B + gap_init

    pA.set_center([-d_init / 2.0, offset_y])
    pB.set_center([+d_init / 2.0, 0.0])
    pA.set_velocity([+v_x,  0.0], omega=omega_A)
    pB.set_velocity([-v_Bx, 0.0], omega=omega_B)

    _set_orientation(pA, theta_A)
    _set_orientation(pB, theta_B)

    sim = CapsuleSim([pA, pB], primitives=[])
    dt_max, _ = sim.estimate_dt_max()
    dt = 0.35 * dt_max

    p0      = _system_momentum([pA, pB]).copy()
    L0      = _system_angular_momentum([pA, pB])
    KE0     = sim.kinetic_energy_rb()

    # Scale for normalising residuals (characteristic magnitudes, NOT p0/L0 which are ~0)
    p_scale = pA.M_disk * v_x + pB.M_disk * v_Bx   # total moving mass × speed
    L_scale = max(pA.I_disk * abs(omega_A), pB.I_disk * abs(omega_B),
                  pA.M_disk * v_x * max(a_A, a_B))

    print(f"    p0={np.linalg.norm(p0):.2e} (≈0)  L0={L0:.2e} (≈0)")
    print(f"    offset_y={offset_y:.3f}  v_Bx={v_Bx:.4f}")
    print(f"    ωA={omega_A:.3f}  ωB={omega_B:.3f}  "
          f"I_A={pA.I_disk:.3f}  I_B={pB.I_disk:.3f}  "
          f"M_A={pA.M_disk:.3f}  M_B={pB.M_disk:.3f}")

    # View: span both initial positions + motion range
    r_maxA = a_A + pA.r_c
    r_maxB = a_B + pB.r_c
    margin = 0.7
    x_half = d_init / 2.0 + max(r_maxA, r_maxB) + margin
    y_bot  = offset_y - max(r_maxA, r_maxB) - margin
    y_top  = max(r_maxA, r_maxB) + margin + 0.3
    xlim   = (-x_half, x_half)
    ylim   = (y_bot, y_top)

    T_cross   = d_init / (v_x + v_Bx + 1e-10)
    n_total   = int(T_cross * 10 / dt) + n_post + 2000
    contact_r = pA.r_c + pB.r_c

    frames     = []
    in_contact = False
    n_sep_ok   = 0

    for step in range(n_total):
        if step % every == 0:
            KE  = sim.kinetic_energy_rb()
            pf  = _system_momentum([pA, pB])
            Lf  = _system_angular_momentum([pA, pB])
            dp  = pf - p0
            dpx = float(dp[0]) / p_scale
            dpy = float(dp[1]) / p_scale
            dL  = float((Lf - L0) / L_scale)
            frames.append(make_ellipse_frame(
                pA, pB, sim.t, KE / max(KE0, 1e-30),
                dpx, dpy, dL, label, xlim, ylim))

        sim.step_rb(dt, alpha_damp=alpha_damp)

        d_cc = float(np.linalg.norm(pA.x_cm - pB.x_cm))
        in_c = d_cc < (a_A + a_B + contact_r + 0.1)
        if in_c:
            in_contact = True
        if in_contact and not in_c:
            n_sep_ok += 1
            if n_sep_ok >= n_post:
                break
        elif in_contact and in_c:
            n_sep_ok = 0

    for _ in range(6):
        KE  = sim.kinetic_energy_rb()
        pf  = _system_momentum([pA, pB])
        Lf  = _system_angular_momentum([pA, pB])
        dp  = pf - p0
        dpx = float(dp[0]) / p_scale
        dpy = float(dp[1]) / p_scale
        dL  = float((Lf - L0) / L_scale)
        frames.append(make_ellipse_frame(
            pA, pB, sim.t, KE / max(KE0, 1e-30),
            dpx, dpy, dL, label, xlim, ylim))

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    imageio.mimsave(out_path, frames, fps=fps)
    print(f"  → {out_path}  ({len(frames)} frames, {len(frames)/fps:.1f}s @ {fps}fps)")


# ── Movie specs ───────────────────────────────────────────────────────────────
#
# omega_A < 0 (CW) + theta_A = +45°: upper-right tip sweeps toward B    ✓
# omega_B > 0 (CCW) + theta_B = -20°: upper-left tip sweeps toward A    ✓
# Both leading faces drive into each other ("into each other" rotation)
#
# offset_y = (I_A*wA + I_B*wB) / (M_A*vx)  — computed inside run_ellipse_movie
#
# EC1: offset_y ≈ (2.724*(−0.40) + 2.724*(+0.25)) / (3.770*0.08)
#              = 2.724*(−0.15) / 0.302 = −1.36
# EC2: offset_y ≈ (2.724*(−0.60) + 2.724*(+0.35)) / (3.770*0.10)
#              = 2.724*(−0.25) / 0.377 = −1.81  (larger offset, more glancing)
# EC3: I_A=2.724, I_B=0.819; M_A=3.770, M_B=2.199; v_Bx = (3.770/2.199)*0.07 = 0.120
#      offset_y ≈ (2.724*(−0.25) + 0.819*(+0.35)) / (3.770*0.07) = −0.594/0.264 = −2.25
#      → use smaller omega: wA=-0.15, wB=+0.25
#      offset_y = (2.724*(−0.15) + 0.819*(0.25)) / (3.770*0.07) = (−0.409+0.205)/0.264 = −0.77
# EC4: I_A=I_B=2.724; M_A=3.770, M_B=4.465; v_Bx = (3.770/4.465)*0.08 = 0.0676
#      offset_y ≈ (2.724*(−0.50) + 2.724*(+0.40)) / (3.770*0.08)
#              = 2.724*(−0.10) / 0.302 = −0.90

_a_A4, _b_A4, _rho_A4 = 1.5, 0.8, 1.0
_I_A4  = _rho_A4 * np.pi * _a_A4 * _b_A4 * (_a_A4**2 + _b_A4**2) / 4.0
_a_B4, _b_B4 = 1.2, 1.0
_rho_B4 = _I_A4 / (np.pi * _a_B4 * _b_B4 * (_a_B4**2 + _b_B4**2) / 4.0)

MOVIES_EC = [
    # EC1: Equal ellipses — asymmetric tilts, unequal |ω|, offset_y ≈ −1.4
    dict(label='EC1: Equal (1.5×0.8), θA=45° θB=−20°, ωA=−0.4 ωB=+0.25, L=p=0',
         a_A=1.5, b_A=0.8, rho_A=1.0,
         a_B=1.5, b_B=0.8, rho_B=1.0,
         v_x=0.08, omega_A=-0.40, omega_B=+0.25,
         theta_A=np.pi/4, theta_B=-np.pi/9,       # 45° and −20°
         out='results/ellipse_benchmarks/movies/ellipse_EC1_asym.gif'),

    # EC2: Equal ellipses — stronger spin, shallower tilt, larger offset
    dict(label='EC2: Equal (1.5×0.8), θA=50° θB=−15°, ωA=−0.6 ωB=+0.35, L=p=0',
         a_A=1.5, b_A=0.8, rho_A=1.0,
         a_B=1.5, b_B=0.8, rho_B=1.0,
         v_x=0.10, omega_A=-0.60, omega_B=+0.35,
         theta_A=5*np.pi/18, theta_B=-np.pi/12,   # 50° and −15°
         out='results/ellipse_benchmarks/movies/ellipse_EC2_asym.gif'),

    # EC3: Different size — unequal mass means v_B ≠ v_A for p=0
    dict(label='EC3: Diff size (A:1.5×0.8, B:1.0×0.7), θA=45° θB=−25°, L=p=0',
         a_A=1.5, b_A=0.8, rho_A=1.0,
         a_B=1.0, b_B=0.7, rho_B=1.0,
         v_x=0.07, omega_A=-0.15, omega_B=+0.25,
         theta_A=np.pi/4, theta_B=-5*np.pi/36,    # 45° and −25°
         out='results/ellipse_benchmarks/movies/ellipse_EC3_diffsize.gif'),

    # EC4: Same I, different mass — both v and offset adjusted
    dict(label=f'EC4: Same I, diff mass (ρ_B≈{_rho_B4:.3f}), θA=45° θB=−20°, L=p=0',
         a_A=_a_A4, b_A=_b_A4, rho_A=_rho_A4,
         a_B=_a_B4, b_B=_b_B4, rho_B=_rho_B4,
         v_x=0.08, omega_A=-0.50, omega_B=+0.40,
         theta_A=np.pi/4, theta_B=-np.pi/9,        # 45° and −20°
         out='results/ellipse_benchmarks/movies/ellipse_EC4_sameI.gif'),
]


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import time
    t_start = time.time()
    os.makedirs('results/ellipse_benchmarks/movies', exist_ok=True)

    print("Phase 1G — Ellipse Collision Movies (asymmetric L=p=0 validation)")
    print("=" * 65)
    print("  p=0: v_B = −(M_A/M_B)*v_x  (momentum balance)")
    print("  L=0: offset_y = (I_A*wA + I_B*wB)/(M_A*v_x)  (angular balance)")
    print("  Both conditions non-trivial due to asymmetric tilts/spins/sizes.")
    print(f"  EC4: ρ_B = {_rho_B4:.4f}  (same I as A, different mass)")

    for spec in MOVIES_EC:
        print(f"\n  {spec['label']}")
        run_ellipse_movie(
            a_A=spec['a_A'], b_A=spec['b_A'], rho_A=spec['rho_A'],
            a_B=spec['a_B'], b_B=spec['b_B'], rho_B=spec['rho_B'],
            v_x=spec['v_x'],
            omega_A=spec['omega_A'], omega_B=spec['omega_B'],
            theta_A=spec['theta_A'], theta_B=spec['theta_B'],
            label=spec['label'], out_path=spec['out'],
            every=400, fps=25, n_post=3000,
        )

    elapsed = time.time() - t_start
    print(f"\nDone. ({elapsed:.0f}s)")
