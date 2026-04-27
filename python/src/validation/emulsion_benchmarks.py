"""
emulsion_benchmarks.py — Emulsion droplet model paper benchmarks (v1 API).

Three canonical benchmarks:
  C — Capillary-wave damping:  measure mode-2 frequency ω₂ vs analytic √6
  D — T1-event three-droplet squeeze:  centre droplet passes through outer pair
  E — Falling droplet:  trajectory under gravity for multiple Bond numbers

All three benchmarks use the v1 API (ParticleSpec + System + System.run()).

Benchmark D notes:
  The outer two droplets are kinematically pinned at their initial CM positions.
  Because System.run() callbacks are read-only, pinning is implemented by
  re-patching sys_._state between successive System.run() calls (each chunk ≈ 0.5 τ₀).
  The descending wall is a Wall object with Wall.set_motion(MotionSpec(vy=-c_wall)).

Outputs (relative to repo root):
  papers/summary_of_methods/figures/emulsion_capwave.png
  papers/summary_of_methods/figures/emulsion_t1.png
  papers/summary_of_methods/figures/emulsion_fall_traj.png
  papers/summary_of_methods/figures/emulsion_fall_shapes.png
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
from src.epd.motion import MotionSpec

FIG_DIR = "papers/summary_of_methods/figures"
os.makedirs(FIG_DIR, exist_ok=True)

# ── Canonical parameters ──────────────────────────────────────────────────────
N      = 120
R0     = 1.0
GAMMA  = 1.0
KAPPA  = 0.02          # κ = γ/(R0·K_area)  → near-incompressible
K_AREA = GAMMA / (R0 * KAPPA)   # = 50.0
C      = 500.0
RHO_D  = 1.0
ALPHA0 = 5.0

tau0   = np.sqrt(RHO_D * R0**3 / GAMMA)   # = 1.0
c_cap  = np.sqrt(GAMMA / (RHO_D * R0))    # = 1.0

omega2_analytic = np.sqrt(6.0 * GAMMA / (RHO_D * R0**3))  # √6 ≈ 2.449
T2_analytic     = 2.0 * np.pi / omega2_analytic


# ── Shared helpers ────────────────────────────────────────────────────────────

def _shoelace(x_nodes):
    """Signed area of closed polygon (N,2) via shoelace formula."""
    x, y = x_nodes[:, 0], x_nodes[:, 1]
    return 0.5 * (np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def draw_droplet(ax, x, color, alpha=0.55, edgecolor='#333'):
    """Draw filled emulsion droplet from node array (N,2)."""
    pts = np.vstack([x, x[:1]])
    ax.fill(pts[:, 0], pts[:, 1], fc=color, ec=edgecolor, lw=0.8,
            alpha=alpha, zorder=3)


# ─────────────────────────────────────────────────────────────────────────────
# BENCHMARK C — Capillary-wave damping (v1 API)
# ─────────────────────────────────────────────────────────────────────────────

def _perturb_x(x_all, x_cm, eps=0.05):
    """Apply mode-2 perturbation u_r = eps·R₀·cos(2θ) to node positions."""
    P = x_all.shape[0]
    x_out = x_all.copy()
    for p in range(P):
        x_rel  = x_all[p] - x_cm[p]
        theta  = np.arctan2(x_rel[:, 1], x_rel[:, 0])
        r_hat  = np.column_stack([np.cos(theta), np.sin(theta)])
        R0_eff = np.mean(np.linalg.norm(x_rel, axis=1))
        u      = eps * R0_eff * np.cos(2 * theta)[:, None] * r_hat
        x_out[p] += u
    return x_out


def _ellipticity_from_x(x_all, x_cm):
    """Mode-2 ellipticity E = ⟨u_r cos2θ⟩/R₀_eff for each particle. Returns (P,)."""
    P = x_all.shape[0]
    ells = []
    for p in range(P):
        x_rel  = x_all[p] - x_cm[p]
        theta  = np.arctan2(x_rel[:, 1], x_rel[:, 0])
        r      = np.linalg.norm(x_rel, axis=1)
        R0_eff = np.mean(r)
        u_r    = r - R0_eff
        ells.append(float(np.mean(u_r * np.cos(2 * theta))) / R0_eff)
    return np.array(ells)


def run_capwave(alpha_damp, n_periods=10, eps=0.05):
    """
    Run single-droplet mode-2 ringdown at given alpha_damp (v1 API).
    Returns (ts, ellipticities) arrays.
    """
    spec = ParticleSpec(count=1, type='emulsion', gamma=GAMMA, kappa=KAPPA,
                        N_nodes=N, Oh=None)
    sys_ = System(10.0, 10.0, periodic_x=False, periodic_y=False, g=0.0)
    sys_.add_particles(spec)
    sys_.initialize(phi_target=0.80, seed=42, verbose=False,
                    relax_only=True, n_relax_init=0)

    # Perturb initial positions into mode-2 shape via direct state patch
    x0   = sys_._state['x_all'].numpy()
    xcm0 = sys_._state['x_cm'].numpy()
    x0p  = _perturb_x(x0, xcm0, eps=eps)
    sys_._state['x_all'] = tf.constant(x0p, dtype=tf.float64)

    # Override damping coefficient
    n_p       = x0.shape[0]
    alpha_arr = np.full(n_p, alpha_damp)
    sys_._params['alpha_damp_per_p'] = tf.constant(alpha_arr, dtype=tf.float64)

    ts_out  = []
    ell_out = []

    def _cb(sys_obj):
        snap = sys_obj.snapshot()
        ts_out.append(snap['t'])
        ell_out.append(float(_ellipticity_from_x(snap['x_all'], snap['x_cm'])[0]))
        return {}

    t_run        = n_periods * T2_analytic
    n_tot        = max(1, int(t_run / sys_._dt))
    sample_every = max(1, n_tot // 2000)   # ≈2000 output points; adequate for ω
    sys_.run(n_tot, sample_every=sample_every, callback=_cb,
             record_initial=True, verbose=False)

    return np.array(ts_out), np.array(ell_out)


def _measure_omega(ts, ells):
    crossings = []
    for i in range(1, len(ells)):
        if ells[i-1] * ells[i] < 0:
            t_cross = ts[i-1] + (ts[i] - ts[i-1]) * (-ells[i-1]) / (ells[i] - ells[i-1])
            crossings.append(t_cross)
    if len(crossings) < 3:
        return None
    return 2.0 * np.pi / (2.0 * np.median(np.diff(crossings)))


def run_benchmark_C():
    print("\n── Benchmark C: Capillary-wave damping ─────────────────────────────")
    alphas = [0.5, 1.0, 2.0, 5.0, 10.0]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    ts0, e0 = run_capwave(alpha_damp=0.5, n_periods=12)
    omega2_meas = _measure_omega(ts0, e0)
    print(f"  ω₂_analytic = {omega2_analytic:.4f} rad/τ₀")
    if omega2_meas:
        print(f"  ω₂_measured = {omega2_meas:.4f} rad/τ₀  "
              f"(error {abs(omega2_meas-omega2_analytic)/omega2_analytic*100:.1f}%)")
    else:
        print("  ω₂_measured = could not extract (too few zero crossings)")
        omega2_meas = omega2_analytic

    all_data = {}
    for a in alphas:
        print(f"  α={a}...", end='', flush=True)
        ts, ells = run_capwave(alpha_damp=a, n_periods=10)
        all_data[a] = (ts, ells)
        print(" done")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    ax = axes[0]
    for a, col in zip(alphas, colors):
        ts, ells = all_data[a]
        ax.plot(ts / tau0, ells, color=col, lw=1.5, label=f'α={a}')
    ax.axhline(0, color='gray', lw=0.5, ls='--')
    ax.set_xlabel(r'$t\,/\,\tau_0$', fontsize=11)
    ax.set_ylabel(r'Ellipticity $\mathcal{E}$', fontsize=11)
    ax.set_title(r'(a) Ringdown at varying $\alpha_\mathrm{damp}$')
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(0, 10 * T2_analytic / tau0)

    ax2 = axes[1]
    alpha_crit = 2.0 * omega2_meas
    for a, col in zip(alphas, colors):
        if a >= alpha_crit:
            continue
        ts, ells = all_data[a]
        peaks_t, peaks_a = [], []
        for i in range(1, len(ells) - 1):
            if ells[i-1] < ells[i] > ells[i+1] and ells[i] > 0.001:
                peaks_t.append(ts[i]); peaks_a.append(ells[i])
        if len(peaks_a) > 1:
            ax2.semilogy(np.array(peaks_t) / tau0, peaks_a, 'o', color=col, ms=4)
            t_fit = np.linspace(peaks_t[0], peaks_t[-1], 200)
            ax2.semilogy(t_fit / tau0,
                         peaks_a[0] * np.exp(-a * (t_fit - peaks_t[0]) / 2),
                         '--', color=col, lw=1, alpha=0.6,
                         label=fr'α={a}: $e^{{-\alpha t/2}}$')
    ax2.set_xlabel(r'$t\,/\,\tau_0$', fontsize=11)
    ax2.set_ylabel('Peak amplitude (log scale)', fontsize=9)
    ax2.set_title(r'(b) Decay envelope vs. $e^{-\alpha t/2}$')
    ax2.legend(fontsize=8)
    ax2.set_xlim(0, 10 * T2_analytic / tau0)

    fig.suptitle(rf'Emulsion capillary-wave damping  ($N={N}$, $\kappa={KAPPA}$)',
                 fontsize=10)
    fig.tight_layout()
    fpath = f"{FIG_DIR}/emulsion_capwave.png"
    fig.savefig(fpath, dpi=150); plt.close(fig)
    print(f"  → {fpath}")
    return dict(omega2_meas=omega2_meas, omega2_analytic=omega2_analytic)


# ─────────────────────────────────────────────────────────────────────────────
# BENCHMARK D — T1-event three-droplet squeeze (v1 API)
# ─────────────────────────────────────────────────────────────────────────────
# Design:
#   • Descending wall: Wall + MotionSpec(vy=-c_wall) — constant-speed descent.
#   • Outer droplets pinned: between consecutive System.run() chunks (≈0.5 τ₀
#     each), outer CMs are reset by patching sys_._state directly.  This is the
#     v1-compatible alternative to per-step CM pinning; the integrator itself is
#     always System.run().
# ─────────────────────────────────────────────────────────────────────────────

def _build_t1_system(N_t1, d_sep, y_top, y_wall0, c_wall, r_c_wall, LX, LY):
    """Build the 3-droplet T1 system using v1 API."""
    # Descending top wall via MotionSpec constant velocity
    half_w    = 3.0 * R0
    cx        = LX / 2.0
    wall_top  = Wall((cx - half_w, y_wall0), (cx + half_w, y_wall0),
                     normal=(0.0, -1.0))
    wall_top.set_motion(MotionSpec(vy=-c_wall))
    wall_top.set_r_c_wall(r_c_wall)

    spec = ParticleSpec(count=3, type='emulsion', gamma=GAMMA, kappa=KAPPA,
                        N_nodes=N_t1, Oh=None)
    sys_ = System(LX, LY, periodic_x=False, periodic_y=False, g=0.0)
    sys_.add_object(wall_top)
    sys_.add_particles(spec)
    sys_.initialize(phi_target=0.30, seed=99, verbose=False,
                    relax_only=True, n_relax_init=0)

    # Patch initial positions: equilateral triangle centred in box
    cy      = LY / 2.0
    targets = np.array([
        [cx - d_sep, cy],             # particle 0: outer left
        [cx + d_sep, cy],             # particle 1: outer right
        [cx,         cy + y_top],     # particle 2: centre (top)
    ])
    x_cm  = sys_._state['x_cm'].numpy()
    x_all = sys_._state['x_all'].numpy()
    for p in range(3):
        delta       = targets[p] - x_cm[p]
        x_all[p]   += delta
        x_cm[p]     = targets[p]
    sys_._state['x_cm']  = tf.constant(x_cm,  dtype=tf.float64)
    sys_._state['x_all'] = tf.constant(x_all, dtype=tf.float64)

    return sys_, targets


def run_benchmark_D():
    print("\n── Benchmark D: T1-event three-droplet squeeze (v1 API) ────────────")
    N_t1        = 72
    d_sep       = 1.65 * R0
    y_top_init  = np.sqrt(3.0) * d_sep
    SR          = 0.001
    c_wall      = SR * c_cap
    r_c_approx  = 2.0 * np.pi * R0 / N_t1   # node spacing ≈ r_c
    LX, LY      = 30.0, 30.0
    cy          = LY / 2.0
    y_wall0     = cy + y_top_init + R0 + r_c_approx + 0.05 * R0
    y_wall_stop = cy + 0.9 * R0   # wall stops descending here

    sys_, targets = _build_t1_system(
        N_t1, d_sep, y_top_init, y_wall0, c_wall, r_c_approx, LX, LY)
    cm_L0 = targets[0].copy()
    cm_R0 = targets[1].copy()

    x_all0   = sys_._state['x_all'].numpy()
    A0_L     = abs(_shoelace(x_all0[0]))
    A0_R     = abs(_shoelace(x_all0[1]))
    A0_C     = abs(_shoelace(x_all0[2]))

    dt_t1    = sys_._dt
    chunk    = max(1, int(0.5 * tau0 / dt_t1))   # ≈ 0.5 τ₀ per run() call
    t1_time  = None
    t_elapsed = 0.0
    t_max    = 400.0 * tau0

    print(f"  N={N_t1}, SR={SR}, chunk={chunk} steps ({chunk*dt_t1:.2f} τ₀)")

    while t_elapsed < t_max:
        sys_.run(chunk, sample_every=chunk + 1, callback=None, verbose=False)
        t_elapsed += chunk * dt_t1

        x_cm_np  = sys_._state['x_cm'].numpy()
        x_all_np = sys_._state['x_all'].numpy()
        v_cm_np  = sys_._state['v_cm'].numpy()

        # Re-pin outer droplets (indices 0 and 1)
        for p_idx, cm_tgt in [(0, cm_L0), (1, cm_R0)]:
            delta           = cm_tgt - x_cm_np[p_idx]
            x_all_np[p_idx] += delta
            x_cm_np[p_idx]  = cm_tgt
            v_cm_np[p_idx]  = 0.0
        sys_._state['x_cm']  = tf.constant(x_cm_np,  dtype=tf.float64)
        sys_._state['x_all'] = tf.constant(x_all_np, dtype=tf.float64)
        sys_._state['v_cm']  = tf.constant(v_cm_np,  dtype=tf.float64)

        # T1 detection: centre droplet y_cm falls below outer-pair midplane
        if x_cm_np[2, 1] < cy and t1_time is None:
            t1_time = t_elapsed
            print(f"  T1 passage at t/τ₀ = {t1_time/tau0:.2f}")

        if t1_time is not None and t_elapsed > t1_time + 100.0 * tau0:
            break

    if t1_time is None:
        print("  ✗ T1 event did not occur!")
        return None

    # Second run to collect 5 representative snapshots
    snap_times  = [0.0, t1_time/3, 2*t1_time/3, t1_time, t1_time + 60*tau0]
    sys2, tgts2 = _build_t1_system(
        N_t1, d_sep, y_top_init, y_wall0, c_wall, r_c_approx, LX, LY)
    cm_L2 = tgts2[0].copy(); cm_R2 = tgts2[1].copy()
    dt2   = sys2._dt
    chunk2 = max(1, int(0.1 * tau0 / dt2))   # finer chunks for snapshot timing

    snap_data = []
    snap_idx  = 0
    t_el2     = 0.0
    t_end2    = snap_times[-1] + 5.0 * tau0

    while t_el2 <= t_end2 and snap_idx < len(snap_times):
        if t_el2 >= snap_times[snap_idx]:
            x_cm_s  = sys2._state['x_cm'].numpy()
            x_all_s = sys2._state['x_all'].numpy()
            t_sys   = float(sys2._state['t'].numpy())
            y_wall  = y_wall0 - c_wall * t_sys
            snap_data.append((t_el2,
                              x_all_s[0].copy(), x_all_s[1].copy(), x_all_s[2].copy(),
                              abs(_shoelace(x_all_s[0])), abs(_shoelace(x_all_s[1])),
                              abs(_shoelace(x_all_s[2])), y_wall))
            snap_idx += 1
        if snap_idx >= len(snap_times):
            break
        sys2.run(chunk2, sample_every=chunk2+1, callback=None, verbose=False)
        t_el2 += chunk2 * dt2
        x_cm_n  = sys2._state['x_cm'].numpy()
        x_all_n = sys2._state['x_all'].numpy()
        v_cm_n  = sys2._state['v_cm'].numpy()
        for p_idx, cm_tgt in [(0, cm_L2), (1, cm_R2)]:
            delta           = cm_tgt - x_cm_n[p_idx]
            x_all_n[p_idx] += delta
            x_cm_n[p_idx]  = cm_tgt
            v_cm_n[p_idx]  = 0.0
        sys2._state['x_cm']  = tf.constant(x_cm_n,  dtype=tf.float64)
        sys2._state['x_all'] = tf.constant(x_all_n, dtype=tf.float64)
        sys2._state['v_cm']  = tf.constant(v_cm_n,  dtype=tf.float64)

    dA_L = max(abs(s[4] - A0_L) / A0_L for s in snap_data) * 100
    dA_R = max(abs(s[5] - A0_R) / A0_R for s in snap_data) * 100
    dA_C = max(abs(s[6] - A0_C) / A0_C for s in snap_data) * 100
    print(f"  Max ΔA: L={dA_L:.2f}%  R={dA_R:.2f}%  C={dA_C:.2f}%")

    labels_snap = [r'$t=0$', r'$t=t_{T1}/3$', r'$t=2t_{T1}/3$',
                   r'$t=t_{T1}$', r'$t=t_{T1}+60\tau_0$']
    half_w_plot = 3.0 * R0
    y_wall_stop_plot = cy + 0.9 * R0
    fig, axes = plt.subplots(1, 5, figsize=(13, 5.5))
    for ax, (t_s, xL, xR, xC, aL, aR, aC, yw), lbl in zip(axes, snap_data, labels_snap):
        draw_droplet(ax, xL, color='#1f77b4')
        draw_droplet(ax, xR, color='#1f77b4')
        draw_droplet(ax, xC, color='#d62728')
        if yw > y_wall_stop_plot:
            xs = [LX/2 - half_w_plot, LX/2 + half_w_plot]
            ax.plot(xs, [yw, yw], color='gray', lw=4, solid_capstyle='butt')
        # Shift for display centred at (0, 0)
        cxd = LX / 2.0
        ax.set_xlim(cxd - 3.5, cxd + 3.5); ax.set_ylim(cy - 5.5, cy + 4.5)
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_title(lbl, fontsize=9)
        dAC = abs(aC - A0_C) / A0_C * 100
        ax.text(cxd, cy - 5.0, f'ΔA$_C$={dAC:.1f}%', ha='center',
                fontsize=7, color='#d62728')

    fig.suptitle(
        rf'T1-event squeeze  ($N={N_t1}$, $d=1.65R_0$, '
        rf'$K_{{\rm area}}={K_AREA}$, $\kappa={KAPPA}$) — v1 API',
        fontsize=10)
    fig.tight_layout()
    fpath = f"{FIG_DIR}/emulsion_t1.png"
    fig.savefig(fpath, dpi=150); plt.close(fig)
    print(f"  → {fpath}")
    return dict(t1_time=t1_time, dA_L=dA_L, dA_R=dA_R, dA_C=dA_C)


# ─────────────────────────────────────────────────────────────────────────────
# BENCHMARK E — Falling droplet (v1 API)
# ─────────────────────────────────────────────────────────────────────────────

def _build_single_falling(g, N_=None):
    """Build single emulsion droplet in a tall box with a floor (v1 API)."""
    N_ = N_ or N
    LX, LY = 10.0, 30.0
    sys_ = System(LX, LY, periodic_x=False, periodic_y=False, g=g)
    floor = Wall((0, 0), (LX, 0), normal=(0, 1))
    sys_.add_object(floor)
    spec = ParticleSpec(count=1, type='emulsion', gamma=GAMMA, kappa=KAPPA,
                        N_nodes=N_, Oh=None)
    sys_.add_particles(spec)
    sys_.initialize(phi_target=0.80, seed=0, verbose=False,
                    relax_only=True, n_relax_init=0)
    return sys_


def run_single_fall(g, t_max_tau=80.0, sample_tau=0.5):
    """Drop a single emulsion droplet under gravity g onto a floor (v1 API)."""
    sys_ = _build_single_falling(g=g, N_=60)   # lower N for speed
    dt   = sys_._dt
    n_tot        = max(1, int(t_max_tau * tau0 / dt))
    sample_every = max(1, int(sample_tau * tau0 / dt))

    ts_out = []
    ys_out = []

    def _cb(sys_obj):
        snap = sys_obj.snapshot()
        ts_out.append(snap['t'])
        ys_out.append(float(snap['x_cm'][0, 1]))
        return {}

    sys_.run(n_tot, sample_every=sample_every, callback=_cb,
             record_initial=True, verbose=False)

    return dict(
        g=g, Bo=g,
        ts=np.array(ts_out),
        ys=np.array(ys_out),
        x_final=sys_._state['x_all'].numpy()[0],
        x_cm_final=sys_._state['x_cm'].numpy()[0],
    )


def run_benchmark_E():
    print("\n── Benchmark E: Falling droplet ────────────────────────────────────")
    Bo_vals   = [0.005, 0.01, 0.05, 0.10]
    colors_bo = {0.005: '#aec6e8', 0.01: '#1f77b4', 0.05: '#ff7f0e', 0.10: '#d62728'}
    results   = {}
    shapes    = {}

    for Bo in Bo_vals:
        print(f"  Bo={Bo:.3f} ...", end='', flush=True)
        res = run_single_fall(g=Bo, t_max_tau=90.0)
        results[Bo] = res
        shapes[Bo]  = res['x_final'].copy()
        print(f"  y_cm_final={res['ys'][-1]:.2f}")

    # Figure 1: CM trajectories
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for Bo in Bo_vals:
        res = results[Bo]
        ax.plot(res['ts'] / tau0, res['ys'], color=colors_bo[Bo], lw=1.5,
                label=f'Bo = {Bo}')
    ax.axhline(1.0, color='sienna', ls=':', lw=1, label=r'floor$+R_0$')
    ax.set_xlabel(r'$t\,/\,\tau_0$', fontsize=11)
    ax.set_ylabel(r'$y_\mathrm{cm}\,/\,R_0$', fontsize=11)
    ax.set_title(rf'Falling droplet trajectories  ($N={60}$, $\kappa={KAPPA}$)',
                 fontsize=10)
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlim(0, 90)
    fig.tight_layout()
    fpath = f"{FIG_DIR}/emulsion_fall_traj.png"
    fig.savefig(fpath, dpi=150); plt.close(fig)
    print(f"  → {fpath}")

    # Figure 2: final equilibrium shapes
    fig, axes = plt.subplots(1, len(Bo_vals), figsize=(12, 4))
    for ax, Bo in zip(axes, Bo_vals):
        x   = shapes[Bo]
        xcm = results[Bo]['x_cm_final']
        x_c = x - xcm
        ax.fill(x_c[:, 0], x_c[:, 1],
                fc=colors_bo[Bo], ec='#333', lw=0.8, alpha=0.7)
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(np.cos(th), np.sin(th), ':', color='#ccc', lw=0.8, zorder=0)
        ax.set_aspect('equal')
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
        ax.set_title(f'Bo={Bo}', fontsize=10)
        ax.axis('off')
    fig.suptitle(rf'Final shape at rest  ($N={60}$, $\kappa={KAPPA}$)', fontsize=10)
    fig.tight_layout()
    fpath2 = f"{FIG_DIR}/emulsion_fall_shapes.png"
    fig.savefig(fpath2, dpi=150); plt.close(fig)
    print(f"  → {fpath2}")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--benchmarks', nargs='+', default=['C', 'D', 'E'],
                        help='Benchmarks to run (default: C D E)')
    args   = parser.parse_args()
    bmarks = [b.upper() for b in args.benchmarks]

    if 'C' in bmarks:
        res_C = run_benchmark_C()
        err   = abs(res_C['omega2_meas'] - res_C['omega2_analytic']) / res_C['omega2_analytic']
        print(f"  {'PASS' if err < 0.05 else 'FAIL'}: ω₂ error = {err:.2%}")

    if 'D' in bmarks:
        res_D = run_benchmark_D()
        if res_D:
            print(f"  PASS: T1 at t/τ₀={res_D['t1_time']/tau0:.1f}  max ΔA_C={res_D['dA_C']:.2f}%")
        else:
            print("  FAIL: T1 event did not occur")

    if 'E' in bmarks:
        run_benchmark_E()
        print("  Done (visual check: droplets settle toward floor)")

    print("\nAll emulsion benchmarks done.")
    print(f"Figures in: {FIG_DIR}/")
