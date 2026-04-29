"""
shell_mode_calibration.py — α_damp calibration from bending-mode theory.

Corrected theory
----------------
The EPD shell has two force contributions that give distinct mode families:

  (a) Membrane/extensional modes  (edge spring El_t):
          ω_ext ≈ 2·c_edge/R0  ≈ 78 rad/s at reference
      These are fast in-plane stretch modes.

  (b) Inextensional bending modes  (bending stiffness EI = S·K_fluid·R0³):
          ω_2_bend = (6/√5) · √(EI / (ρ_L · R0⁴))
                   = (6/√5) · √(S / (ρ_f · τ · R0²))
      At reference (S=1, τ=0.20, R0=1, ρ=1):  ω_2_bend ≈ 6.0 rad/s.

The membrane modes are fast and non-resonant during contact (T_contact >> T_ext),
so they do not persist.  The bending modes are ~13× slower and are strongly
excited during contact (the disk squishes elliptically), so they dominate
the visible post-collision ringing.

Decay law (universal — same for ALL modes):
    T_decay = 2 / α

Quality factor for bending mode:
    Q_bend = ω_2_bend / α

Design rules
------------
  Q_bend < 1  (no visible ringing, overdamped):  α > ω_2_bend
  Q_bend < 3  (few cycles, OK for slow dynamics): α > ω_2_bend / 3

At reference params (S=1, τ=0.20, R0=1, N=32, ρ=1):
    ω_2_bend ≈ 6.0 rad/s
    α = 2.0  →  Q ≈ 3    (current default: ~3 visible cycles)
    α = 6.0  →  Q = 1    (threshold for no visible ringing)

General formula for any params:
    α_target = ω_2_bend / Q_target
             = (6/√5) · √(S / (ρ_f · τ)) / (Q_target · R0)

Tests run here
--------------
  1. Free-vibration with correct n=2 bending initial condition → ω_2_bend_meas
  2. T_decay = 2/α universality sweep
  3. ω_2_bend scaling: vary τ, R0, S independently — verify (6/√5)·√(S/(ρ_f·τ))/R0
  4. COR vs α sweep at reference params
  5. Post-collision ringing: direct measurement from head-on collision

Expected precision:
  - ω_2_bend:  5–10% error acceptable (discrete vs continuum corrections at N=32)
  - T_decay:   <1% error (analytically exact)
  - COR:       physically meaningful relative variation across α sweep
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import json, time
from scipy.optimize import curve_fit

from src.simulation.capsule_shell import CapsuleParticle, CapsuleSim


# ── theory formulas ───────────────────────────────────────────────────────────

def omega2_bend_theory(tau, R0, S=1.0, rho_d=1.0, K_fluid=1.0):
    """
    n=2 inextensional bending mode angular frequency for an elastic ring beam.

    Formula (from Euler-Bernoulli ring mechanics):
        ω_n = n(n²-1)/√(n²+1) · √(EI / (ρ_L · R0⁴))

    For n=2:  n(n²-1)/√(n²+1) = 2·3/√5 = 6/√5 ≈ 2.683

    With:
        EI   = S · K_fluid · R0³  (bending stiffness from capsule_shell.py line 136)
        ρ_L  = ρ_f · τ · R0       (linear mass density = m_node / L0)

    Simplifies to:
        ω_2_bend = (6/√5) · √(S / (ρ_f · τ · R0²))
    """
    EI    = S * K_fluid * R0**3
    rho_L = rho_d * tau * R0
    return (6.0 / np.sqrt(5.0)) * np.sqrt(EI / (rho_L * R0**4))


def make_particle(tau=0.20, R0=1.0, N=32, S=1.0, C=3000.0, rho_d=1.0):
    return CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, rho_d=rho_d,
                           center=(0.0, 0.0))


def get_dt(p, dt_factor=0.35):
    sim = CapsuleSim([p])
    dt_max, _ = sim.estimate_dt_max()
    return dt_factor * dt_max


def bending_mode_init(N, R0, A_init):
    """
    Correct n=2 inextensional bending mode initial condition (body frame).

    For inextensional ring bending, mode n=2:
        u_r(θ) = A · cos(2θ)
        u_θ(θ) = (A/2) · sin(2θ)   ← from inextensional condition du_θ/dθ = -u_r

    In Cartesian (rotating frame, θ from x-axis):
        u_x = u_r·cosθ − u_θ·sinθ
        u_y = u_r·sinθ + u_θ·cosθ
    """
    theta = 2.0 * np.pi * np.arange(N) / N
    u_r = A_init * np.cos(2.0 * theta)
    u_t = (A_init / 2.0) * np.sin(2.0 * theta)   # inextensional condition

    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    u_x = u_r * cos_t - u_t * sin_t
    u_y = u_r * sin_t + u_t * cos_t
    return np.column_stack([u_x, u_y])   # (N, 2) — no net CM displacement, no rotation


# ── Part 1: free-vibration with correct bending mode ─────────────────────────

def free_vibration_bending(tau, R0, N, alpha, S=1.0, rho_d=1.0,
                           A_init=0.01, n_tdecay=5.0):
    """
    Initialize correct n=2 bending mode in body frame, evolve under step_rb,
    measure ω_2_bend and T_decay.

    Returns dict with:
      omega2_theory, omega2_meas, T_decay_theory, T_decay_meas,
      omega2_err_pct, T_decay_err_pct, fit_ok
    """
    p  = make_particle(tau=tau, R0=R0, N=N, S=S, rho_d=rho_d)
    dt = get_dt(p)

    T_decay_theory = 2.0 / alpha
    T_2_bend       = 2.0 * np.pi / omega2_bend_theory(tau, R0, S, rho_d)
    # Record for enough bending cycles and enough decay e-folds
    T_record = max(n_tdecay * T_decay_theory, 6.0 * T_2_bend)
    n_steps  = max(int(T_record / dt), 1000)

    # Set correct bending mode initial condition (body frame, no rigid-body motion)
    p.u     = bending_mode_init(N, R0, A_init)
    p.u_dot = np.zeros((N, 2))
    p.v_cm  = np.zeros(2)
    p.omega = 0.0
    # Reconstruct world positions consistently
    p.x = p.x_cm + p.X_ref + p.u

    sim = CapsuleSim([p])

    # Track u_r projection onto cos(2θ): Q_2 = (2/N) Σ u_r_k cos(2θ_k)
    theta_ref = np.arctan2(p.X_ref[:, 1], p.X_ref[:, 0])
    cos2      = np.cos(2.0 * theta_ref)
    cos_t     = np.cos(theta_ref)
    sin_t     = np.sin(theta_ref)

    Q  = np.empty(n_steps)
    for k in range(n_steps):
        sim.step_rb(dt, alpha_damp=alpha)
        # u_r component projected onto cos(2θ)
        u_r  = p.u[:, 0] * cos_t + p.u[:, 1] * sin_t
        Q[k] = (2.0 / N) * float(np.sum(u_r * cos2))

    t = np.arange(n_steps) * dt

    # ── FFT for frequency ─────────────────────────────────────────────────────
    freqs  = np.fft.rfftfreq(n_steps, d=dt)
    power  = np.abs(np.fft.rfft(Q))**2
    i_peak = int(np.argmax(power[1:])) + 1
    omega2_fft = 2.0 * np.pi * freqs[i_peak]

    # ── Damped-sinusoid fit ───────────────────────────────────────────────────
    fit_ok = False
    T_decay_meas = np.nan
    omega2_meas  = omega2_fft

    try:
        def model(t, A, alpha_fit, omega_fit, phi):
            return A * np.exp(-0.5 * alpha_fit * t) * np.cos(omega_fit * t + phi)

        # Use FFT peak as initial guess; allow ±5× range for discrete corrections
        w0 = max(omega2_fft, 0.1)
        p0 = [float(abs(Q[0]) + 1e-6), alpha, w0, 0.0]
        hi_w = max(10.0 * omega2_bend_theory(tau, R0, S, rho_d), 5.0 * w0)
        bounds = ([0.0, 0.0, 0.1, -np.pi],
                  [20.0 * A_init, 20.0 * alpha + 1.0, hi_w, np.pi])
        popt, _ = curve_fit(model, t, Q, p0=p0, bounds=bounds,
                             maxfev=20_000, ftol=1e-8)
        A_fit, alpha_fit, omega2_fit, phi_fit = popt
        T_decay_meas = 2.0 / alpha_fit
        omega2_meas  = omega2_fit
        fit_ok = True
    except Exception as e:
        # Envelope method for T_decay
        try:
            q1 = float(np.sqrt(np.mean(Q[:n_steps//4]**2) + 1e-30))
            q4 = float(np.sqrt(np.mean(Q[3*n_steps//4:]**2) + 1e-30))
            t_mid = 0.5 * T_record
            if q1 > 1e-15:
                T_decay_meas = -t_mid / np.log(max(q4 / q1, 1e-10))
        except Exception:
            pass

    om2_th = omega2_bend_theory(tau, R0, S, rho_d)
    return dict(
        tau=tau, R0=R0, N=N, alpha=alpha, S=S, rho_d=rho_d,
        omega2_theory=float(om2_th),
        omega2_meas=float(omega2_meas),
        T_decay_theory=float(T_decay_theory),
        T_decay_meas=float(T_decay_meas),
        omega2_err_pct=float(100.0 * (omega2_meas - om2_th) / om2_th),
        T_decay_err_pct=(float(100.0 * (T_decay_meas - T_decay_theory) / T_decay_theory)
                         if np.isfinite(T_decay_meas) else float('nan')),
        Q_bend=float(om2_th / alpha),
        fit_ok=fit_ok,
        dt=float(dt),
        n_steps=n_steps,
    )


# ── Part 5: post-collision ringing ────────────────────────────────────────────

def post_collision_ringing(tau=0.20, R0=1.0, N=32, alpha=0.5, S=1.0,
                           rho_d=1.0, v0=2.0, C=3000.0):
    """
    Head-on collision; measure post-collision ringing from node-0 radial displacement.

    Uses a three-phase state machine to correctly:
      Phase 0: wait until contact (separation < contact_thr)
      Phase 1: wait until clean separation (separation > sep_clear)
      Phase 2: record n=2 bending oscillation

    Larger v0 = 2.0 ensures measurable elastic deformation amplitude.
    Low alpha = 0.5 → Q_bend = 12 → many visible cycles for clean FFT.
    """
    om2_th = omega2_bend_theory(tau, R0, S, rho_d)
    T_2    = 2.0 * np.pi / om2_th
    p_ref  = make_particle(tau=tau, R0=R0, N=N, S=S, C=C, rho_d=rho_d)
    dt     = get_dt(p_ref)
    n_post = max(int(12.0 * T_2 / dt), 2000)

    pA = make_particle(tau=tau, R0=R0, N=N, S=S, C=C, rho_d=rho_d)
    pB = make_particle(tau=tau, R0=R0, N=N, S=S, C=C, rho_d=rho_d)
    d_init = 2.0 * R0 + 4.0
    pA.translate([-d_init / 2, 0.0])
    pB.translate([ d_init / 2, 0.0])
    pA.v_cm = np.array([ v0 / 2.0, 0.0])
    pB.v_cm = np.array([-v0 / 2.0, 0.0])

    sim = CapsuleSim([pA, pB])
    contact_thr = pA.r_c + pB.r_c + pA.R0 + pB.R0 + 0.05
    sep_clear   = 2.0 * R0 + pA.r_c + pB.r_c + 0.40

    # Phase 0→1→2 state machine (same logic as measure_cor)
    phase = 0
    for _ in range(5_000_000):
        sim.step_rb(dt, alpha_damp=alpha)
        d = np.linalg.norm(pA.x_cm - pB.x_cm)
        if phase == 0 and d < contact_thr:
            phase = 1
        if phase == 1 and d > sep_clear:
            phase = 2
            break

    if phase != 2:
        return dict(tau=tau, R0=R0, N=N, alpha=alpha,
                    omega2_theory=float(om2_th), omega_ring=float('nan'),
                    omega_ring_err_pct=float('nan'), amp_init=0.0, dt=float(dt))

    # Record post-separation node-0 radial displacement (body frame)
    theta_ref = np.arctan2(pA.X_ref[:, 1], pA.X_ref[:, 0])
    cos0, sin0 = float(np.cos(theta_ref[0])), float(np.sin(theta_ref[0]))

    r0_trace = []
    for _ in range(n_post):
        sim.step_rb(dt, alpha_damp=alpha)
        u_r0 = float(pA.u[0, 0] * cos0 + pA.u[0, 1] * sin0)
        r0_trace.append(u_r0)

    arr = np.array(r0_trace)
    t   = np.arange(n_post) * dt

    freqs      = np.fft.rfftfreq(n_post, d=dt)
    power      = np.abs(np.fft.rfft(arr))**2
    i_peak     = int(np.argmax(power[1:])) + 1
    omega_ring = 2.0 * np.pi * freqs[i_peak]
    amp_init   = float(np.sqrt(np.mean(arr[:n_post//10]**2) + 1e-30))

    return dict(
        tau=tau, R0=R0, N=N, alpha=alpha,
        omega2_theory=float(om2_th),
        omega_ring=float(omega_ring),
        omega_ring_err_pct=float(100.0 * (omega_ring - om2_th) / om2_th),
        amp_init=float(amp_init),
        dt=float(dt),
    )


# ── COR measurement ───────────────────────────────────────────────────────────

def measure_cor(tau=0.20, R0=1.0, N=32, alpha=2.0, S=1.0, rho_d=1.0,
                v0=0.50, C=3000.0):
    """
    Head-on collision COR.  v_pre measured just before first contact;
    v_post measured once disks are cleanly separated.
    """
    pA = make_particle(tau=tau, R0=R0, N=N, S=S, C=C, rho_d=rho_d)
    pB = make_particle(tau=tau, R0=R0, N=N, S=S, C=C, rho_d=rho_d)
    d_init = 2.0 * R0 + 4.0
    pA.translate([-d_init / 2, 0.0])
    pB.translate([ d_init / 2, 0.0])
    pA.v_cm = np.array([ v0 / 2, 0.0])
    pB.v_cm = np.array([-v0 / 2, 0.0])

    contact_thr = pA.r_c + pB.r_c + pA.R0 + pB.R0 + 0.05
    sep_clear   = 2.0 * R0 + pA.r_c + pB.r_c + 0.40

    sim = CapsuleSim([pA, pB])
    dt  = get_dt(pA)

    # Phase tracking: 0=approach, 1=in_contact, 2=separated
    phase        = 0
    v_pre        = None
    T_contact_steps = 0

    for _ in range(4_000_000):
        sim.step_rb(dt, alpha_damp=alpha)
        d = np.linalg.norm(pA.x_cm - pB.x_cm)

        if phase == 0 and d < contact_thr:
            v_pre = (pA.v_cm.copy(), pB.v_cm.copy())
            phase = 1
        if phase == 1:
            if d < contact_thr:
                T_contact_steps += 1
            elif d > sep_clear:
                phase = 2
                break

    if v_pre is None or phase != 2:
        return dict(alpha=alpha, COR=float('nan'), T_contact=float('nan'),
                    separated=False, Q_bend=float('nan'))

    ex = np.array([1.0, 0.0])
    vr_pre  = float(np.dot(v_pre[0] - v_pre[1], ex))
    vr_post = float(np.dot(pA.v_cm - pB.v_cm, ex))
    COR     = abs(vr_post / (vr_pre + 1e-30))

    om2_bend = omega2_bend_theory(tau, R0, S, rho_d)
    return dict(
        alpha=float(alpha),
        COR=float(COR),
        T_contact=float(T_contact_steps * dt),
        T_decay=float(2.0 / alpha),
        Q_bend=float(om2_bend / alpha),
        separated=True,
        dt=float(dt),
    )


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("═══ Shell Mode Calibration (corrected) ═══\n")
    print("Correct observable mode: n=2 inextensional bending")
    print("    ω_2_bend = (6/√5) · √(S / (ρ_f · τ · R₀²))\n")

    tau_ref, R0_ref, N_ref, S_ref = 0.20, 1.0, 32, 1.0
    om2_ref = omega2_bend_theory(tau_ref, R0_ref, S_ref)
    print(f"Reference params (τ=0.20, R0=1, N=32, S=1, ρ=1):")
    print(f"  ω_2_bend = {om2_ref:.4f} rad/s   T_2 = {2*np.pi/om2_ref:.4f} s")
    print(f"  For Q<3 (few cycles): α > {om2_ref/3:.2f}")
    print(f"  For Q<1 (no ringing): α > {om2_ref:.2f}\n")

    all_results = {}

    # ─────────────────────────────────────────────────────────────────────────
    # Part 1 — Free vibration: correct n=2 bending initial condition
    # ─────────────────────────────────────────────────────────────────────────
    print("── Part 1: Free-vibration with correct bending mode IC (α=0.5) ──\n")
    ref = free_vibration_bending(tau=0.20, R0=1.0, N=32, alpha=0.5, S=1.0)
    print(f"  ω₂_bend theory = {ref['omega2_theory']:.4f}")
    print(f"  ω₂_bend meas   = {ref['omega2_meas']:.4f}  "
          f"  err = {ref['omega2_err_pct']:+.1f}%  "
          f"  fit_ok = {ref['fit_ok']}")
    print(f"  T_decay theory = {ref['T_decay_theory']:.5f}")
    print(f"  T_decay meas   = {ref['T_decay_meas']:.5f}  "
          f"  err = {ref['T_decay_err_pct']:+.2f}%")
    print(f"  dt = {ref['dt']:.2e},  n_steps = {ref['n_steps']}\n")
    all_results['reference_freevib'] = ref

    # ─────────────────────────────────────────────────────────────────────────
    # Part 2 — T_decay = 2/α universality
    # ─────────────────────────────────────────────────────────────────────────
    print("── Part 2: T_decay = 2/α universality (τ=0.20, R0=1, N=32) ──\n")
    print("  Note: T_decay = 2/α holds for underdamped modes (Q > 1).")
    print("  For overdamped (Q < 1): actual decay ~ α/ω² > 2/α (slower).\n")
    print(f"  {'α':>8s}  {'T_dec_th':>10s}  {'T_dec_ms':>10s}  {'err%':>8s}  "
          f"{'Q_bend':>8s}  regime")
    print("  " + "─" * 65)

    alpha_tdecay = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    tdecay_res = []
    for alpha in alpha_tdecay:
        r    = free_vibration_bending(tau=0.20, R0=1.0, N=32, alpha=alpha, S=1.0)
        err  = r['T_decay_err_pct']
        q    = r['Q_bend']
        if q > 1:
            regime = "underdamped (T=2/α valid)"
        else:
            # Overdamped: actual T_actual ≈ α/ω² analytically
            t_actual_theory = alpha / (om2_ref**2)
            regime = f"overdamped (T_act≈{t_actual_theory:.3f})"
        print(f"  {alpha:>8.1f}  {r['T_decay_theory']:>10.5f}  "
              f"{r['T_decay_meas']:>10.5f}  {err:>+8.2f}%  "
              f"{q:>8.2f}  {regime}")
        tdecay_res.append(r)
    all_results['tdecay_universality'] = tdecay_res
    print()

    # ─────────────────────────────────────────────────────────────────────────
    # Part 3 — Scaling: ω_2_bend ∝ √(S/(ρ_f·τ)) / R0
    # ─────────────────────────────────────────────────────────────────────────
    print("── Part 3: ω_2_bend scaling  (α=0.5 throughout) ──\n")

    scaling_configs = [
        # (label, tau, R0, S)
        # τ sweep (R0=1, S=1): ω ∝ 1/√τ
        ("τ=0.05, R0=1, S=1",  0.05, 1.0, 1.0),
        ("τ=0.10, R0=1, S=1",  0.10, 1.0, 1.0),
        ("τ=0.20, R0=1, S=1",  0.20, 1.0, 1.0),  # reference
        # R0 sweep (τ=0.20, S=1): ω ∝ 1/R0
        ("τ=0.20, R0=0.5, S=1", 0.20, 0.5, 1.0),
        ("τ=0.20, R0=2.0, S=1", 0.20, 2.0, 1.0),
        # S sweep (τ=0.20, R0=1): ω ∝ √S
        ("τ=0.20, R0=1, S=0.1", 0.20, 1.0, 0.1),
        ("τ=0.20, R0=1, S=4.0", 0.20, 1.0, 4.0),
    ]

    print(f"  {'config':30s}  {'ω_th':>7s}  {'ω_ms':>7s}  {'err%':>6s}  "
          f"{'α_Q1':>8s}  {'α_Q3':>8s}")
    print("  " + "─" * 76)

    scaling_results = []
    for label, tau, R0, S in scaling_configs:
        t0 = time.time()
        r  = free_vibration_bending(tau=tau, R0=R0, N=32, alpha=0.5, S=S)
        alpha_q1 = r['omega2_theory']
        alpha_q3 = r['omega2_theory'] / 3.0
        print(f"  {label:30s}  {r['omega2_theory']:>7.3f}  {r['omega2_meas']:>7.3f}  "
              f"{r['omega2_err_pct']:>+6.1f}%  "
              f"{alpha_q1:>8.3f}  {alpha_q3:>8.3f}  [{time.time()-t0:.0f}s]")
        r['label'] = label
        r['alpha_for_Q1'] = float(alpha_q1)
        r['alpha_for_Q3'] = float(alpha_q3)
        scaling_results.append(r)
    all_results['scaling'] = scaling_results
    print()

    # ─────────────────────────────────────────────────────────────────────────
    # Part 4 — COR vs α
    # ─────────────────────────────────────────────────────────────────────────
    print("── Part 4: COR vs α (τ=0.20, R0=1, N=32, S=1, v0=0.5, C=3000) ──\n")
    print(f"  ω_2_bend = {om2_ref:.3f}   α_Q1 = {om2_ref:.1f}   α_Q3 = {om2_ref/3:.1f}\n")
    print(f"  {'α':>7s}  {'Q_bend':>7s}  {'T_dec':>7s}  {'T_cont':>7s}  "
          f"{'COR':>6s}  ringing?")
    print("  " + "─" * 62)

    # Span from very underdamped to well overdamped
    alpha_cor_vals = [0.3, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 50.0]
    cor_results = []
    for alpha in alpha_cor_vals:
        t0  = time.time()
        r   = measure_cor(alpha=alpha)
        q   = om2_ref / alpha
        if q > 5:
            note = "persistent ringing"
        elif q > 3:
            note = "moderate ringing"
        elif q > 1:
            note = "few cycles"
        else:
            note = "no visible ringing"
        print(f"  {alpha:>7.1f}  {q:>7.2f}  {2/alpha:>7.4f}  "
              f"{r['T_contact']:>7.4f}  {r['COR']:>6.4f}  ← {note}"
              f"  [{time.time()-t0:.0f}s]",
              flush=True)
        cor_results.append(r)
    all_results['cor_vs_alpha'] = cor_results

    # ─────────────────────────────────────────────────────────────────────────
    # Part 5 — Post-collision ringing (direct)
    # Use v0=2.0 to ensure measurable elastic deformation during collision.
    # ─────────────────────────────────────────────────────────────────────────
    print("\n── Part 5: Post-collision ringing frequency (α=0.5, v0=2.0) ──\n")
    t0  = time.time()
    pcr = post_collision_ringing(tau=0.20, R0=1.0, N=32, alpha=0.5, S=1.0,
                                  v0=2.0, C=3000.0)
    print(f"  ω theory = {pcr['omega2_theory']:.4f} rad/s")
    print(f"  ω ring   = {pcr['omega_ring']:.4f} rad/s  "
          f"  err = {pcr['omega_ring_err_pct']:+.1f}%")
    print(f"  Initial ringing amplitude = {pcr['amp_init']:.3e}")
    print(f"  [{time.time()-t0:.0f}s]")
    all_results['post_collision_ringing'] = pcr

    # ─────────────────────────────────────────────────────────────────────────
    # Summary
    # ─────────────────────────────────────────────────────────────────────────
    print("\n\n═══ DESIGN GUIDELINES (corrected) ═══\n")
    print("  Observable post-collision mode: n=2 inextensional bending")
    print(f"  ω_2_bend = (6/√5) · √(S / (ρ_f · τ · R0²))\n")
    print("  α choice rule (general):")
    print("    α = ω_2_bend / Q_target")
    print("    = (6/√5) · √(S / (ρ_f · τ)) / (Q_target · R0)\n")
    print("  Common targets:")
    print(f"    Q=1 (no ringing):    α = {om2_ref:.2f}  at reference")
    print(f"    Q=3 (few cycles):    α = {om2_ref/3:.2f}  ← current default")
    print(f"    Q=10 (many cycles):  α = {om2_ref/10:.2f}")
    print()
    print("  Dimensionless form:")
    print("    α̃ = α · R0 · √(ρ_f · τ / S) · √5/6")
    print("    α̃ = 1  →  Q = 1 (critical damping for bending mode)")
    print(f"    At reference: α̃ = α / {om2_ref:.2f}\n")

    # Report COR sensitivity
    cor_low  = next((r for r in cor_results if r['alpha'] <= 0.5), None)
    cor_high = next((r for r in reversed(cor_results) if r['alpha'] >= 20.0), None)
    if cor_low and cor_high:
        dcor = abs(cor_high['COR'] - cor_low['COR'])
        print(f"  COR sensitivity: ΔCOR = {dcor:.4f} from α={cor_low['alpha']} → α={cor_high['alpha']}")
        if dcor < 0.05:
            print("  COR is insensitive to α (ΔCOR < 5%) — T_contact << T_2_bend")
        else:
            print("  COR varies significantly with α — T_contact ≈ T_2_bend")

    # ─────────────────────────────────────────────────────────────────────────
    # Save
    # ─────────────────────────────────────────────────────────────────────────
    os.makedirs('results/shell_mode_calibration', exist_ok=True)
    with open('results/shell_mode_calibration/calibration_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print("\nSaved → results/shell_mode_calibration/calibration_results.json")


if __name__ == '__main__':
    main()
