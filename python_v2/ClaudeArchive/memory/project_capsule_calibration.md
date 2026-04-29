---
name: Capsule Shell DEM — Parameter Calibration Findings
description: Empirical calibration results for the Capsule Shell DEM: C, alpha0, squishiness axis q=K_area/El_t, tau vs contact flatness
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Capsule Shell DEM: Calibration Findings (2026-04-18)

## Phase 1B: Squishiness Axis Calibration (COMPLETE)

**Squishiness knob:** q = K_area / El_t  where El_t = 12*S/tau²

**Why:** ΔA (fractional area loss) is the primary squishiness metric. High q → area-preserving
(rubber); low q → area-changing (glass). ν_meas counter-intuitively DECREASES with q because
high K_area keeps disk circular → equatorial nodes don't expand.

**Calibration anchors (tau=0.05, S=0.1):**
- q_glass  = 0.1  → ΔA@8% ≈ 0.94%  (area-changing, glass-like)
- q_rubber = 10.0 → ΔA@8% ≈ 0.10%  (area-preserving, rubber-like)

**Contact flatness vs tau (S=0.1, q=10 rubber-like):**
- tau=0.05 → theta_wall = 9.67°
- tau=0.10 → theta_wall = 7.45°
- tau=0.20 → theta_wall = 5.15°   ← flattest, recommended for flat-contact studies

**Area force sign fix:** `self.f -= Fp` (CORRECT: outward restoring force when A < A₀).
Old code `self.f += Fp` was wrong (inward = destabilizing).

**How to apply:** Use q_glass=0.1, q_rubber=10 as anchors. Use tau=0.20 for flat-contact
studies, tau=0.05 for pointy-contact studies. K_area = q * 12 * S / tau².

---

## Phase 1: C Calibration (N=32, tau=0.05, alpha0=2.0)

**Core result: C=100 works for all S values (0.005–10.0)**
Counter-intuitively, higher C (1000+) causes MORE apparent penetration due to underdamped
contact oscillations. C=100 gives pen_frac < 1.1% across all S.

**Note for S=1.0:** At intermediate q (0.5–5), pen_frac rises to 11–16% even with C=100.
May need higher C for S=1.0 production runs at these q values.

| S     | C_cal | worst_cc_gap/r_c |
|-------|-------|-----------------|
| 0.005 | 100   | 0.001           |
| 0.01  | 100   | 0.002           |
| 0.05  | 100   | 0.007           |
| 0.1   | 100   | 0.011           |
| 0.5   | 100   | 0.010           |
| 1.0   | 100   | 0.000           |
| 3.0   | 100   | 0.000           |
| 10.0  | 100   | 0.000           |

---

## Performance

- dt_factor=0.4 (direct run_squeeze_once): ~27-50s per run at N=32
- Adaptive (run_squeeze with dt_factor=0.05): ~8x slower
- rho_f must be 1.0 (not 0.0) — estimate_dt_max divides by rho_f
- Steps per run: ~15,000 (constant regardless of S)

## Alpha0 Auto-Scaling

alpha_damp = alpha0 / T_wave keeps damping ratio ζ = alpha0/(4π) constant across all S.
- alpha0=2.0 → ζ ≈ 16% (good default)

## Results Files
- results/stage1b/step_A_results.json — full q sweep (7q × 2S × 3eps)
- results/stage1b/step_B_results.json — tau sweep (3tau × 2q × 2S)
- results/stage1b/step_C_four_metrics_tau0.050.png — four-metric summary
- results/stage1b/movies/ — 4 GIFs (glass/rubber × tau=0.05/0.20)
