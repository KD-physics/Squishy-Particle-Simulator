---
name: Capsule Shell DEM — Dimensional Scaling Laws
description: Proven scaling laws for the Capsule Shell DEM: S=force scale (CV=0%), R0=length scale (CV=0% after EI bug fix). Two free parameters are tau and q.
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Capsule Shell DEM: Dimensional Scaling Laws (2026-04-18)

## Result 1: S is a pure force scale (Phase 1C)

**Script:** `src/validation/scaling_verification.py`

At fixed q=1, τ=0.05, C/El_t=0.208, α₀=2.0, R0=1:
Sweeping S ∈ {0.01, 0.1, 1.0, 10.0} gives **CV=0.00%** for F/El_t, ΔA, ν_meas at ε_p=0.05,0.08,0.10.

**Implication:** Fix S=1 with no loss of generality. El_t_base = 12/τ².

## Result 2: R0 is a pure length scale (Phase 1C)

**Script:** `src/validation/R0_N_convergence.py`

At fixed S=1, q=1, τ=0.05, C=1000 (const), K_area=4800 (const), α₀=2.0:
Sweeping R0 ∈ {0.25, 0.5, 1.0, 2.0, 4.0} at N=32, 64, 128:

| N   | F_norm (all R0) | ΔA@10% (all R0) | CV% |
|-----|-----------------|-----------------|-----|
| 32  | 0.00398         | 0.365%          | 0.00 |
| 64  | 0.00432         | 0.378%          | 0.00 |
| 128 | 0.00504         | 0.399%          | 0.00 |

**Implication:** Fix R0=1 with no loss of generality.

## Bug fixed: EI = S × R0³ (was S × R0⁴)

**File:** `src/simulation/capsule_shell.py` line 129

For a thin shell with t = τR0: EI = E_l × t³/12 = S × K_fluid × **R0³**.
Old code had R0⁴, breaking R0 self-similarity (bending grew as R0² relative to membrane).
At R0=1: EI unchanged → all Phase 1/1B calibrations intact.

**How to apply:** Always use the corrected formula. Never revert to R0⁴.

## Correct scaling rules for R0-invariant simulations

- **EI = S × K_fluid × R0³** (fixed)
- **C = const** — contact force: k_c×gap×L0 = (C/R0)×R0×R0 = C×R0 ∝ R0 ✓
- **K_area = Q × El_t_base** (const) — area force ∝ K_area×R0 ∝ R0 ✓
- **Pass alpha0 to run_squeeze_once** — code computes α=α₀/T_wave internally (ζ=const)

## Canonical parameter set (post Phase 1C + C calibration)

| Parameter | Value | Notes |
|-----------|-------|-------|
| S    | 1.0 (fixed) | Force scale |
| R0   | 1.0 (fixed) | Length scale |
| C    | 3000×S×(1+q) | τ-INDEPENDENT; cc_pen<2% all (τ,q) |
| α₀   | 2.0 (fixed) | ζ = α₀/(4π) ≈ 16% |
| **τ** | free | Contact flatness knob |
| **q** | free | Squishiness knob (K_area/El_t_base) |

## q definition

q = K_area / El_t_BASE = K_area / (12S/τ²)

Use El_t_BASE (= 12S/τ² at R0=1), NOT El_t_code (= El_t_base × R0).
With K_area = const, q is R0-independent. ✓

## Why C = const w.r.t. R0 (not C × R0)

k_c = C/R0. Contact force = k_c × gap × L0 ~ (C/R0) × R0 × R0 = C × R0.
Force_unit = El_t_base × R0. So F/force_unit = C/El_t_base = const only if C = const.
If C ~ R0 (wrong): k_c = const, F ~ k_c × R0² → F/force_unit ~ R0 (breaks scaling).

## Why C must be τ-independent (not C = C_RATIO × El_t_base)

EI = S × K_fluid × R0³ is τ-INDEPENDENT at fixed S, R0. The contact spring k_c = C/R0 must
dominate EI/R0² = S. Hence C >> S × R0 → C ∝ S, independent of τ.
Old formula C = C_RATIO × El_t_base ∝ 1/τ² gave C=62 at τ=0.20 (vs C=1000 at τ=0.05),
causing cc_pen=22% at τ=0.20. Fix: C = 3000 × S × (1+q), constant w.r.t. τ.
The (1+q) term accounts for the area spring stiffening.

## C calibration constant

C_0 = 3000 (at S=1, R0=1) was determined by bisection: C_0=3000 gives cc_pen < 2%
for all tested (τ ∈ {0.05,0.10,0.20}, q ∈ {0.1,1,5,10}).
Script: src/validation/C_formula_verify.py
