---
name: Elastic Perimeter Disk Model — Squishiness Axis and Calibration
description: EPD model canonical parameterization: b=0.2 working point, q→ν calibration curve, R0 invariance rules. Phase-space sweep and fixed-b sweep results.
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Elastic Perimeter Disk (EPD) Model: Squishiness Calibration (2026-04-19)

## Model Name
**Elastic Perimeter Disk (EPD) model** — chosen to be descriptive without implying "toy model".

## Two Dimensionless Physical Parameters

- **b = EI/El_t_BASE = τ²/12**  (bending number, at R0=1)
- **q = K_area/El_t_BASE**  (area stiffness ratio)
- **q·b = 1** is the natural midpoint dividing area-squishy (q·b < 1) from shape-squishy (q·b > 1)

## Canonical Working Point: b = 0.2

At b=0.2 (τ = √2.4 ≈ 1.5492), sweeping q gives a clean 1D squishiness axis:

| q    | q·b  | ΔA@8%  | ν@8%  |
|------|------|--------|-------|
| 0.05 | 0.01 | 6.2%   | 0.18  |
| 1.0  | 0.20 | 3.2%   | 0.56  |
| 5.0  | 1.00 | 1.0%   | 0.82  |
| 50.0 | 10.0 | 0.12%  | 0.94  |

ΔA and ν are one-to-one at fixed b → q uniquely determines both.

## Shape Collapse Result
At fixed ΔA, perimeter shape is geometry-determined (independent of τ, q).
Reason: El_t cancels as a prefactor in the E-L equations for perimeter shape at fixed area constraint.
b=0.2 is membrane-dominated enough that shape collapse holds to visual accuracy.

## R0 Invariance Rules
For physics to be R0-independent:
- **EI = S × R0³** (scales as R0³, from thin-shell formula)
- **K_area = q × El_t_BASE** (CONSTANT, does NOT scale with R0)
- **C = 3000 × S × (1+q)** (CONSTANT, does NOT scale with R0)

Verified: CV < 0.01% for ΔA and ν across R0 ∈ {0.25, 0.5, 1.0, 2.0, 4.0}.

## Parameter Recipe: Given (ν_target, R0)
1. q = q(ν_target) from calibration curve (log-linear interpolation of ν vs log(q))
2. τ = √2.4 ≈ 1.5492 (fixed for b=0.2)
3. El_t = 12*S/τ²  (base, R0=1 value)
4. EI = S × R0³
5. K_area = q × El_t  (constant, no R0 factor)
6. C = 3000 × S × (1+q)  (constant, no R0 factor)

## N Sensitivity (COMPLETE — 2026-04-19)
N ∈ {32, 48, 72}, Q: 18 log-spaced values, ε_ref ∈ {0.04,...,0.12}
Script: src/validation/calibration_sweep.py
Output: results/calibration_sweep/calibration_data.json
Result: max Δν=0.0095 across N={32,48,72}; median=0.0064. N=48 calibration transfers to N=72.

## ε_ref Sensitivity (COMPLETE — 2026-04-19)
At N=48: max Δν=0.045 across ε_ref∈{0.04,...,0.12}, median=0.033.
Monotone increase; ε_ref=0.08 chosen as working point (well-developed contact, quasi-linear regime).

## R0 Spot Check (COMPLETE — 2026-04-19)
Script: src/validation/R0_spot_check.py
Tests: ν_target ∈ {0.30, 0.60, 0.85} at R0 ∈ {0.5, 1.0, 2.0}
Result: CV(ν)=0.00%, CV(ΔA)≤0.01% — exact R0 invariance confirmed.

## Canonical Calibration Table (N=48, ε_ref=0.08)
ν range: 0.173 (q=0.05) → 0.939 (q=50)
Key anchor: q=5.0 → ν=0.824, ΔA=1.05% (midpoint q·b=1)
Inversion: log-linear interpolation ν_target → log(q) via scipy.interpolate.interp1d

## Documentation
Model documented in: docs/model_and_numerics.md
Sections: §0 Motivation → §1 Geometry → §2 Energy → §3 EOM → §4 Dimensional analysis →
          §5 Calibration → §6 Phase space + shape collapse → §7 Squishiness axis → §8 Recipe
