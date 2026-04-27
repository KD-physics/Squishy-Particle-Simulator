---
name: alpha_damp_calibration
description: Correct α_damp design rules based on n=2 inextensional bending mode (NOT membrane mode). α-COR tradeoff quantified.
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
## Shell damping α_damp — corrected calibration (Phase 1J, 2026-04-21)

**Why:** Prior entry (Phase 1I) used membrane formula (ω_ext ≈ 78 rad/s) — wrong.
Observable post-collision ringing is the n=2 inextensional bending mode at ω_bend ≈ 6 rad/s.
Verified with 5-part numerical test suite (all parts PASS).

**Correct bending mode formula:**

    ω_2_bend = (6/√5) · √(S / (ρ_f · τ · R0²))

Where EI = S·K_fluid·R0³ (from capsule_shell.py line 136), ρ_L = ρ_f·τ·R0.

At reference (S=1, τ=0.20, R0=1, ρ=1): **ω_2_bend = 6.000 rad/s** (measured 6.0002, +0.0%).

**Scaling verified** (all within 2%):
- ω ∝ 1/√τ  (matches theory)
- ω ∝ 1/R0   (matches theory)
- ω ∝ √S     (matches theory)
- N-independent in continuum limit

**T_decay:**
- T_decay = 2/α for Q > 1 (underdamped); error < 0.1%
- For Q < 1 (overdamped): T_actual ≈ α/ω_bend² >> 2/α (decays SLOWER than 2/α)

**Critical finding — α-COR tradeoff:**
T_contact ≈ 1.57 ≈ 1.5 × T_2_bend → bending mode resonates DURING contact.
COR is STRONGLY sensitive to α (this is physical, not numerical):

| α   | Q_bend | COR  |
|-----|--------|------|
| 0.3 | 20     | 0.894|
| 1.0 |  6     | 0.825|
| 2.0 |  3     | 0.729|  ← current default
| 6.0 |  1     | 0.530|
| 20  | 0.3    | 0.289|

**Design formula:**

    α = ω_2_bend / Q_target = (6/√5) · √(S/(ρ_f·τ)) / (Q_target · R0)

Common choices at reference:
- Q=3 (few ringing cycles, current default): α = 2.0, COR ≈ 0.73
- Q=1 (no visible ringing): α = 6.0, COR ≈ 0.53
- For dense packings needing high COR: accept Q~6-12, use α=0.5-1, but expect ringing

**Dimensionless form:**
    α̃ = α · R0 · √(ρ_f·τ/S) · √5/6
    α̃ = 1 → Q = 1 (critical damping)

**Script:** `src/validation/shell_mode_calibration.py`
**Results:** `results/shell_mode_calibration/calibration_results.json`
