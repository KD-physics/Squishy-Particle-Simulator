---
name: Phase 1G — 3-body collision and squeeze test benchmarks
description: Phase 1G results: machine-precision momentum in 3-body collision, T1-event squeeze test with oscillating wall, Arc geometry correction
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Phase 1G: Extension Benchmarks (2026-04-20)

## Status: COMPLETE ✅

Script: `src/validation/rb_extension_benchmarks.py`
Movies: `src/visualization/rb_extension_movies.py`
Results: `results/rb_benchmarks/extension_benchmarks.json`, `extension_1G1.png`, `extension_1G2.png`

## 1G.1 — 3-body Asymmetric Collision

**Setup:** A at (−b, y_A) and B at (+b, y_B) fall, C at (x_off, 0) rises. b=1.5, x_off=0.5.
ε offsets computed so all three contact simultaneously:
```
y_A = V * t_star + sqrt(contact_thr² − (b+x_off)²)
```
where t_star = (y_B0 - sqrt(contact_thr² − (b−x_off)²)) / V

**Results:** All 5 configs (q ∈ {0.1, 1, 10}, α ∈ {0, 2}, asymmetric v):
- |Δp/p0| ≤ 1.6e-15 (machine precision)
- |ΔL/L0| ≤ 1.3e-14 (machine precision)
- COR_eff: 0.800 (damped) → 0.924 (glass elastic)

**Key:** Central-force contact → exact angular momentum conservation by construction.

## 1G.2 — Rearrangement Squeeze Test (T1-event analogue)

**Setup:** Center EPD particle (N=16) pushed by oscillating wall v(t)=v₀+A·sin(ωt) through
fixed arc obstacles representing two side EPD particles.

**Arc geometry fix (CRITICAL):**
- Use `arc_radius = R0 + r_c` (NOT R0) to match particle-particle contact distance
- `contact_thr_pp = 2*R0 + 2*r_c` (same as EPD-EPD contact threshold)
- `d_side = contact_thr_pp − squeeze_gap` (guaranteed contact at y=0 when squeeze_gap > 0)
- If `arc_radius = R0` is used: all rightmost nodes are deeply inside the arc at d_side > R0

**N=16 for squeeze:** dt is ~3× larger than N=32 (edge-wave limited), making simulations ~3× faster.

**Results (4 configs: q ∈ {0.1, 1, 10}, squeeze_gap ∈ {0.05, 0.10}):**
- All particles pass through gap ✓
- E_elastic_max < W_wall/4 (no blow-up) ✓
- PE at exit << peak PE (clean exit) ✓
- Rubber (q=10): ~6× more wall work than glass (q=0.1) at same squeeze_gap
  (larger r_c → larger contact bead → more force needed to squeeze through)

**Why:** The rubber particle's r_c = L0 = 2πR0/N is larger because r_c = L0 (full edge length, NOT L0/2).

## Contact Shell Geometry Note

r_c = L0 (full edge length, not half). For N=32: r_c ≈ 0.196. For N=16: r_c ≈ 0.393.
contact_thr = 2*R0 + 2*r_c (not 2*R0 + r_c as might be expected for one-sided contact).

## Energy Balance

W_wall_energy = (E_final − E_initial) + W_diss
This is the definition — wall does work that goes into KE + PE + dissipation.
No independent force tracking needed to verify this holds.
