---
name: EllipseParticle — design and benchmark results
description: EllipseParticle class design, equal arc-length parameterisation, and 4-config collision benchmark showing machine-precision angular momentum conservation
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# EllipseParticle (src/simulation/ellipse_particle.py)

## Design: Subclass of CapsuleParticle

EllipseParticle subclasses CapsuleParticle with:
- `R_eff = sqrt(a*b)` passed to parent (sets elastic parameters)
- `X_ref` patched to equal-arc-length ellipse nodes
- `L0 = r_c = perimeter/N` (arc-length per segment)
- `M_disk = rho_f * pi * a * b`
- `I_disk = M_disk * (a^2 + b^2) / 4`
- `R0 = max(a, b)` (broadphase bounding radius)
- First node pinned at semi-major axis tip: (a, 0) in body frame

All of step_rb(), contact forces, and movies work unchanged.

## Equal arc-length parameterisation (ellipse_arclength_nodes)

1. Tabulate ds/dt = sqrt((a sin t)^2 + (b cos t)^2) on fine grid (n_quad=8192)
2. Cumulative arc length via trapezoidal rule → s_cum
3. Invert via linear interpolation: t_of_s = interp1d(s_cum, t_extended)
4. Node positions at s = k * L0 for k = 0, 1, ..., N-1

## Benchmark results (ellipse_collision_benchmark.py)

4 configs all PASS gate |Δp|/p₀ < 1e-8, |ΔL|/L₀ < 1e-8:

| Config | Description | |Δp|/p₀ | |ΔL|/L₀ |
|--------|-------------|---------|---------|
| E1 | Equal, head-on, no spin | 5.7e-17 | 3.6e-16 |
| E2 | Equal, oblique+spinning | 5.6e-17 | 8.1e-13 |
| E3 | Diff size, same ρ=1 | 6.2e-14 | 2.8e-13 |
| E4 | Same I, diff mass (ρ_B≈1.184) | 1.1e-14 | 6.0e-13 |

**Why:** Contact force F ∝ (xq_A − cp_B). Torque = (xq_A − cp_B) × F ∝ v × v = 0 exactly.
Newton 3rd law pair → zero net torque on system for any particle shape.

## E4 density formula

For same I with different ellipses:
```python
rho_B = I_A / (pi * a_B * b_B * (a_B^2 + b_B^2) / 4)
```
For a_A=1.5, b_A=0.8, rho_A=1.0 vs a_B=1.2, b_B=1.0: rho_B ≈ 1.184

**Why:** Demonstrates angular momentum conservation is not an artifact of equal mass/inertia.
**How to apply:** Use EllipseParticle for any shape-dependent rheology studies. Movies in results/ellipse_benchmarks/movies/.
