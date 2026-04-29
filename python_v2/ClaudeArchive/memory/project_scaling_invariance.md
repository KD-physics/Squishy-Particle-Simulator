---
name: scaling_invariance
description: Scaling law for linear elasticity contact solver and its use as a data augmentation operator
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
The 2D linear-elasticity contact solver satisfies exact scale invariance:

**Rule**: Under X → λX with E, ν fixed (same δ/R):
- Displacements: u → λu
- Forces [N/m]: F → λF
- Strains ε, stresses σ, contact pressure P: unchanged
- Contact half-width: a → λa

**Verified**: Exact floating-point equality (rel_err = 0.00e+00) for all quantities,
all λ tested, all geometries (two-disk, wall contact, single-disk).

**Mechanism**: make_disk_mesh uses refine_area ∝ R² → identical topology, exact coordinate
scaling. K = ∫B^T C B dA is scale-invariant (B∝1/R, dA∝R² cancel). Forces f = k_pen·pen
scale as λ. Therefore u = K⁻¹f scales as λ.

**Data augmentation** (exact, zero compute cost):
Given sample (X_peri, u_peri, f_peri, E, ν), augmented at scale λ:
  (λ·X_peri, λ·u_peri, λ·f_peri, E, ν)
This is an exact reproduction of running the simulation at scale λ.
Cost: O(N_peri) multiply vs O(N_nodes²) FEM solve.
Suggested range: λ ∈ [0.5, 2.0] log-uniform for polydisperse particle augmentation.

**Why:** K is scale-invariant; f→λf; ∴ u=K⁻¹f→λu (proven in `scaling_validation.py`).
**How to apply:** Implement in `dataset.py` loader: sample λ, multiply X_peri/u/f by λ.
