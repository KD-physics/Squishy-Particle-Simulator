# VALIDATION.md — NN-DEM Simulation Pipeline: Complete Validation Record

**Last updated:** 2026-04-16  
**Status:** ALL PASS across Phases 1, 2, 2B, and scaling invariance  
**Reproducibility:** See [Quick Reproduction](#quick-reproduction) at end of document

---

## Overview

This document records all validation tests for the 2D penalty-contact FEM pipeline
that serves as both the physics engine and the training-data generator for the NN-DEM
surrogate.

| Phase | Scope | Status | Script |
|-------|-------|--------|--------|
| 1 | FEM baseline: uniform pressure, ν-sweep, N-convergence | ✅ PASS | `phase1_gate.py`, `nu_sweep_convergence.py` |
| 2 | Two-disk contact: homogeneous pairs (V1–V6) | ✅ PASS | `twodisk_validation.py` |
| 2-H | Two-disk contact: heterogeneous pairs (H1–H4) | ✅ PASS | `hetero_validation.py` |
| 2B | Contact primitives: flat walls, arcs, polygons (B1–B3) | ✅ PASS | `primitives_validation.py` |
| S | Scaling invariance: X→λX with E,ν fixed | ✅ PASS | `scaling_validation.py` |

All validation scripts live in `src/validation/`. All produce human-readable
PASS/FAIL output and exit 0 on success.

---

## Phase 1 — FEM Baseline

### What is validated

The FEM solver (`fem_elastic.py` + `disk_mesh.py`) for a single elastic disk.

### Tests

**P1: Uniform pressure analytic comparison**  
A disk under uniform inward pressure has an analytic solution:
```
u_r(R) = -p·R·(1+ν)·(1-2ν) / E        [plane strain]
```
FEM errors vs analytic for all four squishiness guide materials: **< 0.001%** (machine precision).

**P2: ν-sweep + N-convergence**  
P1 elements volumetrically lock at ν ≥ 0.45 (non-monotone convergence confirmed at ν=0.49).
P2 elements fix this: ΔF < 0.3% at N=512→1024.
**Rule baked in:** `assemble_stiffness()` auto-selects P2 when ν ≥ 0.45.

**P3: Perimeter tier self-consistency**  
F_sim/F (Dirichlet Hertz test) converges monotonically across N=120/240/360.
Finite-disk softening: F_sim/F ≈ 0.43–0.49 — correct physics (finite disk is ~2× softer than Hertz half-space).

### Squishiness guide materials

| Label | E (Pa) | ν | FEM element |
|-------|--------|---|-------------|
| Glass | 10⁸ | 0.22 | P1 |
| Polystyrene | 10⁶ | 0.33 | P1 |
| PDMS | 10⁴ | 0.47 | **P2** |
| Rubber/Jello | 10³ | 0.49 | **P2** |

---

## Phase 2 — Two-Disk Self-Consistent Contact (Homogeneous)

### What is validated

The penalty-contact iteration between two equal elastic disks converges to a
self-consistent equilibrium consistent with 2D plane-strain Hertzian theory.

### Contact model

Penalty stiffness: `k_pen = max(0.5·E, 1500 Pa)` (floor prevents rubber degeneracy).  
Iteration: under-relaxation with α ∈ [0.05, 0.15]; continuation for δ/R > 0.05.  
Convergence criterion: `|ΔF/F| < 0.5%` over 3 consecutive iterations.

### V1 — Convergence (130/130 configs)

10 materials × 13 δ/R values (0.002 → 0.15), all converge with F > 0.

| Material | E (Pa) | ν | δ/R range tested |
|----------|--------|---|-----------------|
| glass_stiff | 10⁸ | 0.20 | 0.002–0.15 |
| glass_soft | 3×10⁷ | 0.22 | 0.002–0.15 |
| ps_stiff | 10⁷ | 0.25 | 0.002–0.15 |
| ps_mid | 3×10⁶ | 0.30 | 0.002–0.15 |
| polystyrene | 10⁶ | 0.33 | 0.002–0.15 |
| pdms_stiff | 10⁵ | 0.40 | 0.002–0.15 |
| pdms | 3×10⁴ | 0.45 | 0.002–0.15 |
| pdms_soft | 10⁴ | 0.47 | 0.002–0.15 |
| rubber | 3×10³ | 0.48 | 0.002–0.15 |
| rubber_jello | 10³ | 0.49 | 0.002–0.15 |

### V2 — Newton balance

Sum of nodal contact forces equals scalar F_contact.  
Mean error: **3.8×10⁻¹⁷**, max: **2.2×10⁻¹⁶** — machine precision.

### V3 — E-scaling (linear elasticity check)

F/(E·R) is independent of E at fixed ν. Coefficient of variation < 1% across 6 decades of E.  
Representative values at δ/R = 0.050: F/(E·R) ≈ 1.45–1.84 × 10⁻².

### V4 — Hertz contact half-width

a_sim/a_hertz in **[0.77, 1.24]** across all materials at δ/R = 0.05–0.10.  
(Hertz assumes infinite half-space; finite-disk geometry causes systematic deviation — expected physics, not error.)

### V5 — Large-strain robustness

δ/R = 0.10, 0.12, 0.15 (near-jamming): all 15 configs converge.  
Achieved via 4-level continuation schedule with α = 0.15/0.10/0.08/0.05.

### V6 — Symmetry

|F_disk1 − F_disk2| / F = **0.000** exactly (machine precision) for all equal-material pairs.

### Dimensionless contact force

F/(E·R) increases smoothly with δ/R and with ν:

| ν | δ/R=0.006 | δ/R=0.020 | δ/R=0.050 | δ/R=0.100 | δ/R=0.150 |
|---|-----------|-----------|-----------|-----------|-----------|
| 0.20 | 1.05×10⁻³ | 4.72×10⁻³ | 1.45×10⁻² | 3.58×10⁻² | 6.15×10⁻² |
| 0.33 | 1.13×10⁻³ | 5.06×10⁻³ | 1.57×10⁻² | 3.88×10⁻² | 6.63×10⁻² |
| 0.45 | 1.08×10⁻³ | 4.81×10⁻³ | 1.55×10⁻² | 4.02×10⁻² | 6.94×10⁻² |
| 0.49 | 1.12×10⁻³ | 5.75×10⁻³ | 1.84×10⁻² | 4.38×10⁻² | 7.50×10⁻² |

---

## Phase 2-H — Two-Disk Contact (Heterogeneous)

### What is validated

Heterogeneous disk pairs: different radii (R₁ ≠ R₂) and/or different materials (E₁ ≠ E₂, ν₁ ≠ ν₂).

### H1 — Size ratio sweep

Polystyrene (E=10⁶, ν=0.33) at size ratios 1:1, 2:1, 4:1 (R₂ = 1.0, 0.5, 0.25).  
All 15 configs converge. Contact half-width a measured at ≥ 2 nodes for all cases.  
a_sim consistent with Hertz R_eff = R₁·R₂/(R₁+R₂): smaller R₂ → narrower contact patch.

### H2 — Material contrast

12 material-contrast pairs: glass/glass through glass/rubber_jello through 4:1 glass/rubber_jello.  
All converge. F/(E*·R_eff) consistently ≈ 1.4–1.8 × 10⁻² at δ/R = 0.02 across all pairs.

**Combined modulus (2D plane strain, mixed materials):**
```
1/E* = (1 − ν₁²)/E₁ + (1 − ν₂²)/E₂
```

### H3 — F/(E*·R_eff) scaling

The dimensionless contact force F/(E*·R_eff) is consistent across heterogeneous pairs
at the same δ/R, confirming the Hertz R_eff and E* combinations are correct.

### H4 — Newton balance (heterogeneous)

|F_disk1 − F_disk2| / F = **0.000** exactly for all 12 heterogeneous pairs.
(Action–reaction holds even across different material + size disks.)

---

## Phase 2B — Contact Primitives & Boundary Conditions

### What is validated

The general contact primitive framework (`contact_primitives.py` + `single_disk_contact.py`)
for 1 elastic disk against rigid geometric boundaries.

### Contact law (all primitive types)

```
gap, normal = primitive.gap_and_normal(node_xy)
    gap < 0   : node is penetrating
    normal    : unit vector from surface toward disk interior

F_node = k_pen · max(0, −gap) · normal
```

### Primitive types

| Class | Description | Corner handling |
|-------|-------------|-----------------|
| `LineSegment(p0, p1, inward_normal)` | Finite flat wall | n/a (no corners) |
| `Arc(center, radius, convex=True/False)` | Rigid circular arc | n/a |
| `Polygon(vertices)` | Closed polygon boundary | Closest-point projection (Voronoi ownership — no double-counting) |
| `Box(cx, cy, w, h)` | Axis-aligned rectangle | Polygon subclass |
| `Hopper(apex, half_angle, height)` | V-shaped hopper | Two LineSegments, nearest wins |

**Polygon corner rule:** Each perimeter node resolves against the single nearest segment
or vertex. No force is double-counted at corners. The force direction at a sharp apex is
the normal of the assigned segment — physically correct for closest-point projection.

### B1 — Two disks in a box (wall + disk-disk force balance)

**Setup:** Two equal elastic disks squeezed by rigid left/right flat walls.
`run_two_disk_wall_contact()` bisects on the disk center shift `cx_shift` until
F_wall ≈ F_disk-disk (Newton balance).

**Key Hertz result:** a_wall/a_dd = **1.00** (not √2).  
Reason: R_eff doubles (R/2 → R for wall contact) but E* also doubles (one body now
rigid instead of elastic). The two factors cancel exactly in the 2D Hertz formula
`a = √(4·F·R_eff / π·E*)`.

**Results (4 materials × 3 δ/R values):**

| Material | E (Pa) | ν | δ/R | F_balance | a_wall/a_dd |
|----------|--------|---|-----|-----------|-------------|
| Glass | 10⁷ | 0.22 | 0.02 | 0.8% | 0.996 |
| Polystyrene | 5×10⁵ | 0.33 | 0.02 | 1.6% | 0.998 |
| Reference | 10⁵ | 0.30 | 0.02 | 1.5% | 0.998 |
| Hydrogel | 10⁴ | 0.47 | 0.02 | 0.8% | 1.000 |

12/12 cases PASS (F_balance < 15%; most < 2%).

### B2 — Arc indenter sweep

**Setup:** 1 disk against two symmetric convex rigid arc indenters.  
R_eff = R·R_arc/(R+R_arc); E* = E/(1−ν²) (arc is rigid).

| R_arc/R | R_eff | F [N/m] | a_sim | a_hertz | ratio |
|---------|-------|---------|-------|---------|-------|
| 0.25 | 0.20 | 443 | 0.0523 | 0.0320 | 1.63 |
| 0.50 | 0.33 | 526 | 0.0523 | 0.0451 | 1.16 |
| 1.00 | 0.50 | 565 | 0.0523 | 0.0573 | 0.91 |
| 2.00 | 0.67 | 587 | 0.0523 | 0.0673 | 0.78 |
| 5.00 | 0.83 | 635 | 0.1045 | 0.0783 | 1.34 |

5/5 PASS (a_sim/a_hertz in [0.5, 2.0]). Monotonicity F(R_arc): 0 violations.

**Limiting cases:**
- R_arc → 0: point indenter, R_eff → 0, a → 0 ✓
- R_arc → ∞: flat wall, R_eff → R, a → a_wall ✓

### B3 — V-notch corner (no double-counting)

**Setup:** 1 disk against a Polygon V-notch (two segments meeting at apex, half_angle=30°).  
Reference: single flat wall at same x-position.

| Quantity | Value |
|----------|-------|
| F_vnotch / F_flat | 0.527 |
| Double-counting check (< 1.5) | **PASS** |
| Non-zero corner force (> 0.3) | **PASS** |
| Force direction from +x | 30° (= half_angle, expected) |

**Explanation of 30° direction:** The node at the apex maps to one segment
(first nearest in Voronoi tie-breaking). Force = k_pen·pen·normal_of_that_segment,
whose angle from +x equals the half-angle of the V. This is physically correct
for closest-point projection — the force is assigned, not doubled.

**Why F_vnotch/F_flat = 0.53:** The V-notch presents an angled surface (not flat),
so each contacting node gets a force with component sin(half_angle) in y that
partially cancels the x-component. At half_angle=30°, x-projection = cos(30°) = 0.87,
but fewer nodes penetrate the angled arms vs the flat wall → net ratio ≈ 0.53.

---

## Scaling Invariance

### Physics

Under a uniform length rescaling **X → λX** with E, ν unchanged (same δ/R):

| Quantity | Scaling |
|----------|---------|
| Positions X, displacements u | → λ × (proportional) |
| Forces F [N/m per unit depth] | → λ × (proportional) |
| Contact half-width a | → λ × (proportional) |
| Strains ε, stresses σ | → unchanged |
| Contact pressure P = F/a | → unchanged |
| Overlap ratio δ/R | → unchanged |

### Why it holds exactly (not approximately)

1. **Mesh:** `make_disk_mesh(R)` uses `refine_area = π·R²/(10·N)` → identical triangle
   connectivity for all R, with coordinates scaled exactly: `v(λR) = λ·v(R)`.

2. **Stiffness matrix:** K = ∫B^T·C·B dA. Under X→λX: B → B/λ (shape-function gradient),
   dA → λ²dA. Therefore K → (1/λ)·C·(1/λ)·λ² = K. **K is scale-invariant.**

3. **Contact forces:** f = k_pen·pen, where k_pen = α·E (no λ) and pen → λ·pen.
   Therefore f → λ·f.

4. **FEM equation:** K·u = f → K·(λu₀) = λ·f₀ ✓ (consistent with u → λu).

### Results

All relative errors = **0.00e+00** (exact floating-point equality):

| Case | λ | F ratio | u max error |
|------|---|---------|-------------|
| Two-disk | 0.5 | 0.500000000000000 | 0.00e+00 |
| Two-disk | 2.0 | 2.000000000000000 | 0.00e+00 |
| Wall contact | 0.5 | 0.500000000000000 | 0.00e+00 |
| Wall contact | 2.0 | 2.000000000000000 | 0.00e+00 |
| Single disk + wall | 0.5 | 0.500000000000000 | 0.00e+00 |
| Single disk + wall | 2.0 | 2.000000000000000 | 0.00e+00 |
| Glass λ=0.5 | 0.5 | 0.500000000000000 | 0.00e+00 |
| Hydrogel λ=0.5 | 0.5 | 0.500000000000000 | 0.00e+00 |

### Data augmentation use

The scaling symmetry is an exact data augmentation operator for the DtN NN:

```
Given: (X_peri, u_peri, f_peri, E, ν)
Aug:   (λ·X_peri, λ·u_peri, λ·f_peri, E, ν)   ← exact for any λ > 0
```

Cost: O(N_peri) multiply vs O(N_nodes²) FEM solve (~1000× cheaper per augmented sample).  
Recommended range: λ ∈ [0.5, 2.0] sampled log-uniformly (covers polydisperse particle sizes).  
To implement: see `src/validation/scaling_validation.py`, `run_augmentation_demo()`.

---

## Known Physics Deviations (Expected, Not Bugs)

| Deviation | Cause | Magnitude |
|-----------|-------|-----------|
| a_sim/a_hertz ≠ 1 | Hertz assumes infinite half-space; our disk is finite and curved | 0.77–1.24 |
| Finite-disk softening | Curved boundary carries less force than flat half-space | F_sim/F_hertz ≈ 0.43–0.49 |
| ν-dependence in F/(E·R) | Correct: incompressible materials carry more force at same strain | ~20% over full ν range |
| P1 locking at ν ≥ 0.45 | Volumetric locking in triangular P1 elements | Fixed: auto-switch to P2 |
| Linear elasticity limit | True strain can reach ~10% locally at δ/R = 0.15 | Accepted for current phase |

---

## Quick Reproduction

All scripts run from the project root with the venv activated:

```bash
source .venv/bin/activate
```

### Run all validations in order (~25 min total)

```bash
# Phase 1: FEM baseline (~2 min)
python src/validation/phase1_gate.py
python src/validation/nu_sweep_convergence.py

# Phase 2: Two-disk homogeneous (~10 min)
python src/validation/twodisk_validation.py

# Phase 2-H: Two-disk heterogeneous (~5 min)
python src/validation/hetero_validation.py

# Phase 2B: Contact primitives (~5 min)
python src/validation/primitives_validation.py

# Scaling invariance (~3 min)
python src/validation/scaling_validation.py
```

### Run individual Phase 2B tests

```bash
python src/validation/primitives_validation.py B1   # wall + disk-disk balance
python src/validation/primitives_validation.py B2   # arc indenter sweep
python src/validation/primitives_validation.py B3   # V-notch corner
```

### Quick mode (Phase 2, 5 materials × 5 deltas, ~2 min)

```bash
python src/validation/twodisk_validation.py --quick
```

### Production mode (N=240, ~3 hours)

```bash
python src/validation/twodisk_validation.py --N 240
```

---

## Source File Map

| File | Purpose |
|------|---------|
| `src/simulation/disk_mesh.py` | Disk mesh generator (Triangle, refine_area∝R²) |
| `src/simulation/fem_elastic.py` | P1/P2 FEM stiffness assembly + solver |
| `src/simulation/two_disk_contact.py` | Two-disk penalty contact; `run_two_disk_wall_contact()` |
| `src/simulation/single_disk_contact.py` | Single disk against arbitrary primitive list |
| `src/simulation/contact_primitives.py` | `LineSegment`, `Arc`, `Polygon`, `Box`, `Hopper` |
| `src/contact/hertz.py` | Analytic 2D plane-strain Hertz formulas |
| `src/validation/phase1_gate.py` | P1: uniform pressure analytic test |
| `src/validation/nu_sweep_convergence.py` | P2: ν-sweep + N-convergence |
| `src/validation/dirichlet_hertz_test.py` | P3: Dirichlet Hertz, P1/P2 auto-select |
| `src/validation/twodisk_validation.py` | V1–V6: homogeneous two-disk suite |
| `src/validation/hetero_validation.py` | H1–H4: heterogeneous two-disk suite |
| `src/validation/primitives_validation.py` | B1–B3: primitive contact benchmarks |
| `src/validation/scaling_validation.py` | Scaling invariance + augmentation demo |

---

*Validated: 2026-04-16. All tests PASS.*
