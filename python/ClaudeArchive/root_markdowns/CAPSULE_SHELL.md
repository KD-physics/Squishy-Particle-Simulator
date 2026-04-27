# CAPSULE_SHELL.md — Fluid-Filled Elastic Shell Particle Model

> **Branch:** Phase 5 — Capsule Shell DEM
> **Status:** Design complete, implementation pending
> **Supersedes:** FEM/DtN/NN pipeline for the DEM simulation layer (Phases 3–4 on hold)

---

## 1. Motivation

The FEM interior solve does not scale to N ~ 1000 particles:
- Global K matrix: O((N · N_interior)³) — explodes with particle count
- Per-disk K: O(N_interior²) per timestep — still ~8 ms/particle at N=120
- Penalty contact iteration: 10–120 FEM solves to converge a single strain step
- N=240 solver: non-convergent above δ/R ≈ 0.02 (hits max_iter=120)

The Capsule Shell model replaces the FEM interior entirely. Each particle is a
**fluid-filled elastic shell** discretized as N perimeter nodes. All forces are
computed from the perimeter alone — no interior nodes, no K matrix, no global solve.
Each timestep is explicit Verlet: O(N_particles × N_nodes) with a small constant.

---

## 2. Physical Picture

Each particle consists of:
- **N perimeter nodes** x_i (i = 0…N-1), arranged counterclockwise
- **Fluid interior**: incompressible (or weakly compressible) fill that equilibrates
  pressure instantaneously across the particle
- **Elastic shell**: thin membrane of thickness t and modulus E_l connecting the nodes

The two mechanisms by which a particle absorbs contact stress:
1. **Area compression** (bulk mode): fluid compresses, total area decreases
2. **Shape deformation** (bending mode): shell bends, contact zone flattens,
   free surface bulges at constant perimeter

The relative ease of these two modes is set by the dimensionless squishiness
parameter S (Section 4). Physical particles span from glass (S >> 1, nearly
rigid) to jello (S << 1, floppy).

---

## 3. Geometry

### 3.1 Nodes and Edges

| Symbol | Definition |
|--------|-----------|
| x_i | Position of node i (2D vector) |
| X_i | Reference position (circle of radius R₀, equally spaced) |
| e_i | Edge vector: x_{i+1} − x_i |
| L_i | Current edge length: \|e_i\| |
| L₀ | Reference edge length: 2πR₀/N |
| t̂_i | Unit tangent: e_i / L_i |
| n̂_i | Outward edge normal: rot90_CCW(t̂_i) |
| n̂_i^node | Node outward normal: normalize(n̂_{i-1} + n̂_i) |
| θ_i | Signed turning angle at node i: angle from t̂_{i-1} to t̂_i (CCW positive) |
| θ₀_i | Reference turning angle at node i (assigned at initialization; need not be 2π/N for non-circular shapes) |
| A | Current enclosed area: shoelace formula |
| A₀ | Reference enclosed area: πR₀² |

### 3.2 Capsule Geometry

Each edge i is a **capsule**: a line segment of length L₀ with hemispherical
end caps of radius r_c = L₀/2. Two capsules overlap when the distance between
their line segments is less than 2·r_c = L₀.

---

## 4. Dimensionless Parameters

### 4.1 Reference Scales

Set once for all particles and never varied:

| Quantity | Value | Meaning |
|----------|-------|---------|
| R₀ | 1 | Mean particle radius (length unit) |
| K_fluid | 1 | Fluid bulk modulus (force/pressure unit) |
| A₀ | π | Reference area (πR₀²) |
| L₀ | 2π/N | Reference edge length |

### 4.2 Physical Parameters

| Parameter | Symbol | Definition | Range of interest |
|-----------|--------|-----------|-----------------|
| Shell thickness ratio | τ | t / R₀ | 0.01 – 0.3 |
| Membrane stiffness ratio | β | E_l · t / (K_fluid · R₀) | 10 – 10000 |
| Squishiness | **S** | E_l · t³ / (12 · K_fluid · R₀⁴) | 0.01 – 100 |
| Contact hardness | C | k_c · R₀ / K_fluid | 100 – 2000 |
| Fluid density | ρ_f | (sets mass and timescale) | simulation-dependent |

### 4.3 Derived Relations

```
S = β · τ² / 12          (squishiness from membrane ratio and thickness)
β / S = 12 / τ²          (edge stiffness always >> bending for thin shells)
EI = E_l · t³ / 12 = S · K_fluid · R₀⁴    (bending stiffness)
k_c = C · K_fluid / R₀   (contact spring, material-independent)
m_i = ρ_f · t · L₀       (node mass)
```

### 4.4 The Squishiness Axis

| Regime | S | τ | Physical analog |
|--------|---|---|----------------|
| Nearly rigid | >> 1 | moderate | Glass: barely deforms under contact |
| Stiff shell | 10–100 | 0.1–0.3 | Polystyrene: slight contact flattening |
| Intermediate | 1–10 | 0.05–0.1 | PDMS: visible flattening, moderate bulge |
| Floppy shell | 0.1–1 | 0.02–0.05 | Rubber: strong shape deformation |
| Very floppy | << 1 | 0.01–0.02 | Jello: large deformation, near area-conserving |

### 4.5 The Unphysical Scaling Mode and How τ Prevents It

For extensible shells (small β), a particle under contact can uniformly shrink
(radius scales down) without bending — zero bending energy, no shape change.
This is unphysical: the contact zone should flatten, not the whole particle shrink.

The fix is encoded in τ: for small τ (thin shell):
```
β / S = 12 / τ²  →  ∞  as τ → 0
```
Edge stretching becomes infinitely more expensive than bending. The only
low-energy response to contact is bending (contact flattens), not uniform scaling.
**τ is therefore a physics hypothesis** — choosing τ determines how strongly
inextensibility is enforced. It must be validated (see Section 7).

---

## 5. Force Inventory

All forces are computed per-node. No matrix assembly, no global solve.

### 5.1 Fluid Pressure (global area, normal, instantaneous)

```
A = (1/2) |Σᵢ (xᵢ × x_{i+1})|          # shoelace area
P = K_fluid · (A − A₀) / A₀             # scalar pressure (same for all nodes)
F_press_i = P · L₀ · n̂_i^node           # outward normal force at node i
```

Physical meaning: fluid pressure equilibrates instantly across the particle.
All nodes feel the same P simultaneously — the area signal propagates at infinite
speed (mean-field coupling), appropriate for a fluid fill.

### 5.2 Hydrostatic Gravity (optional)

```
y_top = max(yᵢ)                          # highest node
h_i = y_top − y_i                        # fluid column above node i
P_hydro_i = ρ_f · g · h_i
F_grav_i = P_hydro_i · L₀ · n̂_i^node
F_grav_i −= mean_x(F_grav) · x̂          # regularize: zero net x-force
```

Result: Σ F_y = −ρ_f · g · A (particle weight downward), Σ F_x = 0.
This is Archimedes' principle applied to the enclosed fluid.

### 5.3 Edge Elasticity (local, tangential)

For each edge i (connecting nodes i and i+1):
```
ε_i = (L_i − L₀) / L₀                   # axial strain
F_edge = E_l · t · εᵢ · t̂ᵢ
F_i   += F_edge       (node i pulled toward i+1 if stretched)
F_{i+1} −= F_edge     (node i+1 pulled toward i if stretched)
```

Physical meaning: elastic 1D membrane resists perimeter stretching/compression.
Together with fluid pressure, this prevents the unphysical uniform-scaling mode
when τ is chosen appropriately.

### 5.4 Bending (local, 3-node hinge, energy-consistent)

Bending stiffness: EI = E_l · t³ / 12

For each node i (acting as a hinge between edges i-1 and i):
```
θᵢ = atan2(t̂_{i-1} × t̂ᵢ, t̂_{i-1} · t̂ᵢ)   # signed turning angle (2D)
Mᵢ = (EI / L₀) · (θᵢ − θ₀ᵢ)                # bending moment at hinge i

F_{i-1} += Mᵢ / L₀ · n̂_{i-1}               # outward normal to left edge
Fᵢ      −= Mᵢ / L₀ · (n̂_{i-1} + n̂ᵢ)       # both edges: sum
F_{i+1} += Mᵢ / L₀ · n̂ᵢ                    # outward normal to right edge
```

This is the exact gradient of the bending energy E_bend = (EI/2L₀) Σᵢ (θᵢ − θ₀ᵢ)².
Node i participates in three hinges; the total bending force on i accumulates
contributions from hinges i-1, i, and i+1.

**Reference angles θ₀ᵢ** are set at particle initialization from the initial node
positions — not forced to be 2π/N. Non-circular reference shapes are supported.

### 5.5 Capsule–Capsule Contact (particle–particle, 2-point Gauss quadrature)

Each edge is a capsule of radius r_c = L₀/2.
Contact radius for two capsules touching: 2·r_c = L₀.
Contact stiffness: k_c = C · K_fluid / R₀.

For each pair of edges (i→i+1) on particle A and (j→j+1) on particle B:
```
# 2-point Gauss quadrature on edge A: s ∈ {s₁, s₂}, weights w = 0.5 each
# s₁ = (1 − 1/√3)/2,  s₂ = (1 + 1/√3)/2

for each quadrature point at parameter s:
    x_q  = (1−s) · xᵢ + s · x_{i+1}           # point on edge A
    y_q  = closest_point_on_segment(x_q, xⱼ, x_{j+1})  # closest on edge B
    d    = |x_q − y_q|
    gap  = d − 2·r_c
    if gap < 0:
        n̂ = (x_q − y_q) / d                    # normal from B toward A
        F  = k_c · (−gap) · n̂ · L₀ · w        # force magnitude
        # Distribute to A nodes:
        Fᵢ      += (1−s) · F
        F_{i+1} +=    s  · F
        # Distribute to B nodes (equal and opposite, via B's shape function):
        Fⱼ      −= (1−s_B) · F
        F_{j+1} −=    s_B  · F
        # where s_B = parameter of y_q along edge j
```

### 5.6 Capsule–Primitive Contact (particle–wall/arc/corner)

All primitives (LineSegment, Arc, Polygon, V-notch) expose:
```
primitive.gap_and_normal(x_q) → (gap, n̂)
```
The quadrature is identical to 5.5, but the primitive side is rigid (no force
distribution to primitive nodes). The existing `contact_primitives.py` is
directly reusable — the interface is unchanged.

---

## 6. Time Integration

Pure explicit Verlet on perimeter nodes. No matrix, no global solve.

```
# Per timestep:
for each particle:
    accumulate F_i from forces 5.1–5.6
    a_i = F_i / m_i          # m_i = ρ_f · t · L₀
    v_i += a_i · dt
    x_i += v_i · dt
```

Optionally: add a damping term α_damp · v_i to F_i to reach quasi-static limit.

---

## 7. Stability and Timestep

Three independent timescales constrain dt:

| Mode | Frequency | Constraint on dt |
|------|-----------|-----------------|
| Breathing (area) | ω_breath ~ √(2K_fluid / (ρ_f · t · R₀)) | dt < 2/ω_breath |
| Edge wave | c_edge ~ √(E_l / (ρ_f · L₀)) | dt < L₀/c_edge |
| Contact | ω_contact ~ √(k_c / m_i) | dt < 2/ω_contact |

All three are determined by the parameters already in the model. The stable dt is:
```
dt_max = 2 / max(ω_breath, N · c_edge / R₀, ω_contact)
```

**Open questions (to be measured, not assumed):**
- What is the actual stable dt for S ∈ [0.01, 100] and τ ∈ [0.01, 0.3]?
- Does β (edge stiffness) or breathing mode set the stability limit?
- What damping coefficient α_damp is needed for quasi-static behavior, and
  does the quasi-static limit agree with the FEM two-disk results?
- What value of C (contact hardness) gives accurate Hertz-like contact without
  requiring impractically small dt?

---

## 8. Connection to Existing Infrastructure

| Existing component | Reuse |
|-------------------|-------|
| `src/simulation/contact_primitives.py` | Direct reuse — gap_and_normal() interface unchanged |
| `src/validation/twodisk_validation.py` | Reference target for quasi-static limit validation |
| `src/validation/primitives_validation.py` | Reference for wall/arc/corner contact validation |
| `src/simulation/two_disk_contact.py` | FEM ground truth for force comparison |
| RCPGenerator | Direct reuse for N-body initial packings |

The two-disk test from Phase 2 becomes the **primary validation benchmark**:
run the capsule-shell two-disk simulation to the quasi-static limit (high
damping, small dt) and compare contact force vs δ/R against the FEM results.
This tests: symmetry, force balance, Hertz-like scaling, primitive contact.

---

## 9. Open Hypotheses

| Hypothesis | How to test |
|-----------|------------|
| β = 12/τ² enforces inextensibility (prevents uniform scaling) | Two-disk test: check contact zone flattens (not uniform shrink) vs β sweep |
| Quasi-static limit = FEM two-disk result | Compare F vs δ/R at high damping to `twodisk_validation.py` values |
| τ ~ 0.05–0.1 gives physical Hertz-like pressure distribution | Measure contact arc flattening vs τ |
| C ~ 500 gives accurate contact without excessive stiffness | Measure F vs δ/R vs C; find insensitive plateau |
| Elastic wave speed agrees with analytic c_edge | Measure pulse propagation time around shell |

---

*Created: 2026-04-17 — Design session, architectural pivot from FEM/NN to Capsule Shell DEM*
