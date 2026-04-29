# VISION.md — Science Question, Long-Term Goals, and Design Philosophy

> This is the anchor document. All other files (PLAN.md, CLAUDE.md, HANDOFF.md) point
> here for the *why*. Read this first when onboarding, after a long break, or when
> re-evaluating a design decision.
>
> This document evolves as the science evolves. Update it when the questions sharpen,
> not when the implementation details change.

---

## Model Philosophy and Academic Framing

The capsule shell DEM is a **deliberately constructed toy model**. It is not a
first-principles simulation of any specific physical material. It is not derived from
continuum elasticity, emulsion physics, or foam mechanics, and it makes no claim to
reproduce the quantitative behavior of any of those systems.

The model is designed to have **tunable macroscopic behavioral features** — specifically,
the ratio of shape deformation to area change under contact loading — so that these
features can be varied systematically and their effect on collective dynamics studied in
isolation.

**What this means for how results are communicated:**

Analogies to real materials (rubber balls, glass beads, emulsion droplets) are used only
in the loosest qualitative sense: *"particles whose contact behavior is reminiscent of..."*
or *"a regime that shares features with..."*. These analogies are illustrative, not claims
of faithful simulation. No referee should be armed to reject results on the grounds that
the model fails to reproduce, say, the exact contact mechanics of a 2D elastic disk with
ν = 0.3 — because that is not what the model claims to do.

**The squishiness knob is defined operationally, not derivationally.** It is characterized
by four measurable behavioral metrics (ν_meas, ΔA, contact flatness, penetration) on a
two-disk squeeze test. The physics it interpolates between are described in terms of what
the particles *do*, not what material they are *made of*.

**The science is about emergence.** As the squishiness knob is varied, certain collective
behaviors may appear or disappear — changes in rearrangement statistics, stress
propagation, packing geometry, jamming density. These emergent phenomena are the
scientific contribution. The capsule model parameters are the means, not the end.

This philosophy protects the work from two common failure modes:

1. **Overclaiming** — *"our model correctly simulates rubber"* — vulnerable to any
   quantitative discrepancy with rubber measurements.
2. **Underclaiming** — *"this is a purely abstract model with no physical relevance"* —
   fails to communicate why the results matter.

The correct framing is in between: a toy model with physically motivated behavioral
limits, used to study the role of a specific physical distinction (shape-change
vs. area-change) in collective granular dynamics.

---

## The Science Question

How does particle deformability — the ability of a disk to elastically deform under
contact — influence collective dynamics, rearrangements, and stress relaxation in dense
2D packings near the jamming transition φ_J?

This breaks into several sub-questions, each representing a distinct physical mechanism:

### 1. Volume-preserving deformation and pore filling
An elastic particle under contact can maintain constant volume while redistributing its
shape — bulging laterally to fill the pore space around it. Does this shape change
(distinct from simple Hertz force) alter how particles lock and unlock during
rearrangements?

### 2. Steric-like friction from deformation
When a particle deforms, its perimeter changes shape. This non-circular perimeter can
create geometric interlocking with neighbors — a deformation-induced steric effect —
even without surface friction. Does this mechanism contribute to or inhibit
rearrangements, and how strongly does it scale with squishiness?

### 3. Compressibility and volume loss
A compressible particle (low ν) contracts under isotropic pressure, actually losing
volume. A near-incompressible particle (high ν) cannot — it deforms at constant
volume. Does this distinction matter for local packing geometry and stress transmission?
How does the stress field differ between a jammed system of compressible vs.
incompressible particles?

### 4. Comparison to classical DEM
In classical rigid-body DEM, particles are treated as rigid with overlap-based
repulsion. Soft particles that overlap and slide through each other approximate this
limit. The question is: at what squishiness does the deformable-particle physics become
qualitatively distinct from classical DEM? Where do the two models agree, and where do
they diverge?

### 5. Stress propagators near rearrangements
When a particle rearranges (a T1-like event), the resulting stress pulse propagates
outward through the packing. How does the spatial structure and decay of this
propagator depend on E and ν? Do soft particles act as stress absorbers, damping
propagation, or do they create longer-ranged correlations?

---

## Current Goal

Build and validate the **EPD (Elastic Particle DEM)** capsule shell model — a pure DEM
approach where each disk is a closed elastic membrane with a tunable squishiness knob ν.

1. Model spans the full squishiness spectrum from glass-like (stiff, compressible) to
   rubber-like (soft, area-conserving) with a single user-level parameter ν
2. Validated against Hertzian contact, capillary-wave damping, T1-events, and falling droplet
3. Documented in a methods paper (`papers/summary_of_methods/`)
4. Enables multi-particle studies of how squishiness affects collective dynamics

**Note — axed path:** A FEM/NN-DtN surrogate approach (replacing FEM interior solves with
a neural network) was explored and explicitly abandoned. The EPD capsule shell model is the
chosen path — computationally cheap, physically motivated, and fully validated.

---

## The Squishiness Spectrum

The control variable is **squishiness** — a continuous 1D axis running from rigid/stiff
to very soft/compliant. The real materials below are **guides**, not discrete choices.
Training data, parameter sweeps, and science results should be presented as continuous
functions along this axis.

| Guide material | E range (Pa)  | ν range    | Key character |
|----------------|---------------|------------|---------------|
| Glass beads    | 10⁷ – 10⁸    | 0.18–0.25  | Stiff, compressible, Hertz-like |
| Polystyrene    | 10⁵ – 10⁶    | 0.30–0.35  | Intermediate stiffness |
| PDMS gels      | 10³ – 10⁵    | 0.45–0.48  | Soft, nearly incompressible |
| Rubber/Jello   | 10² – 10⁴    | 0.48–0.49  | Very soft, volume-preserving |

**Note on (E, ν) coupling**: In physical materials, E and ν are not independent.
Rubber-like materials have both low E *and* high ν. It is not physically meaningful
to study (low E, low ν) or (high E, high ν) as real material regimes — though they
may appear as limiting cases in parameter sensitivity tests.

The hardest computational regime — low E, high ν, large strain — is also the least
physically natural. The focus of this study is the physically motivated squishiness
axis, not the full E×ν hypercube.

---

## Physical Limits and Behavioral Targets

### The Core Physical Distinction

The central axis of this study is not "stiff vs. soft" in the naive sense of force
magnitude.  It is about **how stress is absorbed at the particle level**:

> **Squishy (rubber-like):** stress is absorbed through *shape change at constant area*.
> The particle deforms dramatically — flattening at contacts, bulging laterally into pore
> space — but its cross-sectional area is conserved.  All deformation energy goes into
> shear/bending modes.
>
> **Stiff (glass-like):** stress is absorbed through *area change at near-constant shape*.
> The particle barely changes shape but contracts slightly and uniformly under load, and
> the contact zone flattens only over a narrow Hertz-like patch.  Deformation energy goes
> primarily into volumetric compression.

This distinction is what makes squishiness scientifically interesting for jamming: a
rubber-like particle fills pore space (lubricating rearrangements), while a glass-like
particle transmits force chains without redistributing shape.  The DEM model must
faithfully reproduce both limits and everything in between.

---

### Four Operational Metrics

These four metrics are the **north star** for capsule model calibration.  Every parameter
choice is evaluated against them.  They are measured on a two-disk squeeze test with
**top and bottom walls only** (no side walls), so lateral expansion is unconstrained.

| # | Metric | Formula | What it measures |
|---|--------|---------|-----------------|
| 1 | **ν_meas** | ε_lateral / ε_vertical | Fraction of deformation that goes into shape change |
| 2 | **ΔA** | 1 − A/A₀ | Fraction of deformation that goes into area change |
| 3 | **Contact flatness** | Mean turning angle at active contact nodes → 0° | Whether the interface flattens under load |
| 4 | **Penetration** | \|cc_gap\| / r_c | Numerical artifact; controlled by C |

**Geometry note (2D):** ε_lateral is measured at the equator nodes (0° and 180°).
ε_vertical is the per-particle strain already defined (nodes at 90° and 270°).
In 2D, perfect area conservation requires ε_lateral = ε_vertical, giving ν_meas = 1.0
(not 0.5 as in 3D continuum elasticity — that is a 3D number).  ΔA ≈ (1 − ν_meas) × ε_v
to first order, so the two metrics are not independent; both are reported for clarity.

**Contact flatness geometry:** At each node *i* actively in contact, the turning angle
is 180° − (interior angle between edges i−1→i and i→i+1).  For an undeformed N=32
circle this is ~11.25°.  A perfectly flat contact patch has turning angle = 0°.  The
metric is the mean turning angle over all active contact nodes; it decreases toward 0°
as τ decreases and the interface flattens.

---

### Target Values at the Two Limits

| Metric | Soft/rubber limit | Stiff/glass limit | Notes |
|--------|------------------|------------------|-------|
| ν_meas | ≥ 0.8 | ≤ 0.2 | Rubber: shape absorbs load; glass: shape is retained |
| ΔA at ε_v=12% | < 2% | ~15–30% of ε_v | Rubber area-conserving; glass contracts measurably |
| Contact flatness (turning angle) | → 0° (wide flat patch) | → 0° (narrow flat patch) | Both limits should flatten locally; width differs |
| Penetration \|cc_gap\|/r_c | < 10% | < 10% | Numerical target, same for both |

Note: both limits should produce a **locally flat** contact zone (turning angle → 0°) —
this is a physical requirement, not a limit-specific one.  The difference is the *width*
of the flat patch: rubber has many nodes at ~0°, glass has very few.

---

### Control Parameters and Their Roles

| Parameter | Symbol | Primary effect | Secondary effects |
|-----------|--------|---------------|-------------------|
| Edge spring stiffness | S | Shape stiffness (ν_meas, resistance to deformation) | Sets force scale |
| Area stiffness ratio | K_A / S | Area conservation vs. compressibility (ΔA, ν_meas) | Must be coupled to S — see below |
| Bending stiffness | τ | Contact flatness (turning angle at contact) | Also affects global shape at large strain |
| Contact penalty | C | Penetration depth | Interacts with flatness: flatter contact → force distributed → less penetration |

**Critical coupling — K_A must scale with S:**
Setting K_A (the fluid pressure coefficient, `rho_f` in the code) independently of S
breaks both limits.  The physically motivated coupling draws from continuum elasticity:

```
K_A / S  ≈  K_bulk / G  =  2(1+ν) / (3(1−2ν))
```

- ν → 0.49 (rubber): K_A/S ≈ 50  →  area 50× harder to change than shape
- ν → 0.20 (glass):  K_A/S ≈ 1.3  →  area and shape comparably stiff

This ratio is the single "squishiness knob" at the capsule level.  In practice, the
calibration procedure sweeps K_A at fixed edge-case S values to empirically find the
K_A(S) coupling that hits the target (ν_meas, ΔA) values above, rather than
pre-committing to the continuum formula.

---

### Calibration Philosophy

1. **Edge cases first.** Characterize the two limits (S=small, S=large) independently.
   Find K_A and τ values that hit the four metric targets at each limit.
2. **Interpolate.** With anchor points established, define a smooth K_A(S) and τ(S)
   coupling that interpolates between the limits.  This defines the "squishiness axis"
   in capsule parameter space.
3. **Validate continuously.** At each point along the interpolated axis, confirm that
   all four metrics evolve smoothly and monotonically between the two limits.
4. **Only then sweep.** Use the calibrated, physically consistent parameter set for the
   full multi-particle simulations.

---

## EPD Model Architecture

The EPD model represents each particle as a 2D closed elastic membrane (capsule shell)
with N perimeter nodes. No FEM interior solve is needed — the membrane forces act
directly on the perimeter nodes, and rigid-body dynamics (x_cm, v_cm, θ, ω) are
integrated explicitly.

### Per-timestep loop
1. Contact detection: find pairs whose perimeters overlap
2. Contact penalty forces applied to overlapping perimeter nodes
3. Membrane forces (edge springs, bending, area spring) computed from node positions
4. Net force and torque integrated → rigid body step (x_cm, θ updated)
5. Node positions reconstructed from updated rigid-body frame

No NN, no FEM, no interior DOFs. The model is O(N_particles × N_nodes) per timestep.

---

## Perimeter Node Resolution

The perimeter is the fundamental resolution unit. Choices:

| Tier | N_nodes | Angular spacing | Use |
|------|---------|-----------------|-----|
| Development | 120 | 3°/node | Solver verification, early NN prototyping |
| Production | 240 | 1.5°/node | NN training data generation |
| Benchmark | 360 | 1°/node | High-fidelity reference, NN validation |

**Why these numbers**: the contact arc at typical working strains (δ/R ≈ 1–5%) spans
roughly 10–40°. With 3°/node spacing (N=32 at R=1), a 1% strain contact has ~3 nodes —
marginal but usable for calibration. With 1.5°/node (N=240), the same contact has ~6 nodes
— well resolved. Near-jamming contacts (δ/R → 0) approach zero nodes; this is acceptable
as long as the repulsion is smooth and forces go to zero continuously.

---

## Key Assumptions and Their Scope

| Assumption | Status | Scope / Caveat |
|------------|--------|----------------|
| Linear (small-strain) elasticity | **Locked** | Valid for ε < ~5%; adequate near jamming where strains are small. For large-compression rubber, this is approximate. Non-linear extension is a natural future step — the data pipeline and NN architecture are unchanged; only the FEM solver is swapped. |
| 2D plane strain | **Locked** | Tractable first step; captures essential jamming physics. 3D extension is natural but much more expensive. |
| Quasi-static elastic DOFs | **Locked** | Rearrangement timescales >> elastic wave propagation. Valid near jamming; breaks down for fast dynamics. |
| Rigid body dynamics | **Locked** | Full elastodynamics of interior is unnecessary for collective behavior. |
| Circular reference shape | **Locked** | Disks. Polydispersity (different R) is straightforward. Non-circular particles are a separate project. |
| Analytic Hertz ExtForce | **Phase 1–3** | Reassessed in Phase 4; may need learned ExtForce for very soft particles at large overlap. |

---

## What We Expect to Learn as We Go

This document reflects current understanding. The following are expected to change:

- **Required N_nodes**: The minimum perimeter resolution for the DtN NN to generalize
  accurately across (E, ν) is not yet known. May end up needing N=480 or N=720.
- **NN architecture**: Whether a simple MLP, an equivariant architecture, or a neural
  operator (FNO/DeepONet) is needed for 10⁻⁶ precision is an open question.
- **Training data volume**: 50,000 samples is a guess; could be 10× less or more.
- **Linear elasticity limits**: For jello-like particles at φ > 0.90, the linear
  approximation may introduce visible artifacts. We'll know after Phase 4 validation.
- **Science questions**: Sub-questions 1–5 above are hypotheses. Some may turn out to
  be trivial; others may open new directions.

---

## Relation to Existing Files

| File | Role |
|------|------|
| **VISION.md** (this file) | Science question, long-term goal, design philosophy |
| **PLAN.md** | Execution roadmap — phases, waypoints, gate metrics |
| **CLAUDE.md** | Technical instructions for the AI agent |
| **HANDOFF.md** | Session state — where we left off, immediate next step |
| **LOG.md** | Append-only record of what was done and what was measured |
| **REPRODUCIBILITY.md** | How to reproduce key results (created when a branch is locked) |

---

*Last updated: 2026-04-18 — Added "Physical Limits and Behavioral Targets" section: four metrics (ν_meas, ΔA, contact flatness, penetration), target values at rubber/glass limits, K_A/S coupling rationale, calibration philosophy*
