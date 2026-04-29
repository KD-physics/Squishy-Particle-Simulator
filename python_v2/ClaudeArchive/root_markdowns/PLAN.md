# PLAN.md — Execution Roadmap

> **Read VISION.md first** — it defines the science question, model philosophy, and
> behavioral targets that every phase works toward.
> This file is the execution layer: phases, waypoints, and gate metrics.
>
> Waypoints are gates. Do not advance until all gate metrics pass.
> Record measured values next to target metrics as you go.
> Check off waypoints with ✅ when complete, ❌ if blocked, ⬜ pending.

---

## Overall Progress

| Phase | Title | Status |
|-------|-------|--------|
| **0** | Hertz/FEM Infrastructure (scikit-fem) | ✅ COMPLETE |
| **1** | Capsule Shell Core Engine + Two-Disk Validation | ✅ COMPLETE |
| **1B** | Squishiness Axis Calibration | ✅ COMPLETE — q_glass=0.1, q_rubber=10, τ_flat=0.20 |
| **1C** | Dimensional Scaling Verification | ✅ COMPLETE — S=force scale, R0=length scale; EI bug fixed |
| **1D** | C Formula Calibration | ✅ COMPLETE — C=3000×(1+q), τ-independent; cc_pen<2% all (τ,q) |
| **1E** | Rigid-Body Decomposition Integrator | ✅ COMPLETE — step_rb() benchmarked, calibrated |
| **1F** | Collision Physics Validation Suite | ✅ COMPLETE — all gates pass; COR calibration table locked |
| **1G** | Extension Benchmarks (3-body + squeeze) | ✅ COMPLETE — all gates pass; machine-precision momentum in 3-body; T1 squeeze stable |
| **2** | Emulsion Droplet Model | ✅ COMPLETE — all benchmarks pass, movies rendered, §9 in paper |
| **3.1** | TF Intra-Particle Kernels + Integration | ✅ COMPLETE — float64, max_diff ≤ 8.88e-16 |
| **3.2** | Python Candidacy Manager + Inter-Capsule Forces | ✅ COMPLETE — 3-level filter, forces 4.55e-13, traj 1.78e-15 |
| **3.3** | Inter-Capsule Forces on TF, Full Graph Assembled | ✅ COMPLETE — step_full_tf validated, 500-step traj 9.66e-14 |
| **3.4** | C++ Candidacy Manager + pybind11 | ✅ COMPLETE — pip install . works, corpus 18/18, e2e Δx=0 |
| **3.5** | Fast Handshaking + Scale-Up (future) | ⬜ PENDING |
| **4.0** | High-Level API — Design Document & Forward Graph | ✅ COMPLETE — 15/15 patch tests pass |
| **4.1** | `objects.py` — Object hierarchy, resolved(), region_polygon() | ✅ COMPLETE — 32/32 PASS |
| **4.2** | `particles.py` — ParticleSpec, mixed types, material flexibility | ✅ COMPLETE — 29/29 PASS |
| **4.3** | `motion.py` — MotionSpec, TF-native + sampled, time-varying | ✅ COMPLETE — 24/24 PASS |
| **4.4** | `initializer.py` — RSA seeder + adaptive swell, exclusion-aware | ✅ COMPLETE — 14/14 PASS |
| **4.5** | `system.py` + `checkpoint.py` — System class, save/load, full state | ✅ COMPLETE — 9/9 PASS |
| **4.6** | Integration verification — all channels tested, GIF output per case | ✅ CLOSED — valuable groundwork, tests incomplete/visual issues; superseded by 4.7 |
| **4.7** | Systematic test validation → Colab notebook assembly | ✅ COMPLETE — A/B/C/D/F tests pass; Colab notebook built |
| **4.8** | Per-particle R0 scaling, parameter recompute, TF graph fix | ✅ COMPLETE |
| **6** | GitHub Release Package + Reproducibility Notebook | ✅ COMPLETE |
| **Test F** | Emulsion hopper discharge (gravity + Poiseuille) | ✅ COMPLETE — 120k steps, 14 exited, NaN-free |
| ~~3~~ | ~~DtN NN Training Data (FEM)~~ | ⛔ AXED — FEM/NN approach abandoned |
| ~~4~~ | ~~DtN Neural Network Surrogate~~ | ⛔ AXED — FEM/NN approach abandoned |

---

## Phase 0: FEM Infrastructure ✅ COMPLETE

Preserved as ground truth. All validation scripts pass. Do not modify.

Key results: two-disk F vs δ/R curves; contact primitives (LineSegment, Arc, Polygon);
force balance validated to machine precision.

---

## Phase 1: Capsule Shell Core Engine ✅ COMPLETE

Core deliverables:
- `src/simulation/capsule_shell.py` — CapsuleParticle, CapsuleSim, all force terms
- `src/validation/twodisk_capsule.py` — two-disk squeeze, movies, metrics
- `src/validation/stage1_sweep.py` — systematic C/S/SR/alpha0 parameter survey

Key findings:
- C = 100 sufficient for all S ∈ [0.005, 10] at N=32 (pen_frac < 1.1%)
- Steps per run constant regardless of S (~15,000 at N=32, ~70s)
- alpha_damp = alpha0 / T_wave keeps damping ratio ζ constant across S
- Per-particle strain eps_p ≈ 0.25–0.46 × wall strain (wall overestimates compression)

**What Phase 1 revealed that drives Phase 1B:**

The original S parameter alone does not produce the desired behavioral limits.
Specifically:
- Small S: large shape deformation but also uncontrolled area loss (not rubber-like)
- Large S: rigid body behavior with no contact flattening (not glass-like)
The missing ingredient: K_A (area stiffness, `rho_f`) must be coupled to S.
The bending stiffness τ must be coupled to S to control contact flatness independently.
This led to the four-metric framework and the squishiness axis (see VISION.md).

---

## Phase 1B: Squishiness Axis Calibration

**Objective:** Define a single squishiness knob q = K_A / S that produces the correct
behavioral limits at both ends, with τ/S as a secondary knob for contact flatness.
Characterize all four metrics (ν_meas, ΔA, contact flatness, penetration) as functions
of (q, τ/S) at edge-case S values. Validate that a sweep along q gives a smooth,
monotone interpolation between the two physical limits.

**Setup change from Phase 1:** No side walls. Two disks stacked vertically, top and
bottom walls only. This allows free lateral expansion, enabling measurement of ν_meas.

**Pre-numerical expectations (to be confirmed or revised):**
- ν_meas and ΔA are related: ΔA ≈ (1 − ν_meas) × ε_v in 2D → one knob controls both
- ν_meas increases monotonically with q = K_A/S (high q → incompressible rubber limit)
- Contact flatness (turning angle → 0°) controlled by τ/S; approximately independent of q
- Penetration controlled by C; secondary benefit from flatter contact

---

### Waypoint 1B.1 — New measurement infrastructure ✅

Add to `twodisk_capsule.py`:

1. **Lateral strain** ε_l = (x[N//2, 0] − x[0, 0]) / (2 R₀) − 1
   (equator nodes at 0° and 180°; positive when particle widens)

2. **ν_meas** = ε_l / ε_p (per-particle; ε_p already measured at 90°/270° nodes)

3. **Contact turning angle** = 180° − interior angle at each active contact node.
   *Active contact node* = node whose capsule surface is in contact with wall or
   neighbor (contact penalty force > threshold).
   Metric reported: mean turning angle over all active contact nodes.
   Baseline (undeformed N=32 circle): ~11.25°. Target: → 0° as contact flattens.

4. **No-side-wall mode**: make side walls optional (flag `side_walls=False`).
   Without side walls, particles are free to expand laterally. The particle-particle
   capsule contact still keeps the two disks vertically aligned.

**Gate:** New metrics computed and printed each frame in verbose mode.
           No-side-wall run stable (particles don't drift laterally).

---

### Waypoint 1B.2 — K_A/S characterization ✅  q_glass=0.1 (ΔA≈0.9%), q_rubber=10 (ΔA≈0.1%)

**Question:** How do ν_meas and ΔA respond to the ratio q = K_A/S?

**Protocol:**
- Two edge-case S values: S_low = 0.01, S_high = 1.0
- At each S: sweep q ∈ {0.1, 0.5, 1, 5, 10, 50, 200} by varying K_A = q × S
  (K_A is `rho_f` in the code)
- No side walls; squeeze to ε_p ≈ 10%; measure ν_meas and ΔA at ε_p = {5%, 8%, 10%}
- τ fixed at current default; C = 100

**Expected results:**
- ν_meas increases monotonically with q
- ΔA ≈ (1 − ν_meas) × ε_p confirms the 2D geometric relationship
- Saturation: ν_meas → ceiling (< 1.0) at finite q; record where it saturates

**Gate metrics:**
- ν_meas vs q is monotone at both S values
- ΔA and (1 − ν_meas) × ε_p agree within 20% (geometric consistency check)
- Identify q_rubber and q_glass anchor points hitting targets from VISION.md:
  - q_rubber: ν_meas ≥ 0.7, ΔA/ε_p < 30%
  - q_glass:  ν_meas ≤ 0.3, ΔA/ε_p > 50%

---

### Waypoint 1B.3 — τ characterization ✅  θ_wall: 9.67°(τ=0.05)→5.15°(τ=0.20) at q=10, S=0.1

**Question:** How does contact flatness (turning angle) respond to τ/S?

**Protocol:**
- Same two edge-case S values (S_low, S_high)
- At each S: sweep τ/S ∈ {0.001, 0.01, 0.05, 0.1, 0.5}
  (τ = (τ/S) × S; keep K_A at anchor values from 1B.2)
- Measure mean contact turning angle at ε_p = 5% and ε_p = 10%
- Also record penetration |cc_gap|/r_c

**Expected results:**
- Turning angle decreases monotonically as τ/S decreases
- Penetration decreases as turning angle decreases (flatter contact → lower peak force)

**Gate metrics:**
- Turning angle vs τ/S is monotone at both S values
- Identify τ/S value where turning angle < 3° (effectively flat) without excessive penetration
- Confirm flatness and penetration trade-off is acceptable at target τ/S

---

### Waypoint 1B.4 — Define the squishiness axis ✅  Four-metric plot saved; calibration table printed

**Deliverable:** A one-parameter family (q sweep from q_glass to q_rubber) with:
- τ/S = τ/S_flat (value found in 1B.3 that gives flat contact)
- K_A = q × S
- All other parameters (C, alpha0, dt) derived from S as in Phase 1

**Squishiness parameter:** q = K_A / S (dimensionless, range ~0.5 to ~100)

A summary figure: all four metrics vs q in a single 2×2 panel plot.

**Gate metrics:**
- ν_meas: monotone increasing from ≤ 0.3 to ≥ 0.7 over the q sweep
- ΔA: monotone decreasing (more area-conserving at high q)
- Contact flatness: approximately constant across q (controlled by τ/S, not q)
- Penetration: < 10% across all q

---

### Waypoint 1B.5 — Movies along squishiness axis ✅  4 GIFs in results/stage1b/movies/

Update movies to use the calibrated (q, τ/S) parameter set.
Produce 6–8 GIFs spanning q from glass-like to rubber-like.
Each frame shows all four metrics as overlaid text.

**Gate:** Visual inspection confirms:
- Low q (glass): particle keeps shape, area visibly shrinks, narrow contact patch
- High q (rubber): particle bulges laterally, area conserved, wide flat contact patch
- Smooth visual transition between the two

---

### Phase 1B Exit Gate

All of 1B.1–1B.5 pass. Four-metric summary plot produced. Movies reviewed by user.
Proceed to Phase 2 planning.

---

## Phase 1C: Dimensional Scaling Verification ✅ COMPLETE

**Objective:** Prove rigorously that S is a pure force scale and R0 is a pure length scale,
leaving τ and q as the only two free physics parameters.

### Background

After Phase 1B, the parameter set was: S, R0, τ, q, C/El_t, ζ. The question: can any of
these be eliminated by dimensional analysis, reducing the parameter space further?

### Result 1 — S is a pure force scale ✅

**Script:** `src/validation/scaling_verification.py`

Fix q=1.0, τ=0.05, C/El_t=0.208, α₀=2.0 at R0=1. Sweep S ∈ {0.01, 0.1, 1.0, 10.0}.
At each S: K_area = q×El_t, C = (C/El_t)×El_t, α = α₀/T_wave (ζ = const).

**Measured collapse quality (CV% across S values):**

| ε_p  | F/El_t CV% | ΔA CV% | ν CV% |
|------|-----------|--------|-------|
| 0.05 | 0.00      | 0.00   | 0.00  |
| 0.08 | 0.00      | 0.00   | 0.00  |
| 0.10 | 0.00      | 0.00   | 0.00  |

**Conclusion:** S is a pure force scale. Fix S=1 with no loss of generality.

### Result 2 — R0 is a pure length scale ✅ (after bug fix)

**Script:** `src/validation/R0_N_convergence.py`

Fix S=1, q=1.0, τ=0.05, C=1000 (const), K_area=4800 (const), α₀=2.0.
Sweep R0 ∈ {0.25, 0.5, 1.0, 2.0, 4.0} at N ∈ {32, 64, 128}.

**Bug discovered and fixed:**
The original code had `EI = S * K_fluid * R0**4`. For a thin shell of thickness t = τR0,
the correct bending stiffness is EI = E_l × t³/12 = S × K_fluid × **R0³**.
The extra R0 made bending grow faster than the membrane force, breaking R0 self-similarity.
Fix in `src/simulation/capsule_shell.py` line 129. At R0=1 both formulas give EI=S —
all prior calibrations are unchanged.

Also confirmed: C must be **constant** (independent of R0) for self-similar contact forces.
Contact force = k_c × gap × L0 = (C/R0) × R0 × R0 = C × R0 → F/force_unit = C/El_t = const.

**Measured collapse after fix (F/(El_t×R0) and ΔA across all R0 at each N):**

| N   | F_norm (all R0) | ΔA@10% (all R0) | Force CV% | ΔA CV% |
|-----|-----------------|-----------------|-----------|--------|
| 32  | 0.00398         | 0.365%          | 0.00      | 0.00   |
| 64  | 0.00432         | 0.378%          | 0.00      | 0.00   |
| 128 | 0.00504         | 0.399%          | 0.00      | 0.00   |

All five R0 values give identical results at each N to displayed precision.
The F_norm value shifts with N (normal contact-model convergence), but has zero R0 spread.

**Conclusion:** R0 is a pure length scale. Fix R0=1 with no loss of generality.

### Canonical parameter set (post Phase 1C)

| Parameter | Value | Role |
|-----------|-------|------|
| S         | 1.0 (fixed) | Force scale — eliminated |
| R0        | 1.0 (fixed) | Length scale — eliminated |
| C         | 3000×(1+q) (fixed formula) | Contact hardness — τ-INDEPENDENT |
| ζ = α₀/4π | 0.159 (fixed, α₀=2.0) | Damping ratio calibration constant |
| **τ**     | **free** | Shell thickness ratio; controls contact flatness |
| **q**     | **free** | K_area/El_t_base; controls squishiness (glass↔rubber) |

q is defined relative to El_t_BASE = 12S/τ² (not El_t_code = El_t_base × R0).
This definition is R0-independent: q_physical = K_area / El_t_base = const for K_area = const.

**Calibration anchors (τ=0.05):** q_glass=0.1 (ΔA@8%≈0.94%), q_rubber=10 (ΔA@8%≈0.10%)

### Results files

| File | Description |
|------|-------------|
| `results/scaling_verification/scaling_results.json` | S scaling raw data |
| `results/scaling_verification/scaling_collapse.png` | S scaling collapse plot |
| `results/R0_convergence/convergence_cache.json` | R0×N sweep raw data |
| `results/R0_convergence/convergence_CV_vs_N.png` | CV% vs N convergence plot |
| `results/R0_convergence/collapse_by_N.png` | Collapse curves at each N |

---

## Pre-Phase-2 C Calibration ✅ COMPLETE

**Objective:** Verify contact stiffness C is correctly formulated for all (τ, q), not just
the Phase 1 calibration point (τ=0.05, S=0.1).

**Finding:** The old formula C = C_RATIO × El_t_base scales C as 1/τ², because
El_t_base = 12S/τ². This gives C = 62 at τ=0.20 vs C = 1000 at τ=0.05.
But EI = S × K_fluid × R0³ is τ-independent (=1 at S=R0=1). The contact spring k_c = C/R0
must dominate EI/R0² = S, so C must scale with S, NOT El_t_base.

**Correct formula:** `C = 3000 × S × (1 + q)` (τ-independent, q-corrected)

The (1+q) factor accounts for the effective total stiffness (membrane + area spring) that the
contact must overcome. Calibrated by bisection: C_0 = 3000 gives cc_pen < 2% for all (τ, q).

**Verification (cc_pen = worst capsule-capsule penetration / L0 × 100%):**

| τ    | q    | C      | cc_pen% | status |
|------|------|--------|---------|--------|
| 0.05 | 0.1  | 3,300  | 1.02    | ✓      |
| 0.05 | 1.0  | 6,000  | 1.90    | ✓      |
| 0.05 | 10.0 | 33,000 | 0.74    | ✓      |
| 0.10 | 1.0  | 6,000  | 0.80    | ✓      |
| 0.20 | 1.0  | 6,000  | 0.48    | ✓      |

All cc_pen < 2%.  Runtime: ~19s/run (unchanged — edge wave, not contact spring, limits dt).

**Scripts:**
- `src/validation/C_calibration_grid.py` — initial τ×q grid (C0 vs C1 comparison)
- `src/validation/C_formula_verify.py` — formula comparison (C_old, C_fix, C_new, C_0=3000)

---

## Phase 1E: Rigid-Body Decomposition Integrator ✅ COMPLETE

**Objective:** Replace node-only integration with a physically correct split integrator
that damps only elastic deformation modes, preserves rigid-body momentum exactly, and
uses the correct solid-disk mass and moment of inertia.

**Design document:** `docs/rigid_body_integration_benchmark.md`

### Waypoints

| ID | Task | Status |
|----|------|--------|
| 1E.1 | Add rigid-body state to `CapsuleParticle` (x_cm, v_cm, θ, ω, u, u_dot, M_disk, I_disk) | ✅ |
| 1E.2 | Add `step_rb()` to `CapsuleSim` (split integrator) | ✅ |
| 1E.3 | Benchmark 0: free particle — momentum constant to machine precision | ✅ P0a=0, P0b=0, P0c<1e-12 |
| 1E.4 | Benchmark 1: single wall impact — COR measurement, energy closure | ✅ e(0)=0.923→e(8)=0.459, Eclos<1e-4 |
| 1E.5 | Benchmark 2: two-particle collision — momentum conservation, velocity exchange | ✅ P2a=9e-17, P2c=0.000 |
| 1E.6 | Re-run calibration squeeze with step_rb(); verify ν_meas and ΔA unchanged | ✅ |Δν|=0.001, |ΔΔA|=0.15% |

### Exit gate
All three benchmarks pass their numerical metrics. Calibration results unchanged.

---

## Phase 1F: Collision Physics Validation Suite ✅ COMPLETE

**Objective:** Validate that the EPD toy model, while not reproducing exact Hertz contact,
obeys the macroscopic physical laws required for a credible granular simulation:
momentum conservation, correct scattering geometry, and a well-characterised
dissipation axis (α_damp → COR calibration).  Also identify the quasi-static
velocity scale needed for Phase 2 compression protocols.

Three benchmarks, run in one script (`src/validation/rb_collision_benchmarks.py`):

---

### Waypoint 1F.1 — Off-centre collision: scattering geometry + momentum conservation ✅

**Protocol:**
- p2 stationary at origin; p1 approaches from below with velocity v0=0.05
- Impact parameter b ∈ {0, 0.2, 0.4, 0.6, 0.8} × contact_thr
- Three q values: 0.1 (glass), 1.0 (neutral), 10.0 (rubber)
- α_damp = 2.0 (suppress ringing so final state is clean)
- τ=0.2, N=32, n_post=1500

**Measured quantities (before contact vs after separation):**
- Linear momentum: Δp / p0  (should be < 1e-6 — exactly conserved by construction)
- Rigid-body angular momentum: ΔL / |L_ref|  (should be < 0.02)
- Deflection angle θ1 of p1 from original direction
- Final spin |ω1|, |ω2| (expect ~0 for frictionless contact)
- Compare θ1 to rigid-disk analytic: θ1_rigid = π/2 − arcsin(b/contact_thr)

**Gate metrics:**
- |Δp/p0| < 1e-6 for all b, all q  ✅ / ❌
- |ΔL/L_ref| < 0.02 for all b, all q  (rigid-body L conserved to 2%)  ✅ / ❌
- θ1 within 35% of rigid-disk prediction for q=0.1 (glass is closest to rigid)  ✅ / ❌
- Soft particles (q=10) deflect MORE than rigid-disk at same b (extended contact zone)  ✅ / ❌
- Final |ω| < 0.05 × v0/R0 for all cases (frictionless → near-zero spin)  ✅ / ❌

---

### Waypoint 1F.2 — v0 sweep: quasi-static velocity scale ✅

**Protocol:**
- Head-on collision (b=0), q=1.0, τ=0.2, N=32
- v0 ∈ logspace(−2.5, 0, 14)  [~0.003 to 1.0]
- Two α_damp values: 0.0, 2.0
- Measure COR = |v_rel_after| / |v_rel_before|

**Expected behaviour:**
- α_damp=0: COR flat at ~0.92 for all v0 in linear regime; drops at high v0 (nonlinear)
- α_damp=2: COR approximately flat at lower value, determined by viscous dissipation
- Crossover velocity v0_nl: where COR starts to drop sharply (nonlinear onset)

**Gate metrics:**
- α_damp=0: COR variation < 10% over v0 ∈ [0.003, 0.3]  (flat linear regime)  ✅ / ❌
- α_damp=2: COR variation < 10% over same range  ✅ / ❌
- COR(α=2) < COR(α=0) for all v0  ✅ / ❌
- Nonlinear drop identified; record v0_nl  ✅ / ❌

---

### Waypoint 1F.3 — α_damp sweep: COR calibration surface ✅

**Protocol:**
- Head-on, b=0, v0 = 0.05 (linear regime confirmed by 1F.2)
- α_damp ∈ {0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0}
- Three q values: 0.1, 1.0, 10.0
- τ=0.2, N=32, n_post=1500

**Expected behaviour:**
- COR(α=0) ≈ 0.92 for all q (radiative floor — model intrinsic)
- COR decreases monotonically with α_damp for all q
- ln(COR) linear in α_damp → fit COR = COR_0 × exp(−B × α_damp)
- B depends on q (softer particles have longer contact time → more dissipation per unit α)
  → glass (q=0.1) has smaller B; rubber (q=10) has larger B
- This gives a 2D calibration surface: (q, α_damp) → COR

**Gate metrics:**
- Monotone decreasing COR with α_damp for all three q values  ✅ / ❌
- Exponential fit R² > 0.90 for each q  ✅ / ❌
- B(q=10) > B(q=1) > B(q=0.1)  (contact duration ordering)  ✅ / ❌
- Record calibration table: for each q, what α_damp gives COR ∈ {0.9, 0.8, 0.7, 0.5}  ✅ / ❌

---

### Phase 1F Exit Gate

All three benchmarks pass gate metrics.  Produce:
1. Three-panel summary figure saved to `results/rb_benchmarks/collision_benchmarks.png`
2. JSON data file `results/rb_benchmarks/collision_benchmarks.json`
3. Printed calibration table: (q, α_damp) → COR
4. Append results to LOG.md; update HANDOFF.md

---

## Phase 1G: Extension Benchmarks ⬜ PENDING → ✅ COMPLETE (pending gate confirmation)

**Goal:** Two additional multi-particle tests that go beyond two-body collisions and probe
the physics the toy model needs to reproduce for jammed packings:
1. Multi-contact simultaneous collision (3-body)
2. Quasi-static rearrangement under oscillatory forcing (T1-event analogue)

**Script:** `src/validation/rb_extension_benchmarks.py`

### Waypoint 1G.1 — 3-body asymmetric collision ✅ COMPLETE

Three equal EPD particles: A at (−b, y_A) and B at (+b, y_B) fall downward while C at (x_off, 0)
moves upward.  ε offsets are computed so all three reach contact simultaneously.
x_off breaks the symmetry; b > contact_thr/2 prevents A–B overlap.

| Gate metric | Target | Measured |
|-------------|--------|---------|
| |Δp|/p_scale | < 1e-8 | ≤ 1.6e-15 (machine precision ✓) |
| |ΔL|/L_scale | < 1e-8 | ≤ 1.3e-14 (machine precision ✓) |
| Configs: q ∈ {0.1, 1, 10}, α ∈ {0, 2}, asymmetric v | all pass | PASS ✓ |

COR_eff: 0.800 (α=2, damped) → 0.924 (glass-like, α=0).

### Waypoint 1G.2 — Rearrangement squeeze test ✅ COMPLETE

Center particle pushed by oscillating wall v_wall(t) = v₀ + A·sin(ωt) through a tight gap
between two fixed arc obstacles (each equivalent to a second EPD particle of radius R0).
ω << ω_contact.  Tests sustained multi-contact, oscillatory stress loading, energy balance.

| Gate metric | Target | Measured |
|-------------|--------|---------|
| Particle passes through gap | True | True for all configs ✓ |
| W_wall_energy > 0 | True (wall does work) | ✓ W_wall = 2.5k–25k |
| E_elastic_max < 10 × W_wall | True (no blow-up) | ✓ ratio < 0.25 for all |
| PE near zero after exit | True | PE_exit < 14 vs peak 860–3562 ✓ |

Notable: rubber (q=10) requires 5–10× more wall work due to large r_c → wider contact area.

### Phase 1G Exit Gate

All benchmarks pass.  Produce:
1. `results/rb_benchmarks/extension_1G1.png` — 3-body conservation figure
2. `results/rb_benchmarks/extension_1G2.png` — squeeze test energy panels
3. `results/rb_benchmarks/extension_benchmarks.json`
4. Update LOG.md; update HANDOFF.md → advance to Phase 2

---

## Phase 2: Emulsion Droplet Model ✅ COMPLETE

**Branch:** `emulsion-dev`
**Tag at start:** `v1.0-elastic-shell-complete`

### Science goal
Add an emulsion droplet mode: a 2D fluid droplet modelled as a closed curve with
line tension γL, area incompressibility K_area, internal hydrostatic gravity, and the
same contact/damping/rigid-body machinery as the elastic shell.  Validate through a
sequence of three benchmarks culminating in a falling droplet that passes all gate
metrics.  Document the model, calibrations, and results in §9 of the methods paper.

### Dimensionless parameter space

| Symbol | Definition | Role |
|--------|-----------|------|
| κ | γ / (R₀ K_area) | pre-compression ratio; must be < 1 |
| Bo | ρ_d g R₀² / γ | Bond number: gravity vs surface tension |
| C̃ | C R₀ / γ | dimensionless contact stiffness |
| α̃ | α_damp × τ₀ | dimensionless damping (τ₀ = √(ρ_d R₀³/γ)) |
| k̃_reg | k_reg / γ | numerical regularization (not physical) |

Derived: A₀ = π R₀² / (1 − κ).  Working point: κ = 0.2, Bo ∈ [0.1, 2], C̃ = 500.

---

### Bucket A: Foundation — code refactor + new force terms

#### Waypoint 2A.1 — Rename rho_f → rho_d ⬜
- Search-replace `rho_f` → `rho_d` in all src/ files
- Update docstrings; no logic changes
- **Gate:** all existing validation scripts import and run without error

#### Waypoint 2A.2 — Fix hydrostatic gravity in CapsuleParticle ⬜
- Delete placeholder block (lines 282–288 of capsule_shell.py)
- Implement correctly: P_k = ρ_d g (y_cm − y_k), force = +P_k n̂_k ds_k (outward)
- Fx correction: subtract mean after accumulation
- **Gate (analytic):** single stationary particle, g > 0 → |Σ F_y + mg| < 1e-10,
  |Σ F_x| < 1e-10 (machine precision on closed polygon)

#### Waypoint 2A.3 — EmulsionParticle subclass ⬜
File: `src/simulation/emulsion_particle.py`

New parameters vs CapsuleParticle:
- `gamma` (line tension; replaces El_t and EI which are set to 0)
- `rho_d` (internal fluid density)
- `A0` derived from {gamma, K_area, R0}: A₀ = πR₀²/(1−κ)
- `m_node` = ρ_d π R₀² / N (fluid mass, not shell mass)

New force terms in overridden `accumulate_internal_forces()`:
1. **Line tension:** F_k^LT = −γ (t̂_{k,k+1} − t̂_{k−1,k}) — skipped if γ=0
2. **Tangential regularization:** F_reg,k = k_reg (ℓk⁺−ℓk⁻)/2 t̂_k,
   overdamped update inside shell step, NEVER added to RB force sum
3. Inherits area penalty and hydrostatic from parent

**Gate:** isolated EmulsionParticle at rest (g=0):
- Shape remains circular: max radial deviation < 0.1% R₀ after 10⁴ steps
- Area conserved: |A − πR₀²| / πR₀² < 0.5%
- Node spacing uniform: std(ℓk)/mean(ℓk) < 1%

---

### Bucket B: Stage 1 — Two-plate squeeze (no gravity)

File: `src/validation/emulsion_squeeze.py`
Movie: `results/emulsion/squeeze_N120.mp4`

**Setup:** Single EmulsionParticle (N=120, R₀=1) between top and bottom flat walls.
Squeeze to ε_p = 15% at quasi-static SR, then release.

#### Waypoint 2B.1 — Dimensional analysis and parameter recipe ⬜
Derive working-point parameter set from dimensionless groups.
- Fix γ = 1.0 (force scale), R₀ = 1.0 (length scale)
- κ = 0.2 → K_area = γ/(R₀ κ) = 5.0, A₀ = πR₀²/0.8
- C̃ = 500 → C = 500 γ/R₀ = 500
- Derive dt_max from capillary wave speed c = √(γ/(ρ_d R₀))
- k_reg: choose so spacing relaxes in < 0.1 contact time

**Gate:** documented parameter recipe table; dt_max formula verified numerically.

#### Waypoint 2B.2 — Squeeze mechanics validation ⬜
Metrics measured during squeeze cycle:
- Area conservation: |A − A_eq| / A_eq < 1% throughout
- Node spacing: std(ℓk)/mean(ℓk) < 5% throughout
- Force-displacement: smooth, symmetric, monotone loading; clean unloading
- No node crossing (min inter-node distance > 0)

**Gate:** all four metrics pass for κ ∈ {0.1, 0.2, 0.4}.

#### Waypoint 2B.3 — k_reg calibration ⬜
Sweep k_reg / γ ∈ {1, 10, 100, 1000}.
- Too small: std(ℓk) grows during contact → fail spacing gate
- Too large: introduces spurious oscillations → kinetic energy spike
- **Gate:** identify stable band; record recommended k̃_reg with <2% spacing error,
  <1% KE contamination relative to physical modes.

#### Waypoint 2B.4 — κ sweep and A₀ verification ⬜
Sweep κ ∈ {0.05, 0.1, 0.2, 0.3, 0.4, 0.45}.
- Verify A₀ formula: measured rest area matches πR₀²/(1−κ) within 0.5%
- Verify collapse limit: κ → 0.5 should show growing deformation
- **Gate:** A₀ formula confirmed; κ < 0.45 stable; κ ≥ 0.5 shows expected instability.

**Movie (2B):** squeeze cycle at working point (κ=0.2, N=120); show node positions,
area trace, force trace as overlay.

---

### Bucket C: Capillary-wave damping calibration (no gravity)

File: `src/validation/emulsion_damp_calib.py`
**Motivation:** Before testing multi-body contact or gravity, tune α_damp to the
natural oscillation timescale of the emulsion model.  The n=2 elliptical mode is
the relevant resonance (lowest non-trivial capillary wave):
  ω₂ = √(6 γ / (ρ_d R₀³))   →   √6 ≈ 2.449 at working point (γ=ρ_d=R₀=1)
  T₂ = 2π/ω₂ ≈ 2.57 τ₀
  Critical damping: α_crit = 2ω₂ ≈ 4.90

This mirrors the elastic-shell bending-mode calibration (Phase 1J).

#### Waypoint 2C.1 — Measure capillary wave frequency ✅
- Perturb single EmulsionParticle to elliptical (n=2) mode, release from rest
- Measure oscillation period T_measured from ellipticity signal
- **Gate:** |ω_measured/ω₂_analytic − 1| < 20% (discrete-vs-continuum gap due to k_reg)
- **Measured:** ω_measured = 2.5663 rad/τ₀, error = 4.8% ✓

#### Waypoint 2C.2 — α_damp sweep ✅
Sweep α_damp ∈ {0.5, 1.0, 2.0, 5.0, 10.0, 20.0} (all × 1/τ₀, since τ₀=1).
For each: measure Q = ω₂ / (2·decay_rate) from ringdown envelope.
- **Gate:** Q decreases monotonically with α_damp
- **Recommendation:** α_crit = 2ω₂_measured; record for use in Buckets D and E
- **Measured:** α_crit=5.13; underdamped at α=2, overdamped at α=5; transition confirmed ✓

---

### Bucket D: T1-event squeeze (no gravity)

File: `src/validation/emulsion_three_droplet.py`  (rewritten from symmetric squeeze)
Movie: `results/emulsion/t1_squeeze.mp4`

**Setup:** Three equal EmulsionParticles, g=0.  Left and right particle CMs are
pinned (perimeters free to deform).  Centre particle is pushed downward by a moving
top wall, squeezing it through the gap between the outer pair.
This mirrors the Phase 1G T1-event squeeze for the elastic shell.

Physical layout (R₀=1):
  Left  CM: (−d, 0)  d = 1.5 R₀  (surface gap between all pairs = 1.0 R₀, equilateral triangle)
  Right CM: (+d, 0)
  Centre CM: (0, √3·d) — equilateral-triangle start
  Top wall moves downward at SR=0.02 until centre y_cm < 0

Note: d changed from 1.2→1.5 R₀ to keep deformation manageable for the emulsion
model (no bending stiffness). Gap 0.4 R₀ caused polygon self-intersection.

#### Waypoint 2D.1 — T1 passage ✅
- Centre particle passes below the left/right centreline (y_cm_centre < 0)
- **Gate:** passage achieved without node crossing or area collapse
- **Measured:** T1 at t/τ₀=143.66 ✓

#### Waypoint 2D.2 — Area conservation ✅
- Max area loss during squeeze (pre-T1) < 25% for all three particles
  (K_area=5/κ=0.2 → ~20% compression at contact equilibrium is expected)
- **Measured:** Left=13.3%, Right=13.3%, Centre=19.0% ✓

#### Waypoint 2D.3 — Left-right symmetry ✅
- max |x_cm_centre| / R₀ < 3% throughout (geometry is exactly symmetric)
- **Measured:** 0.16% ✓

**Movie (2D):** T1 passage sequence; shape overlay showing deformation.

---

### Bucket E: Falling droplet — physical Bond number

File: `src/validation/emulsion_falling_droplet.py`  (rewritten)
Movie: `results/emulsion/falling_droplet_Bo{Bo}.mp4`

**Physical motivation (250 µm droplet, γ=10 mN/m, ρ_d=1000 kg/m³):**
  P_Laplace = γ/R₀ = 80 Pa;   ΔP_hydrostatic = ρ_d g 2R₀ = 2.45 Pa
  Bo = ρ_d g R₀²/γ ≈ 0.015   →   target range Bo ∈ {0.005, 0.01, 0.02, 0.05}

In simulation (γ=ρ_d=R₀=1): g = Bo directly.
α_damp from Bucket C calibration (near-critical, α ≈ 2ω₂ ≈ 5).

#### Waypoint 2E.1 — Squeeze sanity check, g=0 ✅
Single EmulsionParticle (CM pinned) squeezed to ε_p=10% by a wall, then released.
- (a) F_max > 0.5 γ/R₀ (contact active)
- (b) Max area error < 3% (K_area=5 allows ~2.35% at 10% squeeze)
- (c) F_final/F_max < 1% after release + 5τ₀ settling (reversible)
- **Measured:** F_max=14.0, A_err=2.35%, F_final=0.0 ✓

#### Waypoint 2E.2 — Trajectory accuracy, g>0 ✅
Free fall with g = 0.01 (Bo = 0.01), no floor.
- **Gate:** |y_sim(t) − y_cm0 − ½g t²| / (H·R₀) < 1%
- **Measured:** traj_err = 6.36e-05 ✓

#### Waypoint 2E.3 — Full falling droplet benchmark ✅
Drop from H=5R₀, Bo=0.01, α_damp=5.0 from Bucket C.

| Metric | Target | Measured |
|--------|--------|---------|
| Free-fall trajectory error | < 1% | 6.36e-05 ✓ |
| Node floor penetration | = 0 | 0.00e+00 ✓ |
| Settling: KE < 1e-4 × initial PE | within 40 τ₀ after impact | True ✓ |
| Sag shape | monotone with Bo ∈ {0.005, 0.01, 0.02} | sag={-0.003, 0.010, 0.018} ✓ |

**Movie (2E):** three Bo values; free fall → floor impact → settling with sag.

---

### Bucket F: Documentation

#### Waypoint 2F.1 — Paper section §9 ⬜
Add to `papers/summary_of_methods/main.tex`:
- §9.1: Emulsion model derivation (γL energy, A₀ formula, dimensionless groups)
- §9.2: Tangential regularization — rationale and pseudo-force formula
- §9.3: Internal hydrostatic gravity — derivation, Archimedes check
- §9.4: Parameter recipe (κ, C̃, α̃, k̃_reg, Bo → physical values)
- §9.5: Calibration results: k_reg sweep, κ sweep, ω₂ measurement, α_damp sweep
- §9.6: Benchmark results — Buckets B–E figures; falling-droplet metrics table;
         physical interpretation (250 µm droplet in oil)

#### Waypoint 2F.2 — Compile and verify ⬜
- xelatex twice; no errors; page count increases by ≥ 3 pages
- **Gate:** clean compile, all cross-references resolve

---

### Phase 2 Exit Gate

All of the following must pass before Phase 2 is marked complete:

1. ✅ Waypoint 2A.2 analytic Fy/Fx test (machine precision)
2. ✅ Waypoint 2B.2 squeeze mechanics (area, spacing, F(δ))
3. ✅ Waypoint 2B.3 k_reg calibration table
4. ✅ Waypoint 2C.1 capillary wave frequency — ω_meas=2.566, error=4.8% ✓
5. ✅ Waypoint 2C.2 α_damp sweep — α_crit=5.13, transition confirmed ✓
6. ✅ Waypoint 2D.1 T1 passage — t/τ₀=143.66 ✓
7. ✅ Waypoint 2D.3 T1 symmetry — x_cm_centre max=0.16% R₀ ✓
8. ✅ Waypoint 2E.1 squeeze sanity — F_max=14.0, area=2.35%, reversible ✓
9. ✅ Waypoint 2E.3 falling droplet — all 4 metrics pass at Bo=0.01 ✓
10. ✅ Movies: t1_squeeze.mp4 (482 frames); falling_droplet_Bo{0.005,0.010,0.020}.mp4 (234 frames each)
11. ✅ §9 in paper compiled cleanly — 13 pages, no errors

---

## Phase 3: DtN Training Data Generation ✅ COMPLETE

**Objective:** Generate FEM training data for the DtN operator across the squishiness axis.

### Waypoints

| ID | Task | Status |
|----|------|--------|
| 3.1 | Two-disk sweep at N=240 (10 materials × 13 δ/R values) | ✅ 200 configs, 2400 aug samples |
| 3.2 | `two_disk_contact.py` robust solver with continuation | ✅ |
| 3.3 | Three-disk equilateral-triangle sweep (multi-contact loading) | ✅ 240 configs, 2880 aug samples |
| 3.4 | Assemble combined dataset | ✅ `dtn_combined.h5` — 5280 samples (5256 non-deg.) |

**Dataset summary:**
- `data/processed/dtn_twodisk_aug.h5` — 2400 samples (2376 non-degenerate)
- `data/processed/dtn_threedisk_aug.h5` — 2880 samples
- `data/processed/dtn_combined.h5` — 5280 merged samples

**Physics normalization:** F_input = F/max|F|, u_output = u·E/max|F| (collapses E dependence).

---

## Phase 4A: DtN Architecture Search ✅ COMPLETE

**Objective:** Find the best NN architecture for the DtN operator on the current dataset.

### Key Finding: 5% Linear Floor

After physics normalization, the DtN map is well-approximated as linear.
Best possible linear model achieves ~3-9% residual per material class.
All trained architectures plateau near 5.4–5.9% — near the linear model limit.

**Root cause:** 5% floor = variability due to different contact widths (δ/R) across samples.
Not a model capacity issue — more data with same geometries won't help.

### Architecture Results (300 epochs, 5256-sample combined dataset)

| Architecture | Params | Val rel L2 | Notes |
|-------------|--------|------------|-------|
| FourierDtN | 68K | 5.42% | Fourier-domain per-mode matrices |
| FourierDtNv2 | 160K | **5.39%** | Deeper conditioning MLP, mode-scale param |
| FourierDtNTrunc k≤40 | 26K | 5.91% | Only low modes |
| FourierResidual | 195K | 5.55% | FourierDtNv2 + CircCNN residual |
| DtNMLP baseline | 5M | ~15% | Flat MLP, no rotation equivariance |
| CircCNN | 980K | TBD | |

**Best model:** `results/dtn_fourier_v2_combined.pt` — FourierDtNv2, 5.39%

**Fix:** Dataset scale-up with diverse contact geometries (Phase 4B).
Architecture is near-optimal; more architectures won't help until floor is broken.

---

## Phase 4B: Dataset Scale-Up (diverse contact geometries) ⬜ IN PROGRESS

**Objective:** Add 4-disk and 5-disk contact configurations to break through the 5% floor.
Each new geometry adds different multi-contact loading patterns not represented in 2/3-disk data.

### Waypoints

| ID | Task | Status |
|----|------|--------|
| 4B.1 | Four-disk square contact simulation + data generation | ⬜ Running (PID 1054631) |
| 4B.2 | Assemble expanded combined dataset (2+3+4 disk) | ⬜ After 4B.1 |
| 4B.3 | Retrain FourierDtNv2 on expanded dataset | ⬜ After 4B.2 |
| 4B.4 | Gate: val rel L2 < 4% on expanded dataset | ⬜ |
| 4B.5 | (If 4B.4 fails) Five-disk cross simulation + retrain | ⬜ Optional |

**4-disk geometry:** Square arrangement with orthogonal contacts (0° and 90°/270° per disk).
Complements existing data (3-disk had 60°/120° contacts).

**Expected dataset size after 4B.2:**
- 4-disk sweep: 10 materials × 8 δ values × 4 disks × 12 rotations = 3840 aug samples
- Total combined: 5280 + 3840 = **9120 samples**

**Gate (4B.4):** val rel L2 < 4% on combined 9K-sample dataset.

---

---

## Phase 3: Accelerated TF + C++ Simulation Engine

**Goal:** Replace the pure-Python simulation loop with a TF-based force kernel pipeline
backed by a C++ candidacy manager. End state: same initial conditions run through the
frozen `emulsion-dev` branch or the new branch produce identical results to machine
precision.

### Architecture

```
C++ CandidacyManager  ─── CapCandidates (K×E int32) ──→  TF step graph
     ↑                                                          │
     └──────── x_cm (P,2),  theta (P,) ────────────────────────┘
```

- **C++** owns neighbor topology (cell list + contact-normal registration). Writes
  `CapCandidates` in place. Never computes distances at the edge level.
- **TF** owns all arithmetic (internal forces, contact forces, RB decomposition,
  integration). One `@tf.function(jit_compile=True)` per step. Fixed graph, no
  recompilation.
- **Python** orchestrates: calls C++ manager when update trigger fires, calls TF step,
  pulls metrics when needed.

### Fixed Data Structures (never change shape after init)

```
P   — number of particles (3–4 initially, up to ~10 in testing)
N   — nodes per particle (uniform)
E   — max candidates per edge (skin-buffered, set at sim init)
K   = P * N  — total edges; edge global index k = p*N + e ∈ [1..K]; 0 = ghost

CapCandidates : (K, E)       int32   — each row: E candidate edge indices (0=ghost)
x_all         : (P, N, 2)    DTYPE   — world-frame node positions
x_cm          : (P, 2)       DTYPE   — rigid-body centers
v_cm          : (P, 2)       DTYPE
theta         : (P,)         DTYPE   — orientation angles
omega         : (P,)         DTYPE
u             : (P, N, 2)    DTYPE   — elastic displacement (body frame)
u_dot         : (P, N, 2)    DTYPE
params        : (P,) each    DTYPE   — R0,N,El_t,EI,K_area,k_c,m_node,M_disk,I_disk,L0,A0

DTYPE = tf.float32 by default; switchable to tf.float64 via module-level flag.
C++ buffers use float/double matching DTYPE — no silent truncation.

Ghost element 0: x_ghost >> box size; always gap > contact_r → zero force.
```

### Concepts locked

- If two particles aren't neighbors → their edges can't overlap → filter at particle level
- If neighbors → only contact-facing edges are candidates → Level 2/3 filter
- Refresh infrequently; skin buffer absorbs registration drift
- Update trigger: `max(|Δx_cm| + R0·|Δtheta_A − Δtheta_B|) > skin/2` per pair
- Double-counting contacts is intentional (each row sums independently → pure reduce)
- No scatter anywhere in TF graph; all reductions are `reduce_sum` over a batch axis
- Debug/diagnostic tooling is first-class from day one

---

### Phase 3.1 — TF Intra-Particle Kernels + Integration ⬜

**Goal:** All per-particle physics in TF. No inter-particle forces yet.
Validated against frozen branch to machine precision (single-step state diff).

**Waypoints:**

#### 3.1.1 — TF install + dtype flag ⬜
- Install TensorFlow in venv, verify `@tf.function(jit_compile=True)` works on CPU
- Add `DTYPE` flag (default `tf.float32`; switchable to `tf.float64`) to new
  `src/simulation/tf_sim.py`
- Smoke test: create `(P,N,2)` tensors, run one matmul, confirm XLA compiles

#### 3.1.2 — Batched internal force kernel ⬜
- Port `accumulate_internal_forces()` → `internal_forces_tf(x_all, params)` returning
  `(P, N, 2)`. All ops: roll + arithmetic + reduce_sum. No scatter.
  - Fluid pressure: area via shoelace, outward force via node normals
  - Edge elasticity: strain × tangent, distributed via roll
  - Bending: turning angle hinge, 3-node stencil via roll
  - Hydrostatic gravity (optional, gated by `g != 0`)
- Gate: `single_step_diff` vs frozen branch = 0 ULP on internal forces for a
  saved 3-particle state

#### 3.1.3 — RB decomposition + integration kernel ⬜
- Port `step_rb()` → `step_rb_tf(state, f_ext, dt, alpha_damp)` returning new state
  - F = Σf_contact, T = Σr×f_contact (contact forces only; elastic terms cancel)
  - Undamped Verlet: v_cm, x_cm, omega, theta
  - Deviatoric force f_dev = f_phys − F/N (zero-mean)
  - Damped Verlet in body frame: u_dot, u
  - Perimeter reconstruct: x_i = x_cm + R(θ)(X_ref + u_i)
- Gate: single-step diff vs frozen branch = machine precision on all state tensors

#### 3.1.4 — Metrics kernel + hard wall / periodic BC ⬜
- Metrics inside graph (returned each step at negligible cost):
  `ΔA/A0 (P,)`, `KE_rb (P,)`, `KE_el (P,)`, `PE (P,)`, `pen_max (P,)`
- Hard wall BC: gap/normal as batched TF ops (LineSegment, Arc primitives)
- Periodic BC: minimum-image convention in gap computation; Lees-Edwards stub
  (interface defined, shear physics later)
- Gate: 3-particle wall-bounce, 100 steps, energy conserved to < 1e-10 relative drift

**Milestone:** 3-particle free-relaxation from rcpgenerator initial condition, run 500
steps through frozen branch AND TF branch, final state diff < 1e-12 (float64 mode).

---

### Phase 3.2 — Python Candidacy Manager + Inter-Capsule Forces ⬜

**Goal:** Working candidacy manager and contact forces in pure Python using
`CapCandidates`. Slow but correct. Debug tooling built in. Validated to machine
precision against frozen branch.

**Waypoints:**

#### 3.2.1 — CandidacyManager (Python) ⬜
- `src/simulation/candidacy_manager.py`: `CandidacyManager` class
  - Level 1: center-center cell list (NumPy), O(N_particles)
  - Level 2: contact-normal registration — facing edge index from `x_cm` + `theta`
    via angle arithmetic, O(1) per pair
  - Level 3: index-sliding fill of `CapCandidates` rows with `j ± dj` (ghost padding)
    `dj` scales with `N_B/N_A` for size ratio support
  - Update trigger: `max(|Δx_cm| + R0·|Δtheta_A − Δtheta_B|) > skin/2`
  - `skin` and `E` set at init; warn if any active contact near E capacity

#### 3.2.2 — Debug tooling ⬜
Built into `CandidacyManager` (always compiled in, cheap when inactive):
- `check_missing_contacts(x_all, r_contact)`: brute-force scan for pairs within
  contact distance absent from `CapCandidates` — emits mismatch report
- `candidate_utilization()`: fraction of E slots active this step
- `penetration_log`: list of (pair, gap) for any gap < 0 not in candidates
- `contact_count_histogram()`: contacts per particle
- Input/output logger: when `log=True`, saves `(x_cm, theta) → CapCandidates` for
  every update call to a replay corpus (used for 3.4 C++ validation)

#### 3.2.3 — Inter-capsule force kernel (Python, CapCandidates-driven) ⬜
- `inter_capsule_forces_py(x_all, CapCandidates, params)` → `(P, N, 2)`
  - Gather candidate edge positions from `CapCandidates` indices
  - 2-pt Gauss quadrature gap + force (same algorithm as frozen branch)
  - reduce_sum over E axis → forces on source edges → distribute to nodes via
    roll (no scatter)
  - Double-counts contacts; both sides correct by construction
- Gate: single-step force diff vs frozen branch `_capsule_capsule_forces` = machine
  precision for a 3-particle state with known contacts

#### 3.2.4 — Full Python integration loop ⬜
- Wire `CandidacyManager` + inter-capsule forces + 3.1 TF kernels into a single
  `step()` function in `src/simulation/tf_sim.py`
- Run 3-particle rcpgenerator initial condition through frozen branch and new branch
  for 500 steps; diff all state tensors at every step
- Gate: all diffs < machine precision (float64); zero missing contacts reported

**Milestone:** 4-particle dense packing (rcpgenerator), 200 steps, hard walls,
full state diff vs frozen branch < 1e-12 at every step. Log corpus saved.

---

### Phase 3.3 — Inter-Capsule Forces on TF, Full Graph Assembled ⬜

**Goal:** Complete single-step TF graph — internal + contact + integration — as one
`@tf.function`. Python calls one function per step.

**Waypoints:**

#### 3.3.1 — TF inter-capsule force kernel ⬜
- `inter_capsule_forces_tf(x_all, CapCandidates, params)` → `(P, N, 2)`
  - `tf.gather` on flattened `x_all` using `CapCandidates` indices
    → `(K, E, 2)` candidate edge endpoint positions
  - Gauss quadrature gap+force: pure tensor arithmetic, `reduce_sum` over E axis
  - Reshape `(K, 2)` → `(P, N, 2)` and add to internal forces
  - `CapCandidates` enters as `tf.Variable` (written by C++ in 3.4; Python writes now)
- Gate: single-step diff vs Phase 3.2 Python kernel = machine precision

#### 3.3.2 — Fused step function ⬜
- Combine 3.1 + 3.3.1 into one `@tf.function(jit_compile=True)`:
  `tf_step(state, CapCandidates, dt, alpha_damp, g)` → `(new_state, metrics)`
- `CapCandidates` is a `tf.Variable` updated outside the graph by the candidacy
  manager; the graph reads it as a captured variable (no recompilation on update)
- Gate: graph compiles once, no retracing on subsequent steps

#### 3.3.3 — End-to-end validation ⬜
- 4-particle rcpgenerator packing, 500 steps, hard walls
- Frozen branch vs TF branch: all state diffs < machine precision every step
- Log full candidacy manager input/output corpus for 3.4

**Milestone:** single `tf_step()` call per timestep; 4-particle 500-step run matches
frozen branch to machine precision; replay corpus saved to `results/phase33_corpus/`.

---

### Phase 3.4 — C++ Candidacy Manager + pybind11 ⬜

**Goal:** C++ replaces Python candidacy manager. `pip install .` works on Colab.
Validated against Phase 3.3 logged corpus without re-running simulations.

**Waypoints:**

#### 3.4.1 — Repo structure for pip install ⬜
```
pyproject.toml          ← build-backend: scikit-build-core or setuptools+cmake
                           build-requires: pybind11>=2.11, cmake>=3.18, ninja
CMakeLists.txt          ← top-level; finds Python + pybind11; builds _neighbor_manager
src/cpp/
  neighbor_manager.cpp  ← CandidacyManager: cell list, registration, CapCandidates fill
  neighbor_manager.h
  bindings.cpp          ← pybind11: expose CandidacyManager.step(x_cm, theta, buf)
```
`pip install .` on Colab: cmake auto-invoked by build backend, no manual steps.
DTYPE-aware: template on `float`/`double`; Python binding selects at construction.

#### 3.4.2 — C++ CandidacyManager ⬜
- Same three-level logic as Phase 3.2 Python version:
  - Cell list (cell size = 2·R0_max + skin)
  - Contact-normal registration (angle arithmetic, no distance at edge level)
  - Index-sliding fill of `CapCandidates` buffer (written in place)
- Update trigger identical to Python version
- Debug mode flag: enables brute-force contact scan, emits JSON mismatch reports
- `step(x_cm_np, theta_np, cap_candidates_np)` → writes `cap_candidates_np` in place
  (NumPy array shared with TF variable; zero-copy)

#### 3.4.3 — Corpus replay validation ⬜
- Load Phase 3.3 corpus: sequence of `(x_cm, theta)` inputs and expected
  `CapCandidates` outputs
- Feed each input through C++ manager; diff output vs logged Python output
- Gate: zero mismatches across full corpus
- No full simulation needed for this validation step

#### 3.4.4 — End-to-end integration test ⬜
- Wire C++ manager into full simulation loop (replaces Python manager)
- Run same 4-particle 500-step test as Phase 3.3
- Gate: all state diffs vs frozen branch < machine precision (float64)

**Milestone (Phase 3 complete):** same initial conditions → frozen branch or new
TF+C++ branch → identical final state to machine precision (float64). `pip install .`
succeeds on a clean Colab instance (verified in CI or manually).

---

### Phase 3.5 — Fast Handshaking + Scale-Up ✅ COMPLETE

- Primitive forces (LineSegment, Arc, Polygon/Box, Hopper) — all validated < 5e-13
- Polydisperse inter-capsule forces via (K,) flat arrays
- tf.while_loop runner with C++ candidacy via tf.py_function
- P=10 N=60 benchmark: 1.7× CPU speedup; inter_capsule 90% of cost
- render_utils.py: universal capsule outer-contour renderer + save_gif
- E reduced from 64 to 24 for P=10 benchmark (2.2× extra speedup, no missed contacts)

---

### Phase 3.6 — Flow Primitives: Driven Particles, Periodic BC, Shear Infrastructure ⬜

**Goal:** Build all DEM primitives needed for shear, flow, and jamming experiments.
Not the experiments themselves — the roads they run on.

**Architecture additions:**

```
params['driven_mask']  (P,)     — 0=force-integrated, 1=velocity-prescribed
params['shape_frozen'] (P,)     — 0=elastic shape, 1=rigid (nodes locked to RB frame)
params['traj']         (P, 18)  — DC+AC trajectory: vx(4), vy(4), ω_spin(4), ω_orbit(4), r_ref(2)
params['box']          (2,)     — [Lx, Ly] for periodic BC (default 1e9 = non-periodic)
```

Integration mask (zero extra GPU cost beyond two element-wise multiplies):
```
v_new = v_force*(1-mask) + v_prescribed(t)*mask
```

---

#### Waypoint 3.6.1 — Periodic BC in CandidacyManager ✅

- Add `periodic=False`, `Lx`, `Ly` to `__init__`
- `_level1_pairs`: use minimum-image distance `|Δx_mi| = Δx - L*round(Δx/L)` per component
- `_fill_pair`: wrap `n_AB` direction vector with minimum image before angle arithmetic
- `needs_update`: use minimum-image displacement for trigger condition
- Default `periodic=False` — fully backward compatible
- **Gate:** `check_missing_contacts` finds zero missing pairs in a 25-particle periodic packing

#### Waypoint 3.6.2 — Periodic BC in inter_capsule_forces_tf ✅

- Add `box=None` argument — `(2,)` DTYPE tensor `[Lx, Ly]`; `None` = non-periodic
- Apply minimum image to `diff0 = xq - b0` before computing `t_B` and gap:
  `diff0 -= Lx * tf.math.round(diff0/Lx)` per component
- Pass `box=params['box']` from `step_full_tf`
- **Gate:** two-particle force test across periodic boundary matches analytic Hertz gap to < 1e-12

#### Waypoint 3.6.3 — Driven particle mask + DC+AC trajectory in step_rb_tf ✅

- Add `t` argument (scalar DTYPE tensor) to `step_rb_tf`
- Read `driven_mask (P,)`, `shape_frozen (P,)`, `traj (P,18)` from `params`
- Evaluate prescribed velocity: `v_x(t) = dc + ac*cos(freq*t+phase)` for each DOF
- Orbital contribution: `v_cm += ω_orbit * ⊥(x_cm - r_ref)`
- Integration: `v_new = v_force*(1-mask) + v_prescribed*mask`
- Shape freeze: `u_new = u_elastic*(1-frozen) + 0*frozen` — locks nodes to RB frame
- Add defaults in `make_state`: zeros for all new params
- Add helper: `make_traj(v_dc, v_ac, freq, phase, omega_spin, omega_orbit, r_ref)` → (18,) array
- Add helper: `set_driven(params, indices, traj_rows, frozen=False)` — modifies params in place
- **Gate:** single driven particle follows `v_cm=(1,0.5)` exactly; orbital particle circles reference

#### Waypoint 3.6.4 — Rotating primitives in make_prim_data + primitive_forces_tf ✅

- Add `omega=0, r_ref=(0,0)` per primitive in `make_prim_data` call signature
  - New prim_data keys: `seg_omega (Ls,)`, `seg_r_ref (Ls,2)`, `arc_omega (La,)`, `arc_r_ref (La,2)`
- In `primitive_forces_tf`: compute current endpoint positions at time `t`:
  `p0(t) = r_ref + R(ω*t)*(p0_ref - r_ref) + v_trans*t`; normal also rotates
- Default omega=0 everywhere — backward compatible
- **Gate:** rotating LineSegment sweeps correct arc; phase35_validate.py still passes

#### Waypoint 3.6.5 — rcpgenerator → EPD bridge ✅

File: `src/simulation/rcp_utils.py`

- `rcp_seed(N_particles, N_nodes, polydispersity, seed, walls)` → `(particles, Lx, Ly, phi_J)`
  - Calls `rcpgenerator.Packing` to get positions at φ_J
  - Scales to EPD units: `R_eff = R0 + r_c = 1 + 2π/N`; `scale = 2*R_eff / d_rcp`
  - Returns `CapsuleParticle` list + box dimensions
- `scale_packing(particles, Lx, Ly, phi_target, phi_J)` → `(particles, Lx_new, Ly_new)`
  - Affine rescale: `x_cm *= sqrt(phi_J/phi_target)`, `Lx *= sqrt(phi_J/phi_target)`
  - Returns particles with rescaled CMs and new box dimensions

#### Waypoint 3.6.6 — Phase 3.6 Benchmark: RCP + Swell + Shear ✅

File: `src/validation/phase36_shear.py`

**Protocol:**
1. Seed 25 particles via `rcp_seed(N_particles=25, N_nodes=60, walls=[0,0])` → positions at φ_J
2. Scale box down to φ_start = φ_J - 0.04 (particles just not touching, ~3% linear gap)
3. Run EPD with periodic BC, slow compression (top+bottom walls squeeze) until φ = φ_J + 0.04
4. Identify top row (y_cm > y_threshold) → `driven_mask=1, v_cm=(+V, 0), shape_frozen=0`
5. Identify bottom row (y_cm < y_threshold) → `driven_mask=1, v_cm=(-V, 0), shape_frozen=0`
6. Run shear for 2000 steps, render movie

**Kinematic tests (quick, run before full benchmark):**
- K1: Single particle driven with `v_cm=(1, 0.5)` → verify trajectory to < dt tolerance
- K2: Single particle with `omega_orbit`, circular motion around ref point → position on circle
- K3: Rotating LineSegment (omega=1) → endpoint traces correct arc at time t

**Gate metrics:**
- K1–K3 kinematic tests pass
- Periodic BC: particle crossing x boundary shows no force discontinuity
- Shear movie shows top layer moving right, bottom left, bulk deforming smoothly
- No numerical explosion over 2000 steps (max |v_cm| < 10)

---

### Phase 3.6 Exit Gate

All of 3.6.1–3.6.6 pass.  Documents updated.  Ask user for next direction.

---

## Simulation Environment

| Component | Status |
|-----------|--------|
| Python venv (.venv/) | ✅ Active |
| scikit-fem, h5py, torch, matplotlib, imageio | ✅ Installed |
| rcpgenerator | ✅ Built from source (deps/RCPGenerator) |
| polyfempy | ❌ Build failed (GCC 13) — not needed for capsule DEM |
| GPU | ❌ CPU only — N=32, single runs ~70s |

---

*Last updated: 2026-04-18 — Phase 1 marked complete; Phase 1B (squishiness axis
calibration) defined with four-metric framework and q = K_A/S squishiness knob*

---

## Phase 4: High-Level API (`src/epd/`) ⬜ PENDING

### Governing rules for Phase 4

1. **Engine is frozen.** `capsule_shell.py`, `tf_sim.py`, `candidacy_manager.py`,
   `contact_primitives.py` are read-only except for explicitly listed bug patches.
   Phase 4 is pure orchestration above the engine.

2. **Bug patches are pre-requisites.** All patches listed in Phase 4.0 must be
   applied and verified before any 4.1+ work begins.

3. **Backward compatibility.** Every unspecified material/force parameter defaults
   to zero/None. New force kernels are no-ops when their parameter is absent.
   Existing `ParticleSpec` call sites never break.

4. **No Python in the step loop.** Everything time-critical (frozen shape, motion
   application, force computation) happens inside TF via masks and tensors.
   Python only runs outside the loop (setup, candidacy update, bookkeeping).

5. **Each phase has one runnable test.** Run `python src/epd/tests/test_4X.py`.
   It produces an image or movie. Developer checks output visually.

6. **Autonomous execution.** Once Phase 4.0 is confirmed, phases 4.1–4.6 execute
   autonomously in order. Each phase ends with a passing test before proceeding.

---

## Phase 4.0: Design Document & Pre-Requisite Patches ✅ COMPLETE

### Forward graph

```
User layer                         Translation                    Engine
──────────────────────────────────────────────────────────────────────────
ParticleSpec                       CapsuleParticle list           capsule_shell.py
  .type = 'elastic'|'emulsion'     → force kernel selector        tf_sim.py
  .nu / .q, .tau_b                 → q, TAU, K_area, C, alpha     (derived)
  .N_nodes, .R0_mean, .poly_dist   → per-particle R0_arr          (sampled)
  .extra_forces = {drag: γ, ...}   → extra force tensors          (future)
  .motion = MotionSpec             → driven_mask + traj_fn(t)     step_rb_tf
  .frozen_shape = bool             → frozen_mask (P,) bool        step_full_tf

Object (Wall/Arc/Composite)        prim_data dict                 primitive_forces_tf
  .exclusion = 'interior'|'exterior'  → seeder domain queries     initializer.py
  .motion = MotionSpec             → vel_fn(t), omega_fn(t), r_ref

System(Lx, Ly, bc_x, bc_y)
  state = {                        ← complete resume state
    x_all    (P, N, 2) float64     node positions
    x_cm     (P, 2)    float64     CM positions
    v_cm     (P, 2)    float64     CM velocities
    omega    (P,)      float64     angular velocities
    theta    (P,)      float64     rotation angles
    u        (P, N, 2) float64     elastic displacements
    u_dot    (P, N, 2) float64     elastic velocities
    X_ref    (P, N, 2) float64     reference node positions
    frozen_mask (P,)   bool        shape-frozen particles
    t        scalar    float64     absolute simulation time
    step     scalar    int64       total steps taken
    Lx, Ly   scalar    float64     current box dimensions
    params   dict                  TF params (periodic box, dt, alpha)
    CapCandidates (P*N, E) int32   candidacy table (rebuilt on load)
  }
```

### State completeness rule

A checkpoint is complete if and only if `system.load(path)` followed by
`system.step(N)` produces bit-identical output to a reference run that never
stopped. This is verified in Phase 4.5.

### Required patches to existing engine files

#### Patch A — `candidacy_manager.py`: per-pair R0 threshold

**Problem:** `_level1_pairs` and `_fill_pair` use `2 * self.R0 + self.skin` as a
global threshold. For polydisperse or mixed-type systems, the correct threshold
per pair (i, j) is `R0_arr[i] + R0_arr[j] + skin`.

**Fix:** Pass `R0_arr` (shape `(P,)`) to `__init__`. Store it. In `_level1_pairs`
and `_fill_pair`, use `R0_arr[pA] + R0_arr[pB] + skin` for each pair.
`self.R0` becomes `np.mean(R0_arr)` for display only.

**Test:** Two particles with R0=0.8 and R0=1.2. Old code uses threshold 2×1.0+skin;
new code uses 0.8+1.2+skin=2.0+skin (same in this case). Test with R0=0.5 and
R0=1.5: old threshold=2.0+skin, new=2.0+skin (still same). Test with R0=1.0 and
R0=2.0: old=2.0+skin, new=3.0+skin — new correctly extends range.

#### Patch B — `tf_sim.py`: frozen_shape mask in step_full_tf

**Problem:** No mechanism to clamp elastic DOFs to zero for rigid obstacle particles.

**Fix:** Add `frozen_mask` parameter (shape `(P,)` bool tensor, default all False)
to `step_full_tf`. After elastic displacement update:
```python
u     = tf.where(frozen_mask[:, None, None], tf.zeros_like(u), u)
u_dot = tf.where(frozen_mask[:, None, None], tf.zeros_like(u_dot), u_dot)
```
This is inside TF — zero Python overhead.

**Test:** Single particle with `frozen_mask=True`, initial u≠0. After one step,
verify u=0, u_dot=0. Second particle with `frozen_mask=False` is unaffected.

#### Patch C — `tf_sim.py`: absolute time `t` threaded through motion

**Problem:** `step_full_tf` currently accepts `t` as a parameter but only passes
it to `primitive_forces_tf`. Driven particle trajectory interpolation also needs `t`.

**Fix:** Ensure `t` (absolute sim time, float64 scalar tensor) is passed to all
motion-dependent functions. `step_full_tf` signature already has `t`; verify it
reaches both `primitive_forces_tf` and the driven-particle interpolation path.
No change needed if already correct — this is a verification step.

#### Patch D — `candidacy_manager.py`: mixed-type skin per particle type

**Problem:** Emulsion particles have a different effective contact radius than
elastic particles (surface tension replaces edge spring; r_c still defined but
physically different). Skin should be per-particle-type.

**Fix:** Accept optional `skin_arr` (shape `(P,)`) alongside scalar `skin`.
If provided, threshold for pair (i,j) = `R0_arr[i] + R0_arr[j] + 0.5*(skin_arr[i] + skin_arr[j])`.
If not provided, falls back to scalar skin (backward compatible).

---

## Phase 4.1: `objects.py` — Object Hierarchy ⬜ PENDING

### File: `src/epd/objects.py`

#### Class hierarchy

```
SimulationObject          base class, dict-like properties
  Wall                    single line segment primitive
  Arc                     single arc primitive
  CompositeObject         list of primitives + origin/rotation transform
    Box                   4 Wall primitives, inward normals
    Channel               2 Wall primitives (top+bottom)
    CouetteCell           2 Arc primitives (inner excluded, outer container)
    CircleObstacle        1 Arc primitive (interior excluded)
    RegularPolygon        N Wall primitives
    CustomObject          user-assembled composite
```

#### Key methods

`SimulationObject`:
- `__init__(kind, **props)` — stores kind + properties dict
- `set_exclusion(mode)` — `'interior'` or `'exterior'` or `None`
- `set_motion(motion_spec)` — attaches a `MotionSpec` (Phase 4.3)
- `resolved(t)` → `list[dict]` — returns flat list of primitive dicts at time t,
  with motion applied. Each dict has keys: `kind`, `p0/p1/normal` or `center/radius`,
  `vel`, `omega`, `r_ref`, `exclusion`.

`CompositeObject` adds:
- `add_primitive(obj)` — append to primitive list
- `set_origin(x, y)` — translation offset applied to all children in `resolved()`
- `set_rotation(theta)` — static orientation (not motion)
- `region_polygon(t)` → `{'vertices': [...], 'exclusion': str}` or `None`
  Traces the boundary by concatenating primitive path vertices in order.
  Used by seeder for inside/outside point tests.
  Returns `None` if object has no closed boundary (e.g. a single wall).

#### Pre-built objects

`Box(width, height, x0, y0, theta, exclusion='exterior')`:
- 4 Wall primitives at correct local positions with inward normals
- `region_polygon()` returns the 4-corner polygon

`Channel(width, height, x0, y0, theta, exclusion='exterior')`:
- Top and bottom walls only (periodic in x assumed by user)
- `region_polygon()` returns the bounding rectangle

`CouetteCell(inner_radius, outer_radius, x0, y0, exclusion='exterior')`:
- Outer arc (particles inside, exclusion='exterior') +
  inner arc (particles outside inner cylinder, exclusion='interior')
- `region_polygon()` returns annulus approximated as polygon (64 pts)

`CircleObstacle(radius, x0, y0, exclusion='interior')`:
- Single arc, particles excluded from interior

`RegularPolygon(sides, radius, x0, y0, theta, exclusion='interior')`:
- N Wall primitives, polygon region

#### `resolved(t)` translation contract

For a `CompositeObject` with origin `(ox, oy)`, rotation `theta`, and motion
`MotionSpec` that gives `(dx, dy, dtheta)` at time `t`:
1. Apply motion: effective origin = `(ox + dx(t), oy + dy(t))`,
   effective theta = `theta + dtheta(t)`
2. For each child primitive, rotate its local `(x0, y0)` by effective theta,
   add effective origin → world-frame position
3. Propagate `_native_motion` (velocity + omega + r_ref) from composite to children

The returned dicts map directly to `contact_primitives.py` constructors and
`make_prim_data()` — no further translation needed.

### Test: `src/epd/tests/test_4_1.py`

```
1. Box(width=10, height=8, x0=5, y0=4, exclusion='exterior').resolved(t=0)
   → 4 primitives, normals pointing inward (verified by dot product with inward vector)
   → region_polygon() has 4 vertices, all points inside box pass point-in-polygon test
   → 10 random exterior points all fail

2. CouetteCell(inner_radius=2, outer_radius=5, x0=0, y0=0).resolved(t=0)
   → 2 arc primitives, correct radii and centers
   → region_polygon() annulus: point at r=3 is inside, r=1 and r=6 are outside

3. Box with set_motion(MotionSpec(vx=1.0)).resolved(t=2.0)
   → origin shifted by 2.0 in x
   → all 4 wall primitives shifted accordingly

OUTPUT: PNG showing Box and CouetteCell boundary polygons with test points
colored green (allowed) / red (excluded).
```

---

## Phase 4.2: `particles.py` — ParticleSpec ⬜ PENDING

### File: `src/epd/particles.py`

#### `ParticleSpec`

```python
ParticleSpec(
    count       : int,            # number of particles in this population
    type        : str = 'elastic', # 'elastic' | 'emulsion' | 'rigid'
    # -- high-level physical params (preferred) --
    nu          : float = None,   # effective Poisson ratio → derives q
    Bo          : float = None,   # Bond number (emulsion only) → derives gamma/K_area
    # -- low-level params (override if specified, else derived from nu/Bo) --
    q           : float = None,   # K_area / El_t; if None, derived from nu
    tau_b       : float = 0.2,    # bending working point
    alpha_damp  : float = None,   # if None, auto from T_wave
    C           : float = None,   # contact hardness; if None, 3000*S*(1+q)
    # -- size distribution --
    N_nodes     : int   = 60,     # perimeter nodes
    R0_mean     : float = 1.0,    # mean radius (always normalized to 1.0)
    poly_dist   : str|dict = None, # None=mono, 'gaussian_0.05', or
                                  # {'type':'gaussian','sigma':0.05} etc.
    # -- motion (set after construction) --
    motion      : MotionSpec = None,
    frozen_shape: bool = False,
    # -- extensible force parameters (all default zero = no-op) --
    extra_forces: dict = None,    # {'drag': gamma, 'activity': v0, ...}
)
```

#### Material derivation rules

```
type='elastic':
  if nu given:   q = log_interp(nu, calibration_table)
  elif q given:  use q directly
  else:          q = 2.0 (default)
  TAU      = sqrt(12 * tau_b)
  El_t     = 12 / TAU**2              (= S=1, R0=1)
  K_area   = q * El_t
  C        = C if given else 3000 * S * (1 + q)
  alpha    = alpha_damp if given else 2.0  (Q≈3 default)

type='emulsion':
  Bo given → surface tension gamma = rho*g*R0^2 / Bo
  K_area from Laplace pressure
  C, alpha as above

type='rigid':
  frozen_shape forced True
  C = large (e.g. 1e6) for hard contact
  K_area, El_t irrelevant but set to defaults
```

#### Size distribution

```python
poly_dist = None            → all R0 = 1.0
poly_dist = 0.05            → Gaussian sigma=0.05 (shorthand)
poly_dist = {'type': 'gaussian', 'sigma': 0.05}
poly_dist = {'type': 'bimodal', 'ratio': 0.5, 'delta': 0.1}
poly_dist = {'type': 'explicit', 'values': [0.9, 1.0, 1.1, ...]}  # length must == count
```
Always normalized so mean=1.0 exactly after sampling.

#### `build(seed=42)` → `list[CapsuleParticle]`

Instantiates `CapsuleParticle` objects with derived parameters. Called by
`System.initialize()` — user never calls this directly.

#### `extra_forces` extensibility

`extra_forces={'drag': 0.1}` stores `{'drag': 0.1}` in the spec. When
`step_full_tf` is called, Phase 4 builds a `force_extras` dict and passes it in.
The TF kernel checks `force_extras.get('drag', 0.0)` and applies viscous drag
`F_drag = -gamma * v_node` if non-zero. When `drag` key is absent or zero: no-op.
**No existing call sites break.**

### Test: `src/epd/tests/test_4_2.py`

```
1. spec = ParticleSpec(count=5, type='elastic', nu=0.5)
   → verify q matches calibration table at nu=0.5 (within 1%)
   → verify TAU = sqrt(12*0.2), K_area = q * El_t, C = 3000*(1+q)

2. spec = ParticleSpec(count=5, type='elastic', q=2.0, tau_b=0.2)
   → verify same result as above (q given directly)

3. spec = ParticleSpec(count=20, poly_dist=0.05, seed=42)
   → R0_arr mean = 1.0 exactly
   → std/mean ≈ 0.05 (within 10%)
   → all R0 > 0

4. spec1 = ParticleSpec(count=10, nu=0.3)   # stiff
   spec2 = ParticleSpec(count=10, nu=0.8)   # soft
   → spec1.q < spec2.q (monotone)

OUTPUT: printed table of derived parameters for each test case.
```

---

## Phase 4.3: `motion.py` — MotionSpec ⬜ PENDING

### File: `src/epd/motion.py`

#### `MotionSpec`

Wraps a time-dependent motion description. Two modes:

**Mode 1 — TF-native callable (preferred):**
```python
MotionSpec(vx=lambda t: tf.sin(omega*t),
           vy=0.0,
           omega=0.0,
           r_ref=(x0, y0))
```
Each component is either a scalar constant or a TF-traceable callable `f(t) → scalar`.
`resolve_tf(t)` → `(vx, vy, omega, r_ref)` as TF tensors, evaluated at time `t`.

**Mode 2 — pre-sampled (Python lambda fallback):**
```python
MotionSpec.from_samples(vx_fn=lambda t: np.sin(omega*t),
                        dt=0.001, duration=10.0)
```
Pre-samples at `t = 0, dt, 2dt, ...`. `resolve_tf(t)` interpolates using
`tf.gather` + linear interpolation. Wraps/clamps at end of sample array.

**Constant shorthand:**
```python
MotionSpec(vx=1.0, vy=0.0)            # constant translation
MotionSpec(omega=0.5, r_ref=(5, 5))   # constant rotation about pivot
```

#### `resolve_tf(t)` contract

Returns `(vx_tf, vy_tf, omega_tf, r_ref_tf)` as `tf.float64` scalars/tensors.
All downstream code (prim_data assembly, driven particle injection) uses this method.
`t` is the absolute simulation time tensor from the state dict.

#### Integration with `make_prim_data`

`make_prim_data` currently takes static `vel` and `omega` arrays. Phase 4 wraps
this: before each step (outside TF loop), `prim_data` is rebuilt from the current `t`
if any object has time-varying motion. If all motions are constant (or TF-native),
`prim_data` is built once and reused.

**Optimization rule:** If all object motions are TF-native callables, they are
baked into the TF graph and `prim_data` never needs Python-side rebuild.
If any motion is pre-sampled, `prim_data` is rebuilt from Python at each
candidacy-update interval (every ~50 steps, already happening).

#### Driven particle motion

`ParticleSpec.motion = MotionSpec(vx=1.0)` with `frozen_shape=True` maps to:
- `driven_mask[i] = True` for all particles in this spec
- At each step, `trajectory[i] = motion.resolve_tf(t)` → `(vx, vy, omega_cm)`
- Fed to `step_rb_tf` as the driven-particle channel

### Test: `src/epd/tests/test_4_3.py`

```
1. MotionSpec(vx=1.0).resolve_tf(t=5.0) → vx=1.0, vy=0.0, omega=0.0

2. MotionSpec(vx=lambda t: tf.sin(t)).resolve_tf(t=pi/2) → vx≈1.0

3. MotionSpec.from_samples(vx_fn=lambda t: np.cos(t), dt=0.01, duration=10.0)
   .resolve_tf(t=pi) → vx≈-1.0 (within 1% for dt=0.01)

4. Wall with MotionSpec(omega=1.0, r_ref=(0,0)).resolved(t=1.0)
   → wall center has rotated 1 radian about origin

OUTPUT: plot of vx(t) for a sine MotionSpec vs analytic, showing interpolation error.
```

---

## Phase 4.4: `initializer.py` — RSA + Adaptive Swell ⬜ PENDING

### File: `src/epd/initializer.py`

#### RSA seeder

```python
rsa_seed(specs, objects, Lx, Ly, phi_init=0.30, seed=42, max_attempts=50000)
  → list[CapsuleParticle], placed_positions
```

Places particles one by one:
1. Sample candidate CM position uniformly in `[0,Lx] × [0,Ly]`
2. Reject if CM is in an excluded region (`_point_excluded(x, y, objects)`)
3. Reject if outer perimeter overlaps any already-placed particle
   (CM distance < `R0_i + R0_j + 2*r_c_i + 2*r_c_j`)
4. Reject if outer perimeter is too close to any object primitive
5. Accept → record position

**Accessible area** computed from object exclusions (same logic as Reference Code):
- `exterior` exclusion: particles must be inside; accessible area = enclosed area
- `interior` exclusion: particles excluded from inside; accessible area -= enclosed area
- Used to compute target `N_particles` from `phi_init` if not specified

**R0 normalization:** After sampling all sizes from poly_dist, normalize so
`mean(R0_arr) = 1.0` exactly. Applied before seeding.

#### Adaptive swell

Promotes the logic from `swell_adaptive.py` into a callable:

```python
adaptive_swell(state, params, cm_mgr, prim_data,
               phi_target, Lx, Ly,
               dphi_init=0.002, dphi_max=0.018, dphi_min=0.0001,
               n_relax=1500, max_extra_relax=3000,
               f_warn=0.45, f_crit=0.65,
               V_CRIT_FRAC=50.0, max_restores=6,
               verbose=True)
  → state, Lx, Ly   (updated)
```

Returns when `phi_outer >= phi_target`. Sets `state['t'] = 0.0` after completion.

#### `_point_excluded(x, y, objects)` → bool

For each object:
- If `exclusion='interior'`: reject if point is inside `region_polygon()`
- If `exclusion='exterior'`: reject if point is outside `region_polygon()`
- For objects with no `region_polygon()` (single walls): use primitive clearance check

### Test: `src/epd/tests/test_4_4.py`

```
1. Monodisperse P=16, fully periodic, phi_target=0.75
   → phi_outer >= 0.74
   → all CMs in [0,Lx]×[0,Ly]
   → no particle pair overlapping (min CM dist > 2*R_eff - tolerance)
   → render snapshot → image saved

2. P=16 inside Box(width=12, height=12, exclusion='exterior')
   → all CMs inside box region_polygon
   → phi_outer measured relative to box area

3. P=16 around CircleObstacle(radius=2, exclusion='interior') in periodic box
   → no CM within radius 2 of obstacle center

OUTPUT: PNG snapshots for all 3 cases.
```

---

## Phase 4.5: `system.py` + `checkpoint.py` ⬜ PENDING

### File: `src/epd/system.py`

#### `System` class

```python
System(Lx, Ly,
       periodic_x=True, periodic_y=True,
       dt_factor=0.4,      # multiplier on auto-stable dt
       alpha_damp=None,    # if None, auto from T_wave of first spec
       E_candidates=128,   # candidacy table width
       skin=1.0)           # candidacy skin factor × R0_mean
```

Public methods:

```python
.add_particles(spec: ParticleSpec) → self
.add_object(obj: SimulationObject) → self
.initialize(phi_target, seed=42, verbose=True) → self
  # RSA + adaptive swell; sets state['t']=0, state['step']=0
.step(N=1) → self
  # N dynamics steps; updates state in place; candidacy updated as needed
.run(N, record_every=None, output_dir=None, progress_fn=None) → self
  # N steps; if record_every: save snapshot every K steps to output_dir
.save(path) → Path
.load(path) → self    # fully self-contained; no re-registration needed
.freeze(zero_velocity=True) → self
  # zero v_cm, omega, u_dot; useful before production run
.render(ax=None, output_path=None, t_label=True) → Figure
.make_gif(output_path, frame_paths) → Path
```

Properties:
```python
.t         → float (current absolute time)
.step_count → int
.phi_outer  → float (current actual packing fraction)
.phi_box    → float (Lx*Ly-based packing fraction)
.particles  → list[CapsuleParticle] (read-only)
.state      → dict (full state; read-only reference)
```

#### `step(N)` implementation

```python
def step(self, N=1):
    for _ in range(N):
        # 1. Rebuild prim_data if any object motion is time-varying and t changed
        if self._prim_needs_rebuild():
            self._prim_data = self._build_prim_data(self.state['t'])
        # 2. Candidacy update (Python, every ~50 steps or on threshold)
        if self._cm_mgr.needs_update(x_cm, theta):
            self._cm_mgr.update(x_cm, theta)
        # 3. TF step
        self.state, _ = step_full_tf(
            self.state, self._cm_mgr.CapCandidates,
            self._dt_tf, self._alpha_tf, self._g_tf, self._params,
            t=tf.constant(self.state['t'], dtype=DTYPE),
            prim_data=self._prim_data,
            frozen_mask=self.state['frozen_mask'],
            driven_mask=self._driven_mask,
            trajectories=self._eval_trajectories(self.state['t']),
        )
        # 4. Wrap periodic coordinates
        self.state = self._wrap(self.state)
        # 5. Advance time
        self.state['t'] += self._dt
        self.state['step'] += 1
```

### File: `src/epd/checkpoint.py`

#### Checkpoint format

```
path/
  state.npz      ← all numpy arrays from state dict
                    (x_all, x_cm, v_cm, omega, theta, u, u_dot, X_ref,
                     frozen_mask, CapCandidates)
  config.json    ← {Lx, Ly, periodic_x, periodic_y, dt, alpha_damp,
                     E_candidates, skin, t, step}
  particles.json ← list of ParticleSpec dicts (serialized, includes motion)
  objects.json   ← list of Object dicts (serialized, includes motion)
```

`load()` reconstructs:
1. `CapsuleParticle` list from `particles.json` (calls `spec.build()`)
2. `CandidacyManager` from positions in `state.npz` (not from saved CapCandidates —
   candidates are rebuilt fresh to avoid stale data; saved only for inspection)
3. Motion specs from `particles.json` / `objects.json`
4. `prim_data` from objects at `t=state['t']`
5. All TF constants (dt, alpha, params) from `config.json`

#### Motion serialization

TF-native callables cannot be pickled. Serialization rules:
- Constant motion → `{'mode': 'constant', 'vx': ..., 'vy': ..., 'omega': ..., 'r_ref': ...}`
- TF-native → user must provide a string tag + parameters:
  `{'mode': 'sine', 'axis': 'x', 'amplitude': 1.0, 'frequency': 0.5}`
  Supported built-in profiles: `constant`, `sine`, `cosine`, `linear`, `square`
- Pre-sampled → samples array saved as part of `state.npz`

### Test: `src/epd/tests/test_4_5.py`

```
1. Initialize P=25 polydisperse, phi=0.80. Run 200 steps. Save.
   Load. Run 200 more steps. Compare to reference that ran 400 steps uninterrupted.
   Gate: max|Δx_cm| < 1e-12 (bit-identical resume)

2. Save with a sine-wave wall motion. Load. Verify wall velocity at t=T/4
   matches expected value from MotionSpec.

3. system.t before and after save/load cycle → identical float64 values.

OUTPUT: printed diff table (step, max|Δx|) confirming bit-identical resume.
```

---

## Phase 4.6: Integration Verification ⬜ PENDING

### File: `src/epd/tests/test_4_6.py`

Six sub-tests, each producing a GIF or PNG:

**4.6.A — Dense bulk packing (baseline)**
```python
sys = System(Lx=12, Ly=12)
sys.add_particles(ParticleSpec(count=25, nu=0.5, poly_dist=0.05))
sys.initialize(phi_target=0.85)
sys.run(2000)
sys.render(output_path='4_6A_dense_packing.png')
```
Gate: phi_outer ≥ 0.84, max_f < 0.05, min_circ > 0.95.

**4.6.B — Oscillating wall pushes particles**
```python
sys = System(Lx=12, Ly=12, periodic_x=False, periodic_y=False)
wall = Wall(...)
wall.set_motion(MotionSpec(vy=lambda t: tf.sin(2*t)))
sys.add_object(wall)
sys.add_particles(ParticleSpec(count=10, nu=0.5))
sys.initialize(phi_target=0.60)
sys.run(500, record_every=20, output_dir='4_6B_frames/')
sys.make_gif('4_6B_oscillating_wall.gif', ...)
```
Gate: wall position at t=π/2 has moved by expected amount; particles show
corresponding displacement.

**4.6.C — Driven frozen particle (rigid pusher)**
```python
spec = ParticleSpec(count=1, type='rigid', frozen_shape=True)
spec.set_motion(MotionSpec(vx=0.5))
sys.add_particles(spec)
sys.add_particles(ParticleSpec(count=15, nu=0.5))
sys.initialize(phi_target=0.65)
sys.run(400, record_every=20, output_dir='4_6C_frames/')
sys.make_gif('4_6C_rigid_pusher.gif', ...)
```
Gate: pusher CM moves at vx=0.5 exactly (verify displacement = 0.5*t);
pushed particles are displaced; pusher shape u=0 throughout.

**4.6.D — Couette cell**
```python
sys = System(Lx=14, Ly=14)
cell = CouetteCell(inner_radius=2, outer_radius=6, exclusion='exterior')
sys.add_object(cell)
sys.add_particles(ParticleSpec(count=20, nu=0.5))
sys.initialize(phi_target=0.70)
# Drive top/bottom particle layers
top = sys.select_particles('top_fraction', fraction=0.15)
bot = sys.select_particles('bottom_fraction', fraction=0.15)
sys.particles_set_motion(top, MotionSpec(vx=0.3, frozen_shape=True))
sys.particles_set_motion(bot, MotionSpec(vx=-0.3, frozen_shape=True))
sys.run(600, record_every=20, output_dir='4_6D_frames/')
sys.make_gif('4_6D_couette.gif', ...)
```
Gate: all CMs remain inside annulus throughout; top particles move at vx=0.3.

**4.6.E — Save/load mid-run (regime change)**
```python
sys.initialize(phi_target=0.75)
sys.run(300)
sys.save('checkpoint_E/')
sys2 = System.from_file('checkpoint_E/')
sys2.run(300)
sys.run(300)   # continue original
# sys and sys2 must agree
```
Gate: max|Δx_cm| < 1e-12 between sys and sys2 after 300 more steps.

**4.6.F — Mixed elastic + emulsion**
```python
sys.add_particles(ParticleSpec(count=10, type='elastic', nu=0.5))
sys.add_particles(ParticleSpec(count=10, type='emulsion', Bo=0.01))
sys.initialize(phi_target=0.70)
sys.run(500)
sys.render(output_path='4_6F_mixed.png')
```
Gate: no crash; both particle types coexist; emulsion particles show
characteristic deformation distinct from elastic.

---

## Phase 4.7: Systematic Test Validation → Colab Notebook Assembly 🔄 ACTIVE

### Purpose and motivation

This phase has two intertwined goals:

1. **Validate the C++/TF backend and System API systematically**, test by test, with
   a clear chain of ground truth: NumPy reference → TF direct → API → Colab.

2. **Build a "Getting Started" Google Colab notebook** (`notebooks/getting_started.ipynb`)
   that documents working, confirmed test cases as self-contained sections. Each section
   has a markdown explanation cell followed by a runnable code cell.

**Why the Colab notebook matters:** Development permanently lives on this CPU-only machine.
A future phase (est. Phase 5 or 6) will make the repo live and installable. At that point,
a Colab session can pull the repo, install, and run the notebook on GPU — capturing results,
logs, and GIFs that are uploaded back to the repo for review here. This closes the
CPU↔GPU feedback loop without requiring a local GPU. The Colab notebook is the artifact
that makes GPU testing possible.

A second Colab notebook (for the methods paper) will eventually collect all
validation/calibration scripts so readers can reproduce every experiment. That is a
later deliverable; do not start it now.

### Workflow (repeats for each new test)

```
1. NumPy reference   — establishes physics ground truth; user approves GIF
2. TF direct         — port to C++/TF; patch tf_sim.py until results agree; user approves
3. API version        — re-implement via System API; patch API layer as needed; user approves
4. Colab section      — add approved API code as a section in getting_started.ipynb
                        (markdown explanation cell + code cell)
```

The user dictates tests one at a time. There is no predefined test list for this phase.

### Approved tests (cumulative)

| Test | Description | NumPy | TF direct | API | Colab |
|------|-------------|-------|-----------|-----|-------|
| A-elastic | Two-disk oscillatory squeeze, N=32 | ✅ | ✅ `test_A_elastic_tf.py` (~85s) | ⬜ | ⬜ |
| A-emulsion | Two-droplet oscillatory squeeze, N=32, κ=0.02, k_reg=10 | ✅ | ✅ `test_A_emulsion_tf.py` (~80s) | ⬜ | ⬜ |

### Patches applied this phase

| File | Change | Reason |
|------|--------|--------|
| `tf_sim.py` | Added `k_reg_forces_tf()` + `k_reg_per_p` param; excluded from RB sums in `step_rb_tf` | Tangential regularization keeps emulsion nodes evenly spaced |

### Immediate next step

User will dictate the next test, or proceed to API version of Test A.

---

## Phase 4 Exit Gates (4.1–4.6 — legacy)

These gates applied to the 4.1–4.6 API scaffolding work, which is complete.
Phase 4.7 has its own incremental approval workflow above.

- [x] All 4.1–4.6 tests pass and produce expected images/GIFs
- [x] `test_4_5` checkpoint resume is bit-identical (max|Δx| < 1e-12)
- [x] `test_4_6C` driven particle displacement matches 0.5*t to < 0.1%
- [x] `test_4_6D` all CMs remain inside CouetteCell annulus throughout


---


## Phase 4.8: Per-Particle R0 Scaling, Parameter Recompute, TF Graph Fix 🔄 ACTIVE

### Motivation

Two related correctness gaps, both now addressed in one phase:

1. **Emulsion `build()` scaling bugs** — K_area was scaled by R0² (should be R0⁰) and C by
   R0¹ (should be R0⁰) to preserve κ and C̃ as R0 varies. These are invisible at 5%
   polydispersity but significant at size ratios ≥1.2×.

2. **alpha_damp / xi_drag are per-system scalars** — to preserve COR (α ∝ R0⁻¹) and terminal
   velocity (ξ ∝ R0⁺¹) across polydisperse particles they must be per-particle tensors of
   shape (P,) in the TF graph. This also subsumes the known 10-retraces/batch retrace issue.

Elastic particles (`CapsuleParticle`) were already R0-correct: the constructor scales El_t,
EI, k_c from the passed R0. Only the emulsion path and the per-system damping globals needed
fixing.

### Waypoints

| # | Task | Status |
|---|------|--------|
| 4.8.1 | Fix emulsion `build()` exponents: K_area R0²→R0⁰, C R0¹→R0⁰ | ⬜ |
| 4.8.2 | Stamp `_kappa_target`, `_Oh_target`, `_nu_target` on particles at build time | ⬜ |
| 4.8.3 | `p.R0` property setter → `_recompute_params(r0)` on CapsuleParticle & EmulsionParticle | ⬜ |
| 4.8.4 | `System.adjust_params_for_size()` — reads actual R0 from geometry, fires setter | ⬜ |
| 4.8.5 | Per-particle alpha_damp & xi_drag: store on particles; assemble as (P,) tensors; update integrator; scalar input broadcasts | ⬜ |
| 4.8.6 | Fix retrace issue: 10 retraces/batch → 1 (graph signature stable across `run_fast` calls) | ⬜ |
| 4.8.7 | Fix scaling-relation errors in `papers/summary_of_methods/` if any | ⬜ |
| 4.8.8 | Validation spot-check: two-population system (R0=0.714 and R0=1.4) — confirm κ, Oh, ν preserved after `adjust_params_for_size()` | ⬜ |

### Design rules

- **Backwards compatible:** scalar `alpha_damp` / `xi_drag` in any existing call is broadcast
  to shape (P,) internally; no existing test or notebook cell changes.
- **TF graph shape contract:** `alpha_damp_p` and `xi_drag_p` always assembled as fixed-shape
  `(P,)` tensors before graph entry so the signature never changes and retrace count stays at 1.
- **Idempotent recompute:** calling `adjust_params_for_size()` twice gives identical results.
- **Elastic already correct:** CapsuleParticle constructor scales El_t ∝ R0, EI ∝ R0³,
  k_c ∝ R0⁻¹ at construction. Property setter only needs to redo those same lines in-place.

### Exit gate

- [ ] 4.8.1–4.8.4: `adjust_params_for_size()` spot-check passes (κ, Oh, ν within 1% of target
  for both R0=0.714 and R0=1.4 particles after recompute)
- [ ] 4.8.5: per-particle alpha_damp/xi_drag — existing tests (A-elastic, A-emulsion, D-elastic,
  D-emulsion, F-hopper) all pass unchanged
- [ ] 4.8.6: retrace count = 1 per `run_fast` call (not 10)
- [ ] All existing Colab notebook cells produce identical output

---

## Phase 5: Stokes Drag + Parameter Update Infrastructure ⬜ PENDING

**Branch:** `stokes-drag-dev` (branched fresh from `tf-fast` tip before any edits)
**Freeze:** `tf-fast` is tagged/frozen at Phase 4.7 completion before this branch is created.
All existing test scripts and the colab notebook must continue to pass unchanged throughout
this phase. Stokes drag is strictly opt-in (Oh=None → zero drag, backward compatible).

---

### Science goal

Add physically grounded external Stokes drag to the EPD/emulsion model as a new
node-level force, parameterised by the dimensionless Ohnesorge number Oh. Simultaneously
introduce a clean parameter update mechanism so any dimensionless parameter (Bo, κ, ν, Oh,
or lower-level equivalents) can be changed post-initialisation at per-particle granularity,
with automatic bidirectional propagation through the derivation chain. Both additions are
demonstrated and verified through a suite of four benchmark tests (Tests E1–E4).

This phase proves the architecture's extensibility claim: a new force and a new
parameter can be added without restructuring the integrator, without breaking backward
compatibility, and without requiring users to reason about internal TF tensor layout.

---

### Dimensionless parameter — Ohnesorge number

The Ohnesorge number Oh compares viscous drag to the particle's dominant restoring force,
independent of gravity:

| Type     | Reference velocity           | ξ from Oh           | Physical meaning         |
|----------|------------------------------|---------------------|--------------------------|
| Emulsion | v_cap = √(γ / ρ_d R₀) = 1   | ξ = Oh              | drag vs capillary forces |
| Elastic  | v_el  = √(El_t / ρ_d R₀)    | ξ = Oh · √(El_t)    | drag vs membrane stiffness |

With γ = ρ_d = R₀ = 1 in our units, Oh = η_ext (external fluid viscosity per unit length).
Oh=1 means viscous drag comparable to the restoring force in both model types.

Terminal velocity is derived, not set directly:

```
v_terminal = g / (2 · ξ)   [emulsion, R₀=1, ρ_d=1]
           = g / (2 · Oh)   [emulsion working point]
```

Convenience helpers allow the user to work in whichever direction is natural:

```python
# Physicist path: set Oh, read terminal velocity
spec = ParticleSpec(..., Oh=0.5)
print(spec.terminal_velocity(g=0.05))      # → 0.05

# Engineer path: target terminal velocity, get Oh
from src.epd.drag import oh_from_terminal_velocity
Oh = oh_from_terminal_velocity(v_t=0.05, g=0.05, particle_type='emulsion')

# Override path: set terminal velocity, Oh updated automatically
spec.set_terminal_velocity(v_t=0.05, g=0.05)   # back-computes ξ → Oh
print(spec.derived['Oh'])                        # consistent
```

---

### Drag force formulation

Resistive Force Theory (Lauga & Powers 2009; Pozrikidis 1992) applied to the
closed membrane perimeter. No projection — full vector drag opposing all relative
motion between node and background flow:

```
f_drag,i = −ξ · (v_node,i − U_bg(x_i, t)) · ΔL_i      [force on node i]
```

where:

- `v_node,i = v_cm + ω × r_i + R(θ) · u_dot_i`
  Reconstructed each step from existing state (v_cm, omega, u_dot, theta).
  No new state variables required.

- `ΔL_i = |x_{i+1} − x_{i-1}| / 2`
  Arc-length weight at node i. For near-circular particles ΔL_i ≈ L₀/N.

- `U_bg(x_i, t)` — background flow evaluated at node position (see presets below)

The summed drag decomposes automatically through the existing rigid-body solver into:
- Net CM force (from the v_cm component of v_node)
- Torque (from the ω × r component — larger particles experience more rotational drag)
- Elastic damping of internal modes (from the R·u_dot component)

This is orthogonal to the existing α_damp term, which damps body-frame elastic
oscillations only and contributes zero net CM force.

---

### Background flow — TF-native presets

To remain inside the `tf.while_loop` and preserve the tf-fast optimisation,
background flow is encoded as an integer type flag plus a parameter vector
in `params`. No Python callbacks per step. All flows are pure TF tensor ops.

| Preset          | U_bg(x, y, t)                        | Params stored          |
|-----------------|--------------------------------------|------------------------|
| `'zero'`        | (0, 0)                               | none                   |
| `'constant'`    | (U₀x, U₀y)                          | U_vec (2,)             |
| `'shear'`       | (γ̇ · y, 0)                          | shear_rate scalar      |
| `'parabolic'`   | (U_max · (1 − (y/H)²), 0)           | U_max, H scalars       |
| `'extensional'` | (ε̇ · x, −ε̇ · y)                    | ext_rate scalar        |

User-facing API (set after initialisation, takes effect immediately):

```python
sys.U_background = None                                       # quiescent (default)
sys.U_background = ('constant',    {'U': (0.1, 0.0)})
sys.U_background = ('shear',       {'rate': 0.05})
sys.U_background = ('parabolic',   {'U_max': 0.2, 'H': 10.0})
sys.U_background = ('extensional', {'rate': 0.01})
```

Arbitrary Python callables are explicitly out of scope for this phase (would require
tf.py_function per step, destroying the while_loop optimisation). Documented as a
future extension in the code.

---

### Parameter hierarchy and bidirectional propagation

All material parameters sit in a two-level hierarchy with diagnostics below:

```
High-level dimensionless (user-facing):
    Bo, κ (emulsion)  |  ν, q (elastic)  |  Oh (both)
           ↕  bidirectional — setters walk chain in both directions
Mid-level physical (internal, derived):
    g, K_area, El_t, EI, C, ξ (drag per arc length)
           ↕
Diagnostic read-only (computed from mid-level):
    v_terminal, dt, alpha_eff, c_cap
```

Rules:
- Setting a high-level param recomputes all mid-level dependents downward,
  then recomputes all other high-level params upward for consistency.
- Setting a mid-level param recomputes dependent mid-level params, then
  recomputes all high-level params upward.
- Diagnostic properties are always computed on-the-fly, never cached.
- All existing param derivation paths are preserved unchanged (backward compat).

Example cascade when user calls `spec.set_kappa(0.01)`:

```
κ = 0.01
  → K_area = γ / (κ · R₀)          [mid-level down]
  → q = 1/κ = 100                   [mid-level down]
  → C = 500 · γ / R₀                [mid-level, C independent of κ — unchanged]
  → v_terminal = g / (2·ξ)          [diagnostic, recomputed]
  → if system initialized: patch params['K_area'][particle_ids] in TF tensors
```

Example cascade when user calls `spec.set_xi(0.025)`:

```
ξ = 0.025
  → Oh = ξ / v_ref                  [high-level up]
  → v_terminal = g / (2·ξ)          [diagnostic]
  → if system initialized: patch params['xi_drag_per_p'][particle_ids] in TF tensors
```

---

### Per-particle update granularity

All material params in the TF `params` dict are already (P,) tensors (K_area, El_t,
gamma_lt, r_c_per_p, etc.). The new `xi_drag_per_p (P,)` tensor follows the same
pattern. Updates target specific particle indices:

```python
sys.set_param('Oh',    0.5,  particles='all')        # all particles
sys.set_param('kappa', 0.01, particles=[0, 1, 5])    # specific indices
sys.set_param('nu',    0.6,  particles=range(10))    # slice
sys.set_param('xi',    0.02, particles='all')        # mid-level override
```

Mixed-population example (some emulsion, some elastic):

```python
spec_em  = ParticleSpec(count=20, type='emulsion', kappa=0.02, Oh=0.3)
spec_el  = ParticleSpec(count=10, type='elastic',  nu=0.5,     Oh=0.5)
sys.add_particles(spec_em)
sys.add_particles(spec_el)
# After swell, update just the emulsion population:
sys.set_param('Oh', 0.1, particles=range(20))
```

`System` owns the particle-index bookkeeping and the TF params dict, making it
the natural owner of `set_param`. It determines the particle type (emulsion vs elastic)
for each index to apply the correct ξ derivation.

---

### Pre-init vs post-init behaviour of set_param

**Before `sys.initialize()`:**
- `set_param` updates `ParticleSpec` internal state only.
- No TF tensors exist yet; values are consumed when `initialize()` builds params.
- Useful for assembling a system before running.

**After `sys.initialize()`:**
- `set_param` updates `ParticleSpec` internal state AND immediately patches the
  affected (P,) TF tensor slice in `self._params` in-place.
- Simulation resumes at the new parameter values on the very next step.
- Particle positions, velocities, and all other state are preserved exactly.
- Useful for parameter continuation (e.g. swell at κ=0.02, then resume at κ=0.005
  from the same packed configuration without re-running RSA + swell).

---

### Branch and backward-compatibility strategy

1. **Freeze `tf-fast`:** before any edits, confirm all existing tests pass on
   `tf-fast` (Tests A, B2, C, D-emulsion, D-elastic). Tag the tip as `v2.0-tf-fast-frozen`.

2. **Create `stokes-drag-dev`** branched from that tag. All Phase 5 work happens here.

3. **Backward compatibility invariant:** Oh defaults to None → xi_drag=0 →
   drag force term evaluates to exactly zero. The `if xi_nonzero` TF branch is
   short-circuit zero-cost. All existing scripts run without modification.

4. **Merge policy:** `stokes-drag-dev` merges back into `tf-fast` only after all
   four Test E benchmarks pass and the exit gate is satisfied.

---

### Files changed

| File | Nature of change | Backward compatible |
|------|-----------------|---------------------|
| `src/simulation/tf_sim.py` | Add `xi_drag_per_p (P,)` + `U_bg_type/U_bg_params` to params; add `v_node` reconstruction + arc-length-weighted drag force inside `step_rb_tf`; add TF-native preset evaluator function | ✅ new params default to zero/zero-preset |
| `src/epd/particles.py` | Add `Oh` kwarg to `__init__`; extend `_derive_material` to compute `ξ`; extend `derived` dict; add `terminal_velocity(g)` method; add `set_terminal_velocity(v_t, g)` method; add setters for cascading updates | ✅ Oh=None default |
| `src/epd/system.py` | Add `U_background` property with setter; add `set_param(name, value, particles)` method; pass `xi_drag_per_p` and U_bg params into `step_full_tf` and `run_simulation_tf` | ✅ default None/zero |
| `src/epd/drag.py` (new) | `oh_from_terminal_velocity(v_t, g, particle_type, **material)` standalone helper | ✅ new file |
| `src/validation/test_E_stokes_drag.py` (new) | Four benchmark tests E1–E4 | ✅ new file |

---

### Bucket A: Core drag force in tf_sim.py

#### Waypoint 5A.1 — Add xi_drag_per_p to params ⬜
- Add `xi_drag_per_p: (P,)` tensor to `make_state()`, defaulting to zeros
- Add `U_bg_type: scalar int32` (0=zero, 1=constant, 2=shear, 3=parabolic, 4=extensional)
- Add `U_bg_params: (4,) float64` storing preset parameters (padded to fixed length)
- **Gate:** `make_state()` produces identical results to before for all existing tests
  (xi_drag_per_p=0, U_bg_type=0 → no change to any force)

#### Waypoint 5A.2 — TF-native background flow evaluator ⬜
New function `_eval_U_bg(x_nodes, t, params)` → `(P, N, 2)` TF tensor.
Dispatches on `params['U_bg_type']` using `tf.switch_case` or nested `tf.cond`:
- type 0: return zeros
- type 1: broadcast constant U_vec to all nodes
- type 2: shear — U_x = rate * y, U_y = 0
- type 3: parabolic — U_x = U_max * (1 − (y/H)²), U_y = 0
- type 4: extensional — U_x = rate * x, U_y = −rate * y
- **Gate:** unit tests for each preset at fixed node positions; values match
  analytic formula to machine precision

#### Waypoint 5A.3 — Node velocity reconstruction ⬜
Inside `step_rb_tf`, before force decomposition, compute:
```
r_i     = x_all − x_cm[:, None, :]
v_rot_i = omega[:, None, None] * stack([-r[:,:,1], r[:,:,0]])
v_el_i  = R(theta) applied to u_dot  (lab-frame elastic velocity)
v_node  = v_cm[:, None, :] + v_rot_i + v_el_i
```
- **Gate:** for a rigidly translating particle (u_dot=0, omega=0), v_node = v_cm
  for all nodes (machine precision). For a purely rotating particle (v_cm=0,
  u_dot=0), |v_node,i| = omega * |r_i| for all nodes.

#### Waypoint 5A.4 — Arc-length weights + drag force ⬜
```
dL_i    = |x_{i+1} − x_{i−1}| / 2            (P, N) tensor
f_drag  = −xi_drag_per_p[:, None, None]
          * (v_node − U_bg) * dL_i[:, :, None]
```
Add f_drag to f_contact before passing into the force decomposition.
- **Gate (zero drag):** xi_drag_per_p=0 → f_drag = 0 exactly → all existing
  test metrics unchanged to machine precision
- **Gate (quiescent):** single particle, constant v_cm, no other forces.
  Net drag force = −ξ · v_cm · ∮ dL_i = −ξ · v_cm · L_perimeter (within 0.1%)

---

### Bucket B: ParticleSpec — Oh parameter + cascading setters

#### Waypoint 5B.1 — Add Oh to ParticleSpec ⬜
- Add `Oh=None` kwarg to `__init__`
- In `_derive_material`: if Oh is not None, compute ξ = Oh (emulsion) or Oh·√(El_t) (elastic)
- Store `self._xi_drag` internally
- Extend `derived` dict: add `'Oh'`, `'xi'`
- **Gate:** `ParticleSpec(count=5, type='emulsion', gamma=1.0, kappa=0.02, Oh=0.5).derived`
  returns correct Oh=0.5, xi=0.5; `ParticleSpec(..., Oh=None).derived['xi'] == 0.0`

#### Waypoint 5B.2 — terminal_velocity method ⬜
```python
spec.terminal_velocity(g)            # returns g / (2 * xi) for emulsion
spec.terminal_velocity(g)            # returns g / (2 * xi / sqrt(El_t)) for elastic
```
Raises ValueError if Oh=None (no drag) or g=0.
- **Gate:** for emulsion Oh=0.5, g=0.05 → v_t = 0.05

#### Waypoint 5B.3 — set_terminal_velocity ⬜
```python
spec.set_terminal_velocity(v_t=0.05, g=0.05)
```
Back-computes ξ → Oh, stores both, cascades to xi_drag in derived.
- **Gate:** round-trip: set_terminal_velocity(v_t, g) → terminal_velocity(g) == v_t
  to machine precision

#### Waypoint 5B.4 — Cascading setters ⬜
Implement `set_Bo`, `set_kappa`, `set_nu`, `set_Oh`, `set_xi` on ParticleSpec.
Each setter:
1. Updates the named parameter
2. Recomputes all mid-level dependents downward
3. Recomputes all high-level dimensionless params upward for consistency
4. If `self._system` is set (post-init), calls `self._system._push_params(self)`
   to patch live TF tensors immediately

- **Gate:** `spec.set_kappa(0.01)` → `spec.derived['kappa'] == 0.01`,
  `spec.derived['K_area'] == 100`, `spec.derived['q'] == 100`
- **Gate:** `spec.set_xi(0.025)` → `spec.derived['Oh'] == 0.025` (emulsion, v_ref=1)

---

### Bucket C: System — set_param + U_background

#### Waypoint 5C.1 — set_param method on System ⬜
```python
sys.set_param(name, value, particles='all')
```
- Resolves `particles` to a list of integer indices
- Determines particle type (emulsion/elastic) for each index
- Calls appropriate setter on the owning `ParticleSpec` for those particles
- If post-init: immediately patches the corresponding (P,) slice of `self._params`
  using `tf.tensor_scatter_nd_update` or numpy array reassignment
- If pre-init: stores in ParticleSpec only; applied at initialize() time
- **Gate (pre-init):** set_param before initialize → correct params after initialize
- **Gate (post-init):** set_param after initialize → `self._params['K_area'][i]`
  reflects new value immediately; next step uses updated value

#### Waypoint 5C.2 — U_background property ⬜
```python
sys.U_background = None                                      # default
sys.U_background = ('constant', {'U': (0.1, 0.0)})
sys.U_background = ('shear',    {'rate': 0.05})
sys.U_background = ('parabolic', {'U_max': 0.2, 'H': 10.0})
```
Setter validates preset name, encodes type int + params vector, updates
`self._params['U_bg_type']` and `self._params['U_bg_params']` in-place.
Works pre- and post-init.
- **Gate:** setting U_background='zero' → no change to any test metric

#### Waypoint 5C.3 — Pass drag params through run pipeline ⬜
`xi_drag_per_p` and `U_bg_type/U_bg_params` already in params dict, which is
already passed to `step_full_tf` and `run_simulation_tf`. Verify no signature
changes are needed.
- **Gate:** existing `run()` and `run_fast()` calls unchanged; drag active only
  when xi_drag_per_p > 0

---

### Bucket D: Convenience helper

#### Waypoint 5D.1 — drag.py ⬜
New file `src/epd/drag.py`:
```python
def oh_from_terminal_velocity(v_t, g, particle_type='emulsion',
                               gamma=1.0, El_t=None, rho_d=1.0, R0=1.0):
    """
    Compute Oh such that terminal velocity equals v_t under gravity g.
    particle_type: 'emulsion' or 'elastic'
    """
```
- Raises ValueError if g == 0
- **Gate:** `oh_from_terminal_velocity(0.05, 0.05, 'emulsion') == 0.5`

---

### Bucket E: Test benchmarks

File: `src/validation/test_E_stokes_drag.py`

#### Test E1 — Terminal velocity (emulsion + elastic) ⬜
Setup: single particle, quiescent fluid (U_bg=zero), constant g.
Run until v_cm_y plateaus (|dv/dt| < 1e-4 for 500 steps).

Three (g, Oh) pairs for emulsion:
- (0.05, 0.25) → v_t_theory = 0.10
- (0.05, 0.50) → v_t_theory = 0.05
- (0.10, 0.50) → v_t_theory = 0.10

One pair for elastic (ν=0.5):
- (0.05, 0.5) → v_t_theory = g / (2 · Oh · √El_t)

**Gate:** |v_t_measured − v_t_theory| / v_t_theory < 5% for all cases.

Note: The existing α_damp also damps v_cm slightly (it feeds back through the
elastic modes). The test confirms the dominant drag comes from xi_drag as expected.

#### Test E2 — Particle tracking constant background flow ⬜
Setup: single emulsion particle, g=0, Oh=0.5.
U_background = ('constant', {'U': (0.1, 0.0)}).
Particle starts at rest.

Expected: v_cm_x → 0.1 exponentially, time constant τ_drag = M / ζ = ρ_d·π·R₀² / (2π·ξ·R₀).

**Gate:** |v_cm_x_final − 0.1| / 0.1 < 5% after 5 × τ_drag steps.

#### Test E3 — Shear flow drift and rotation ⬜
Setup: single emulsion particle, g=0, Oh=0.5.
U_background = ('shear', {'rate': γ̇}).
Particle starts at y_cm = H/2.

Expected (leading-order Jeffery orbit for a sphere):
- CM drifts in x at rate U_bg(y_cm) = γ̇ · y_cm
- Particle rotates at Ω ≈ γ̇/2

**Gate:**
- |dx_cm/dt − γ̇ · y_cm| / (γ̇ · y_cm) < 10% (drift rate)
- |omega_measured − γ̇/2| / (γ̇/2) < 20% (rotation rate; looser because elastic
  modes contribute to angular momentum budget)

#### Test E4 — Post-init parameter update ⬜
Setup: single emulsion particle, Oh=0.5, g=0.05.
Run 2000 steps → v_cm_y plateaus at v_t1 = g/(2·Oh) = 0.05.
Call `sys.set_param('Oh', 1.0, particles='all')`.
Run 2000 more steps → v_cm_y should plateau at v_t2 = g/(2·1.0) = 0.025.

**Gate:**
- v_t1 within 5% of 0.05 before update
- v_t2 within 5% of 0.025 after update
- No re-initialization called between the two runs

---

### Phase 5 Exit Gate

All of the following must pass before Phase 5 is marked complete and merged to `tf-fast`:

- [ ] All existing tests (A-elastic, A-emulsion, B2, B2-emulsion, C-emulsion,
      D-emulsion, D-elastic) pass unchanged on `stokes-drag-dev` with Oh=None default
- [ ] Waypoint 5A.4 zero-drag gate: xi=0 → identical physics to pre-Phase-5 code
- [ ] Test E1 PASS: terminal velocity within 5% for all 4 (g, Oh) pairs
- [ ] Test E2 PASS: constant background flow tracking within 5%
- [ ] Test E3 PASS: shear drift and rotation within tolerances
- [ ] Test E4 PASS: post-init set_param updates live TF params, new terminal velocity reached
- [ ] `spec.terminal_velocity(g)` and `spec.set_terminal_velocity(v_t, g)` round-trip exact
- [ ] `sys.set_param` pre-init and post-init gates both pass
- [ ] Colab notebook cells added for Test E — only after user approves benchmark results


---

## Test F — Emulsion Hopper Discharge  🔄 ACTIVE

**Objective:** Showcase simulation combining all EPD/emulsion capabilities:
gravity + drag + background flow + custom geometry (HopperRegion) + polydispersity.

**Geometry:** 30° funnel walls, 5 R₀ outlet, 12 R₀ wide reservoir, open top and bottom.
Two vertical side walls prevent particles escaping the sides.

**Cases:**
- F1: gravity g=0.05, Oh=0.15, no background flow
- F2: Poiseuille vertical (g=0, U_max=2×v_t≈0.33), Oh=0.15

**Parameters:** N=20 emulsion droplets, κ=0.02, N_nodes=36, 5% Gaussian polydispersity.

**Code changes for Test F:**
- `src/simulation/tf_sim.py`: added preset 5 (poiseuille_v): U_y(x)=-U_max*(1-(x-x_c)²/H²)
- `src/epd/objects.py`: added HopperRegion class (4 walls, exclusion='exterior', rescale support)
- `src/epd/initializer.py`: RSA draws from container bounding box (not full Lx×Ly)
- `src/epd/system.py`: 'poiseuille_v' preset; _accessible_area() uses container polygon;
  set_param('U_max') support; skip phi auto-expand for relax_only=True

**Waypoints:**
- [x] Infrastructure: HopperRegion, Poiseuille-v preset, RSA bounding box, accessible_area
- [x] Quick mode PASS (QUICK=1): no NaN, particles move downward, F1+F2 both PASS
- [ ] Full run (default mode): ≥3 particles exited in each case, GIF movies saved

**Exit gate:**
- QUICK mode: no crash, mean_y decreasing for both F1 and F2  ✅ PASS
- Full mode (pending): n_exited ≥ 3 for F1 and F2  ⬜ PENDING

---

## Phase 6: GitHub Release Package + Reproducibility Notebook

**Objective:** Package v1.0.0 for public upload to GitHub and produce a single
Colab reproducibility notebook that lets any reader run every script from the
methods paper using the v1 API.

**Waypoints:**

- [ ] **6.1** `sys.eval_forces()` — add force-breakdown method to System; returns
      per-particle per-node numpy dict: `{f_elastic, f_contact, f_drag, f_total}`
- [ ] **6.2** Upload folder — create `Upload2Repo/python/` with clean file tree
      (live src/, notebooks/, papers/ minus movies/aux, calibration data, requirements.txt)
- [x] **6.3** `getting_started.ipynb` — prepend clone/install/import header cells;
      add "Logging & Force Analysis" section using two-disk squeeze as worked example
- [x] **6.4a** Port + test `twodisk_squeeze.py` → v1 System/ParticleSpec API ✅
- [x] **6.4b** Port + test `calibration_sweep.py` → v1 API ✅
- [x] **6.4c** Port + test `fixed_b_q_sweep.py` → v1 API ✅
- [x] **6.4d** Port + test `R0_spot_check.py` → v1 API ✅
- [x] **6.4e** Port + test `contact_law_fd.py` + `contact_law_extended.py` → v1 API ✅
- [x] **6.4f** Port + test `emulsion_benchmarks.py` (capwave C, T1 D, falling E) → all v1 API ✅
      Benchmark D: Wall+MotionSpec descent + between-chunk CM repin (no old CapsuleSim)
- [x] **6.4g** Polish `test_E_stokes_drag.py` + `test_E_figures.py` + `test_E_movies.py` ✅
      All confirmed v1 API (ParticleSpec, System, System.run)
- [x] **6.5** Build `papers/summary_of_methods/reproduce_paper.ipynb` ✅
      Valid Jupyter JSON; 4 sections; inline quick demos + full-script run instructions
      Copied to `Upload2Repo/python/notebooks/`; `getting_started.ipynb` also converted to JSON
- [x] **6.6** HANDOFF / LOG / PLAN updates; tag v1.1.0 after upload package verified ✅

**Reproducibility notebook structure:**
- Section 0: Setup (clone, pip install, imports)
- Section 1: Elastic Calibration (§6–7) — fixed-b sweep, full calib sweep, sensitivity, R0 spot-check
- Section 2: Contact Force Law (§8) — F-D curves, N-convergence, power-law fits
- Section 3: Emulsion Droplet Model (§9) — capillary wave, T1 event, falling droplet
- Section 4: Viscous Drag (§10) — E1–E4 benchmarks, figures, movies S12–S13

**Exit gate:**
- `Upload2Repo/python/` contains all needed files, nothing dead/spurious
- `getting_started.ipynb` opens cleanly in Colab from the repo
- `reproduce_paper.ipynb` runs end-to-end in QUICK mode with no errors
- Every ported script has a passing quick-run test before being added to notebook
- All figures match (visually) those in `papers/summary_of_methods/main.pdf`
