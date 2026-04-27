# HANDOFF.md — Session State

> **Every new Claude Code session must read this file first.**
> Then read the relevant phase section of PLAN.md.
> Then read the last 5–10 entries of LOG.md.
> Then execute the **Immediate Next Step** listed below.

---

## CURRENT STATE

| Field | Value |
|-------|-------|
| **Current Branch** | `emulsion-dev` |
| **Current Phase** | Phase 6 — GitHub Release Package + Reproducibility Notebook |
| **Phase Status** | ✅ COMPLETE — all waypoints done; ready to commit + tag v1.1.0 |
| **Last Successful Step** | Phase 6.6: all scripts v1 API compliant; reproduce_paper.ipynb built; getting_started.ipynb converted to proper JSON; old-API scripts removed from Upload2Repo |
| **Immediate Next Step** | Commit Phase 6 work on emulsion-dev; merge or tag v1.1.0 |

### Phase 6 Context

**Repo URL:** https://github.com/KD-physics/Squishy-Particle-Simulator/

**Goal:** Clean upload package (`Upload2Repo/python/`) for GitHub + one big Colab
reproducibility notebook (`papers/summary_of_methods/reproduce_paper.ipynb`) that
ports ALL paper scripts to v1 API and lets readers regenerate every figure/table.

**Key design decisions:**
- Every simulation script uses `ParticleSpec`, `System`, `sys.run(N, sample_every, callback)`
- No script calls `CapsuleSim` directly — users verify the exact v1 code
- `sys.eval_forces()` needed first (pre-req for logging example in getting_started)
- Port one script at a time: convert → quick-run test → add to notebook
- Upload folder excludes: movies, LaTeX aux files, __pycache__, dead NN/FEM code

### Phase 4.8 Context

**Why this phase exists:**
- Emulsion `build()` has two wrong R0 exponents: K_area ∝ R0² (should R0⁰), C ∝ R0¹ (should R0⁰).
  Invisible at 5% polydispersity; breaks at size ratios ≥1.2×.
- alpha_damp (∝ R0⁻¹) and xi_drag (∝ R0⁺¹) must be per-particle tensors of shape (P,) for
  correct polydisperse behaviour, and to fix the 10-retraces/batch TF graph issue.
- Elastic CapsuleParticle was already R0-correct (constructor scales El_t, EI, k_c).

**Design:**
- `p.R0` property setter fires `_recompute_params(r0)` on each particle.
- `System.adjust_params_for_size()` reads actual R0 from geometry → fires setter per particle.
- Scalar alpha_damp / xi_drag inputs are broadcast to (P,) tensors; API is backwards-compatible.
- Retrace fix: ensure (P,) tensors have stable shape across `run_fast` calls.

### Phase 3.5 API Reference

```python
# Build wall primitives tensor dict (called once before simulation)
prim_data = make_prim_data([
    (LineSegment(p0, p1, normal), k_pen_mult, vel_or_None),
    (Box(cx, cy, w, h),           k_pen_mult, vel_or_None),
    (Arc(center, radius, convex), k_pen_mult, vel_or_None),
    (Hopper(apex, angle, height), k_pen_mult, vel_or_None),
])

# Run simulation with tf.while_loop (recommended for GPU)
final_state = run_simulation_tf(
    state0, dt, alpha_damp, g, params, n_steps,
    mgr_cpp, skin, prim_data=prim_data, R0_max=max_R0)
```

Key: `k_pen_mult=1.0` matches frozen branch exactly. `vel` is (2,) constant translation velocity.

### Benchmark Results (CPU, P=10 N=60 polydisperse, 500 steps)
| Branch | Wall time | ms/step | Speedup |
|--------|-----------|---------|---------|
| Frozen Python | 9.46 s | 18.93 | 1.0× |
| TF Python loop | 5.56 s | 11.11 | 1.7× |
| TF while_loop | 6.16 s | 12.32 | 1.5× |

Profile: inter_capsule 90.5%, primitives 9.0%, step_rb 0.4%, cpp_candidacy 1.2%.
GPU expected to change TF numbers significantly; while_loop benefit increases on GPU.

### Phase 3 Architecture Summary (COMPLETE)

**Files created / updated:**
- `src/simulation/tf_sim.py` — Full pipeline: `internal_forces_tf`, `step_rb_tf`, `inter_capsule_forces_tf` (polydisperse), `primitive_forces_tf`, `step_full_tf`, `make_prim_data`, `run_simulation_tf` (tf.while_loop)
- `src/simulation/candidacy_manager.py` — Python 3-level candidacy manager
- `src/simulation/contact_forces_py.py` — Python reference contact kernel
- `src/cpp/candidacy_manager.hpp` — C++ CandidacyManager (header-only)
- `src/cpp/bindings.cpp` — pybind11 bindings
- `CMakeLists.txt`, `pyproject.toml` — pip install . works (scikit-build-core)
- `src/validation/phase31_validate.py` through `phase35_validate.py`
- `src/validation/phase35_benchmark.py` — P=10 N=60 polydisperse wall-squeeze benchmark
- `results/candidacy_corpus/candidacy_corpus.json` — 18-entry C++ validation corpus

**Validation results (all float64, all PASS):**
- Phase 3.1: internal forces max_diff=9.33e-15, step_rb max_diff=4.44e-16, 100-step trajectory max|Δx|=8.88e-16
- Phase 3.2: contact forces max_diff=4.55e-13, full step max_diff=4.44e-16, 50-step trajectory max|Δx|=1.78e-15
- Phase 3.3: TF forces max_diff=4.55e-13, full step max_diff=4.44e-16, 500-step trajectory max|Δx|=9.66e-14
- Phase 3.4: C++ vs Python PASS (all rows identical), corpus 18/18 PASS, e2e max|Δx|=0.00
- Phase 3.5: LineSegment/Box/Arc/Hopper/moving-wall/polydisperse ALL PASS; while_loop runner max|Δx|=0.00

**Key bug found:** `tf.cast(python_float, tf.float64)` silently rounds through float32.
Fix: use `tf.constant(value, dtype=dtype)` for all non-trivial literal scalars in TF functions.

### ⛔ AXED: Phases 3 and 4 (FEM/NN approach)

The FEM solver (PolyFEM) and NN-DtN surrogate were explicitly axed. The following
files and directories are dead weight — do not modify or run them:

- `src/nn/` — neural network model, dataset, training
- `src/data_gen/` — FEM training data generation
- `src/simulation/disk_mesh.py`, `fem_elastic.py`, `two_disk_contact.py`,
  `three_disk_contact.py`, `four_disk_contact.py`, `single_disk_contact.py`,
  `polyfem_shear_benchmark.py`
- `src/validation/dtn_linear_floor_analysis.py`
- `data/processed/` — FEM/NN training HDF5 files
- `results/dtn_*.pt` — trained NN model weights

### F(δ) Contact Law Findings (Phase 1K / 1K+)

Script: `src/validation/contact_law_fd.py`  (SR=0.001, ALPHA0=2.0, b=0.2)
Outputs: `results/contact_law_fd/`  — `fd_summary.json`, `fd_curves.png`, `fd_exponent.png`, `fd_N_convergence.png`

**Gates (N=32, 18 ν values):** ν verification PASS (max Δν=0.027 < 0.04); equilibrium PASS (max 4.3% < 10%).

**N-convergence of power-law exponent n (measured from contact_law_extended.py):**
  N=32: n≈{0.33:1.547, 0.56:1.492, 0.83:1.452}
  N=48: n≈{0.33:1.290, 0.56:1.235, 0.83:1.264}
  N=72: n≈{0.33:1.158, 0.56:1.172, 0.83:1.149}
  N=120: n≈{0.33:1.100, 0.56:1.085, 0.83:1.054}
  N=240: n≈{0.33:1.090, 0.56:1.047, 0.83:1.006}  ← ALL COMPLETE
  5-pt fit: n∞≈{0.33:0.966, 0.56:0.954, 0.83:0.923}, a≈{17.3, 16.1, 16.7}

**Physical interpretation:**
- True EPD contact law (N→∞): sub-linear n∞≈0.92–0.97 (mildly ν-dependent)
- At production N=240: n≈1.00–1.09 (nearly linear to slightly super-linear)
- N=32 development tier gives n≈1.5 (artifact of coarse node contact discretization)

**CRITICAL**: Fit on F_dd (disk-disk contact force), NOT F_wall. F_wall inflated by viscous damping.
At SR=0.005 (vs 0.001), n was spuriously measured as 0.20 — damping artifact.

**Calibration N-drift (eps_ref=0.08):**
  Within N∈{32-72}: max Δν=0.010 (q≈0.75)
  Across N∈{32-240}: max Δν=0.046 (q≈0.75-1.0)
  → Production-tier simulations at N=240 MUST use N=240 calibration table (tab:calib240)
  → N=48 calibration table underestimates ν by ~0.04 at N=240

**New files (2026-04-20+):**
- `src/validation/calibration_sweep.py` — extended to N=[32,48,72,120,240]
- `src/validation/contact_law_extended.py` — F(δ) N-convergence with per-N calibration
- `src/validation/contact_law_fd.py` — core F(δ) measurement, SR=0.001, F_dd (not F_wall)
- `src/validation/replot_calibration.py` — standalone re-plot from JSON
- `run_calib_and_contact_law.sh` — pipeline: wait for cal → run contact law
- `results/calibration_sweep/calibration_data.json` — N=32,48,72,120 (18q each) + N=240 (17q)
- `results/contact_law_fd/fd_N_convergence_extended.{json,png}` — pending (N=240 running)
- `papers/summary_of_methods/main.tex` — §7.2 expanded + new §8 (contact law)
- `papers/summary_of_methods/figures/calib_N_{sensitivity,drift}.png` — updated with N=120,240

### New files (2026-04-20)
- `src/simulation/ellipse_particle.py` — EllipseParticle: equal arc-length, CapsuleParticle subclass
- `src/simulation/square_particle.py` — SquareParticle: N nodes, 4 at corners, M=4ρa², I=2Ma²/3
- `src/validation/ellipse_collision_benchmark.py` — 4 configs, all PASS
- `src/validation/collision_kinematics.py` — COR table, oblique disk impulse validation
- `src/validation/spin_conservation.py` — Spin conservation gold-standard (12 cases)
- `src/validation/square_corner_validation.py` — Square corner-wall impulse test (exploratory)
- `src/validation/direct_impulse_test.py` — Gold-standard unit test for step_rb decomposition (ALL PASS)
- `src/validation/contact_force_analysis.py` — Torque audit + friction/geometry analysis
- `src/validation/ellipse_normal_collision.py` — Designed normal-incidence collision test
- `src/visualization/ellipse_collision_movies.py` — 4 ellipse GIFs (EC1–EC4)
- `src/visualization/rb_extension_movies.py` — upgraded: N_SQ=32, gap=0.20, 3-body momentum display
- `results/ellipse_benchmarks/ellipse_conservation.json` — benchmark data
- `results/spin_conservation/spin_results.json` — 12 spin-conservation cases
- `results/direct_impulse/direct_impulse_results.json` — 48 impulse test cases, all pass
- `results/contact_force/contact_force_results.json` — Torque audit + friction ratio data

---

## IMMEDIATE NEXT STEP

**Phase 4.7 — API version of Test A**

### What Phase 4.6 was
Phase 4.6 built the System API scaffolding (objects, particles, motion, initializer,
system, checkpoint) and ran 12 integration tests (37/37 PASS, 272s). It is closed as
**complete but superseded**: the tests were functional but visually incomplete, and the
workflow was not yet aligned to the NumPy→TF→API→Colab pipeline.

### What Phase 4.7 is
Systematic test-by-test validation of the C++/TF backend and System API, building toward
a "Getting Started" Colab notebook. User dictates tests one at a time. Workflow per test:

```
NumPy reference → TF direct (patch tf_sim.py as needed) → API (patch API as needed) → Colab section
```

### Approved so far
- ✅ Test A elastic TF: `src/validation/test_A_elastic_tf.py` — N=32, ~85s, PASS
- ✅ Test A emulsion TF: `src/validation/test_A_emulsion_tf.py` — N=32, k_reg=10, ~80s, PASS
  - Patch applied: `k_reg_forces_tf()` added to `tf_sim.py`; excluded from RB sums
- ✅ Test A elastic API: `src/validation/test_A_elastic_api.py` — machine-precision PASS (max|Dx|=0.00e+00)
- ✅ Test A emulsion API: `src/validation/test_A_emulsion_api.py` — machine-precision PASS (max|Dx|=0.00e+00)
  - Bug fixed: `OscillatingWall.is_time_varying()` returns True; `System._has_moving_objects()` uses it
  - Added: `OscillatingWall`, `initialize_from_particles()`, `make_squeeze_gif()`, `_capsule_outline_polygon()`
- ✅ Colab notebook: `notebooks/getting_started.ipynb` — both API tests as clean Colab sections
- ✅ Test B1 seed: `src/validation/test_B1_seed.py` — Channel + SquareObstacle, 12 particles, seeding PASS
  - Bug fixed: SquareObstacle walls must be submitted as Polygon (not 4 independent LineSegments) for group masking
  - Added: SquareObstacle, exclusion_area(), r_ref_offset, CompositeObject.to_make_prim_list(), CustomObject.region_polygon()
  - Added: accessible_area correction to compute_phi_outer and System.phi_outer

### API layer patches applied (objects.py / system.py)
- `OscillatingWall(Wall)`: bit-identical y_at(t) = y0 + sign*A*sin(omega*t); overrides `resolved(t)`
- `SimulationObject.is_time_varying()`: delegates to `_motion`; overridden in `OscillatingWall` → True
- `System._has_moving_objects()`: `any(obj.is_time_varying() ...)`
- `System.initialize_from_particles(particles, alpha_damp, dt)`: bypasses RSA/swell for exact geometry
- `System.make_squeeze_gif(frames, ...)`: 2-panel animated GIF matching TF test visual style
- `System._capsule_outline_polygon(x, r_c, n_arc=6)`: arc-rounded capsule outline

### Next action
User dictates next test. Workflow: TF direct first (if not already approved), then API, then Colab section.

**Phase 2 is COMPLETE** (all 11 exit gates pass, movies rendered, §9 in paper).

---

## PHASE 2 EMULSION BENCHMARK RESULTS

### Bucket C — Capillary-wave damping calibration  ✅ ALL PASS
Script: `src/validation/emulsion_damp_calib.py`

- ω₂_analytic = √6 = 2.4495 rad/τ₀;  ω₂_measured = 2.5663;  error = 4.8% ✓
- α_crit = 2ω₂ = 5.13;  underdamped at α≤2, overdamped at α≥5 ✓
- Recommended α_damp = 5.1 (near-critical); collisions use α ≈ 2.5 (Q≈1)

### Bucket D — T1-event three-droplet squeeze  ✅ ALL PASS
Script: `src/validation/emulsion_three_droplet.py`
Geometry: d=1.5 R₀ (equilateral, all gaps = 1.0 R₀); N=60; α_damp=5; SR=0.02

NOTE: gap changed from original d=1.2 R₀ (0.4 R₀ gap) to d=1.5 R₀ (1.0 R₀ gap).
Without bending stiffness, extreme compression causes polygon self-intersection.
1.0 R₀ gap requires ~2:1 aspect ratio at passage — manageable without self-intersection.

Gates redesigned for emulsion physics:
- 2D.2: max area loss during squeeze < 25% (K_area=5/κ=0.2 → 20% compression expected)
- 2D.3: x_cm_centre < 3% R₀ (was F-symmetry; contact force sums unreliable in contact equilibrium)

Results: T1 passage at t/τ₀=143.66 ✓; area max=[13.3%, 13.3%, 19.0%] ✓; x_sym=0.16% ✓

### Bucket E — Falling droplet (physical Bo)  ✅ ALL PASS
Script: `src/validation/emulsion_falling_droplet.py`
Bo=0.01 (physical: 250 µm droplet, γ=10 mN/m → Bo=0.015); α_damp=5.0

Gate 2E.1 redesigned: squeeze sanity check instead of F(δ) load/unload
(CM must be pinned during squeeze; F measured directly from gap formula, not nodal sum)

Results: traj_err=6.36e-05 ✓; pen=0 ✓; settled=True ✓; sag_mono=[−0.003, 0.010, 0.018] ✓

**α_damp selection for dense packings (Phase 1J — CORRECTED):**

The observable post-collision ringing is the n=2 **inextensional bending mode** (NOT the
membrane mode). Verified numerically: ω_2_bend measured = 6.0002 rad/s, theory = 6.000 rad/s.

  ω_2_bend = (6/√5) · √(S / (ρ_f · τ · R0²))
  At reference (S=1, τ=0.20, R0=1, ρ=1): ω_2_bend = 6.0 rad/s

  Design formula:  α = ω_2_bend / Q_target = (6/√5)·√(S/(ρ_f·τ)) / (Q_target · R0)

  Common operating points (reference params):
    Q=3  (few cycles):    α = 2.0  ← current default; COR = 0.73
    Q=1  (no ringing):    α = 6.0;  COR = 0.53
    Q=10 (many cycles):   α = 0.6;  COR = 0.87

  **CRITICAL:** COR is STRONGLY sensitive to α because T_contact ≈ 1.57 ≈ 1.5 × T_2_bend.
  The bending mode resonates during contact. Higher α → less ringing but lower COR.
  This is a physical design tradeoff, not a numerical artifact.

  COR data at reference (τ=0.20, R0=1, N=32, S=1, v0=0.5, C=3000):
    α=0.3 → Q=20 → COR=0.894;  α=1.0 → Q=6 → COR=0.825;  α=2.0 → Q=3 → COR=0.729
    α=6.0 → Q=1  → COR=0.530;  α=20  → Q=0.3 → COR=0.289

  T_decay = 2/α for underdamped (Q ≥ 1).
  For overdamped (Q < 1): T_actual ≈ α/ω_2_bend² > 2/α (decays MORE slowly).

  Script: `src/validation/shell_mode_calibration.py`
  Results: `results/shell_mode_calibration/calibration_results.json`

**Note on calibration parameters:** For production runs use:
  - TAU = sqrt(12*B_TARGET) where B_TARGET = 0.2 (the "b=0.2" working point)
  - K_area = q * El_t (passed explicitly to CapsuleParticle)
  - C = 3000 * S * (1+q)
  NOT tau=0.2 raw (that is a different, stiffer shell at b=0.003).

---

## USER-LEVEL KNOBS

The EPD model has **two user-level parameters**. Everything else is derived.

| Knob | Symbol | Range | Meaning |
|------|--------|-------|---------|
| Squishiness | **ν** | 0.18 – 0.94 | Effective Poisson ratio of the disk |
| Damping | **Q** | 0.3 – 20 | Post-collision ringing quality factor; Q=3 (α₀=2.0) is default |

**To use ν:** interpolate q from `results/calibration_sweep/calibration_data.json`
(18-point table, N=32, b=0.2, ε_ref=0.08), then derive all other params:

```python
q = np.exp(np.interp(nu_target, nu_arr, np.log(q_arr)))   # log-interp
TAU    = np.sqrt(12.0 * 0.2)   # ≈ 1.5492  — fixed at b=0.2 working point
K_area = q * (12.0 * S / TAU**2)
C      = 3000.0 * S * (1.0 + q)
alpha_damp = 2.0 / T_wave      # ALPHA0=2.0 keeps ζ≈16% constant
```

q is an internal detail, not a user knob. See memory/project_nu_knob.md for full recipe.

---

## CANONICAL PARAMETER RECIPE — FULLY LOCKED

For any (τ, q) with S=1, R0=1, N=32 (dev) or N=240 (production):

| Parameter | Formula | Notes |
|-----------|---------|-------|
| El_t_base | 12 / τ² | Membrane stiffness at S=R0=1 |
| K_area    | q × El_t_base | Area spring |
| **C**     | **3000 × (1 + q)** | Contact hardness — τ-INDEPENDENT |
| α₀        | 2.0 | Damping; ζ = α₀/(4π) ≈ 16% |
| dt_factor | 0.4 | Auto-stable for all (τ, q) |

**Key insight:** C must be τ-independent. C = C_RATIO × El_t_base (old formula) scales as 1/τ²,
giving C=62 at τ=0.20 vs C=1000 at τ=0.05 with the same k_c/EI ratio — but EI=S is τ-independent,
so contact must also be τ-independent. At τ=0.20 with old formula: cc_pen=22% (unacceptable).
With C = 3000×(1+q): cc_pen < 2% for all tested (τ, q). See results/C_formula_verify/.

---

## VERIFIED PARAMETER LOCKS

### S = 1 (pure force scale, CV=0.00%)
Script: `src/validation/scaling_verification.py`
Sweep S ∈ {0.01, 0.1, 1.0, 10.0}: F/El_t, ΔA, ν all CV=0% at ε_p = 5%, 8%, 10%.

### R0 = 1 (pure length scale, CV=0.00% after EI fix)
Script: `src/validation/R0_N_convergence.py`
Sweep R0 ∈ {0.25, 0.5, 1.0, 2.0, 4.0} at N=32, 64, 128:

| N   | F_norm (all R0) | ΔA@10% (all R0) |
|-----|-----------------|-----------------|
| 32  | 0.00398         | 0.365%          |
| 64  | 0.00432         | 0.378%          |
| 128 | 0.00504         | 0.399%          |

Values shift with N (N-convergence) but have zero R0 spread.

### α₀ = 2.0 (universal damping)
T_wave = 2πR0/sqrt(El_t_base×R0/rho_f/L0) — depends on τ but NOT q or K_area.
So ζ = α₀/(4π) ≈ 16% is constant for all (τ, q) with fixed α₀=2.0. ✓

### dt_factor = 0.4 (auto-stable for all (τ, q))
estimate_dt_max() checks breathing (∝1/sqrt(K_area)), edge wave (∝1/sqrt(El_t_base)),
and contact spring (∝1/sqrt(C/m_node)). With C=3000×(1+q), dt_contact dominates only
at q>5, but edge wave is already limiting at τ=0.05 (no slowdown vs old formula observed).
Runtime: ~19s/run at N=32 across all tested cases.

### C = 3000 × S × (1+q) (contact hardness — NEW formula)
Script: `src/validation/C_formula_verify.py`
Verification at S=1, R0=1, N=32:

| τ    | q    | C      | cc_pen% |
|------|------|--------|---------|
| 0.05 | 0.1  | 3300   | 1.02    |
| 0.05 | 1.0  | 6000   | 1.90    |
| 0.05 | 5.0  | 18000  | 1.06    |
| 0.05 | 10.0 | 33000  | 0.74    |
| 0.10 | 1.0  | 6000   | 0.80    |
| 0.20 | 1.0  | 6000   | 0.48    |

All cc_pen < 2%.  (Wall_pen = 2–4%; acceptable for penalty DEM.)

---

## BUG FIXES APPLIED

1. **EI = S × K_fluid × R0³** (was R0⁴)
   File: `src/simulation/capsule_shell.py` line 136
   At R0=1: EI unchanged. All Phase 1/1B calibrations intact.

2. **C = 3000 × S × (1+q)** (was C_RATIO × El_t_base × (1+q))
   The old formula's El_t_base ∝ 1/τ² made C collapse to 62 at τ=0.20, causing 22% cc_pen.

---

## SQUISHINESS AXIS ANCHORS (τ=0.05, S=1, R0=1)

| q | ΔA@8% | θ_wall° | character |
|---|-------|---------|-----------|
| 0.1 | 0.94% | 6.5° | glass-like |
| 10  | 0.10% | 10.3° | rubber-like |

τ sweep: θ_wall = 9.67° (τ=0.05) → 7.45° (τ=0.10) → 5.15° (τ=0.20)

---

## PHASE 1G RESULTS

| File | Description |
|------|-------------|
| `results/rb_benchmarks/extension_benchmarks.json` | 1G.1 3-body + 1G.2 squeeze data |
| `results/rb_benchmarks/extension_1G1.png` | 3-body: momentum conservation + COR_eff |
| `results/rb_benchmarks/extension_1G2.png` | Squeeze: trajectory + energy panels per config |
| `results/rb_benchmarks/movies/extension_3B{1,2,3}.gif` | 3-body collision movies with Δp display |
| `results/rb_benchmarks/movies/extension_SQ{1,2,3}.gif` | Squeeze movies (N=32, gap=0.20) |
| `results/ellipse_benchmarks/ellipse_conservation.json` | Ellipse 4-config benchmark data |
| `results/ellipse_benchmarks/ellipse_conservation.png` | Conservation plots (linear + angular) |
| `results/ellipse_benchmarks/movies/ellipse_EC{1,2,3,4}.gif` | Ellipse collision movies |

**3-body key findings:**
- Machine-precision: |Δp/p0| < 2e-15, |ΔL/L0| < 1.4e-14 for all 5 configs
- COR_eff from 0.80 (α=2) to 0.924 (glass, α=0)
- Movies: per-frame Δpx/p₀ and Δpy/p₀ visible text box

**Squeeze key findings (N=32, gap=0.20):**
- All 4 configs pass through gap ✓; energy balance ✓
- Rubber (q=10) ~6× more wall work than glass due to larger r_c (r_c ∝ 1/N, N=32 now)

**Ellipse collision findings:**
- EllipseParticle inherits step_rb() unchanged; angular momentum exact for any shape
- |ΔL|/L₀ ≤ 8e-13 for all 4 configs (equal/diff size, same-I/diff-mass, spinning)
- Theory: contact force F ∝ (xq − cp) gives torque ∝ v × v = 0 → exact cancellation

---

## RESULTS FILES

| File | Description |
|------|-------------|
| `results/stage1b/step_A_results.json` | q sweep data |
| `results/stage1b/step_B_results.json` | τ sweep data |
| `results/stage1b/step_C_four_metrics_tau0.050.png` | Four-metric summary |
| `results/stage1b/movies/` | 4 GIFs: glass/rubber × τ=0.05/0.20 |
| `results/scaling_verification/scaling_results.json` | S scaling data |
| `results/R0_convergence/convergence_cache.json` | R0×N sweep data |
| `results/C_formula_verify/verify_results.json` | C formula calibration results |
| `results/C_calibration/C_calibration_results.json` | C grid (C0 vs C1 comparison) |

---

## ENVIRONMENT

| Component | Status |
|-----------|--------|
| Python venv (.venv/) | ✅ Active |
| scikit-fem, h5py, torch, matplotlib, imageio | ✅ Installed |
| rcpgenerator | ✅ Built from source (deps/RCPGenerator) |
| polyfempy | ❌ Build failed (GCC 13) — not needed |
| GPU | ❌ CPU only |

---

*Last updated: 2026-04-25 — Phase 4.8: Test C (Couette cell) Part 1 + Part 2 PASS; save/restore state implemented.*
