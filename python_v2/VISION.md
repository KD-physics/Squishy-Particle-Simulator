# DPM — Vision & Phase Roadmap

## Goal
Transform the CPU-developed EPD simulator into a GPU-optimized, publication-ready codebase.
Primary targets: all paper validation scripts running correctly on RTX 3080, then systematic
speedup targeting 10x per experiment locally and 10–20x additional gain via batch scaling on A100.

---

## Phase 0 — Environment Setup
**Goal:** Working DPM conda environment on the local RTX 3080 machine.

- Python 3.10, TensorFlow 2.15, CUDA 12.x
- All Python dependencies installed (numpy, scipy, matplotlib, pandas, imageio, pybind11)
- C++ build tools in place (cmake, gcc)
- C++ extension built and importable (`_candidacy_cpp`)
- GPU detected and verified (`tf.config.list_physical_devices('GPU')`)
- Smoke test passing (basic TF GPU op + one unit test)

**Success criterion:** `import tensorflow as tf; tf.config.list_physical_devices('GPU')` returns
the RTX 3080, and a basic simulation step runs without error.

---

## Phase 1 — Notebook Validation (reproduce_paper.ipynb)
**Goal:** Every cell in `notebooks/reproduce_paper.ipynb` runs correctly on GPU, producing
outputs that match expected physics. This notebook is the full validation suite for the
academic publication — correctness here is the foundation for all subsequent work.

Cells are worked through sequentially. When a cell fails, the relevant class/method is patched
— never the cell itself. New functionality discovered via broken cells gets built into the classes
so it's available to all callers.

**Scope:**
- Section 0: Setup (cells 3–4) — imports, env
- Section 1: Elastic Calibration (cells 6–7) — squeeze demo, calibration sweep commands
- Section 2: Contact Force Law (cells 9–10) — F(δ) demo, power-law fit, full suite commands
- Section 3: Emulsion Droplet Model (cells 12–14) — capillary wave, falling droplet, T1 event
- Section 4: Viscous Drag (cells 16–17) — terminal velocity demo, full drag suite commands

**Success criterion:** All 10 code cells execute without error on RTX 3080; numerical outputs
match expected physics (ν values, contact exponents, ω₂ frequency, terminal velocity).

---

## Phase 2 — Candidacy Optimization
**Goal:** Reduce the size and rebuild cost of the CapCandidates matrix (K×E contact neighbor list).

Attack vectors:
- Tighter skin distance → smaller E → fewer inter-capsule force evaluations per step
- Skin-distance-triggered rebuild: only rebuild when max particle displacement since last rebuild
  exceeds a threshold (rather than every N steps fixed)
- Profile C++ BVH traversal for cache/branch inefficiencies
- Multi-thread the CandidacyManager across particles

**Target:** 2–4x speedup in the contact force computation step.
**Benchmark:** phase35_benchmark.py before and after.

---

## Phase 3 — Integration Upgrade
**Goal:** Replace the current time integrator with a higher-order symplectic scheme
(velocity Verlet) to allow larger dt without energy drift or instability.

Since integration is cheap relative to force evaluation, the gain comes from taking fewer steps
to cover the same physical time — directly multiplying everything else.

**Target:** 2–3x reduction in steps for equivalent simulation time (i.e., 2–3x larger stable dt).
**Benchmark:** Compare energy conservation and stability at increasing dt for current vs. new integrator.

---

## Phase 4 — Batch Experiments (N_exp Dimension)
**Goal:** Add a leading N_exp batch dimension to all state tensors so that multiple independent
experiments run in a single TF graph pass on the GPU.

Architecture:
- All tensors gain shape `(N_exp, ...)` leading dimension
- TF operations batch naturally; force calculations unchanged in structure
- C++ CandidacyManager extended to index by experiment
- Staggered candidacy rebuilds: experiment i rebuilds on step `i % N_exp`, amortizing CPU cost
- Multi-thread C++ manager across experiments (embarrassingly parallel)

**Target:** Near-linear scaling with N_exp on A100 — e.g., 12 experiments for ~2x the cost of 1,
yielding ~6x per-experiment throughput gain. Combined with Phases 2–3: 10–20x per experiment.

**Memory note:** A100 (40–80GB) can comfortably hold 12+ experiments. RTX 3080 (10GB) likely
supports 2–4 experiments for local development.

---

## Combined Target
| Phase | Mechanism | Expected Gain |
|-------|-----------|---------------|
| Phase 2 | Candidacy optimization | 2–4x |
| Phase 3 | Integration upgrade (larger dt) | 2–3x |
| Phase 4 | Batch experiments (N_exp=12, A100) | 10–20x per experiment |
| **Combined** | | **40–240x theoretical; 10x conservative local target** |

The 10x local target (RTX 3080, single experiment) is achievable through Phases 2+3 alone.
Phase 4 unlocks the per-experiment gains at Colab A100 scale.
