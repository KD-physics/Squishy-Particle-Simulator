# ClaudeArchive — Version 1 Snapshot

This directory is a point-in-time archive of all Claude-authored and Claude-maintained
markdown files from the EPD project, captured at the completion of Phase 6 (emulsion
droplet model, all benchmarks passing, v1 API finalized).

It is intended as a reference for future LLMs or collaborators to understand the full
development context, design decisions, and rationale behind the codebase without needing
to reconstruct it from git history.

---

## Contents

### `root_markdowns/`
Top-level project governance and tracking files:

| File | Purpose |
|------|---------|
| `CLAUDE.md` | Master instructions for Claude Code — scope, autonomy rules, architecture overview |
| `VISION.md` | Science question and design philosophy |
| `PLAN.md` | Full execution roadmap with completed/pending waypoints |
| `HANDOFF.md` | Session handoff state — current phase, last command, next step |
| `LOG.md` | Append-only milestone log of every significant action taken |
| `CAPSULE_SHELL.md` | Technical reference for the CapsuleParticle / CapsuleSim model |
| `VALIDATION.md` | Benchmark results and gate metrics for all validation phases |
| `REPRODUCIBILITY.md` | Instructions to reproduce all paper results from scratch |
| `FUTURE_DIRECTIONS.md` | Performance and scalability roadmap (candidacy optimization, vectorized multi-experiment architecture) |

### `docs/`
Technical documentation written during development:

| File | Purpose |
|------|---------|
| `model_and_numerics.md` | Full model description: forces, parameters, calibration, numerics |
| `rigid_body_integration_benchmark.md` | Rigid-body decomposition integrator design and benchmark results |

### `memory/`
Claude's persistent session memory — key facts, calibration results, and design
decisions accumulated across all sessions:

| File | Purpose |
|------|---------|
| `MEMORY.md` | Index of all memory entries |
| `project_nu_knob.md` | ν→q calibration recipe (the single user-level squishiness knob) |
| `project_epdm_calibration.md` | EPD canonical parameter recipe and calibration results |
| `project_capsule_calibration.md` | Capsule shell calibration: C, alpha, scaling laws |
| `project_scaling_laws.md` | Proven S and R0 scaling invariance rules |
| `project_ellipse_particle.md` | EllipseParticle equal-arc-length node benchmark |
| `project_alpha_damp_calibration.md` | Damping calibration and COR tradeoff |
| `project_1g_benchmarks.md` | Phase 1G: 3-body momentum conservation, T1-event, arc geometry |
| `project_test_F_hopper.md` | Hopper discharge benchmark setup and status |
| `project_dtn_phase4_findings.md` | AXED PATH: FEM/NN DtN findings (kept for history) |
| `project_dataset_status.md` | AXED PATH: FEM training dataset (kept for history) |
| `project_scaling_invariance.md` | AXED PATH: FEM/NN scaling law (kept for history) |

---

## Project Status at Archive Time

- **Phase:** 6 complete — emulsion droplet model, all benchmarks passing
- **API:** v1 finalized (`ParticleSpec`, `System`, `System.run()`)
- **Paper:** EPD methods paper draft complete (`papers/summary_of_methods/`)
- **Axed:** FEM solver and NN-DtN surrogate paths (documented in memory for history)
