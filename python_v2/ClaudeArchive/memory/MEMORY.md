# MEMORY.md — Index of saved memories

## Active (EPD / Emulsion model)
- [project_test_F_hopper.md](project_test_F_hopper.md) — Test F hopper discharge: HopperRegion, Poiseuille-v preset, RSA bbox seeding, quick PASS; full run pending
- [project_1g_benchmarks.md](project_1g_benchmarks.md) — Phase 1G: 3-body machine-precision momentum, T1-event squeeze test, Arc geometry fix (arc_radius=R0+r_c)
- [project_capsule_calibration.md](project_capsule_calibration.md) — Capsule Shell DEM calibration: C=100 works all S, alpha0=2.0, per-particle strain vs wall strain, Python DEM timing
- [project_scaling_laws.md](project_scaling_laws.md) — Proven scaling laws: S=force scale (CV=0%), R0=length scale (CV=0% after EI=S×R0³ fix). Correct C/K_area/EI rules. Two free params: τ and q.
- [project_epdm_calibration.md](project_epdm_calibration.md) — EPD model canonical: b=0.2, q→ν calibration, R0-invariance rules, parameter recipe, docs in docs/model_and_numerics.md
- [project_ellipse_particle.md](project_ellipse_particle.md) — EllipseParticle: equal-arc-length nodes, subclasses CapsuleParticle, 4-config benchmark PASS (|ΔL|/L₀ ≤ 8e-13)
- [project_alpha_damp_calibration.md](project_alpha_damp_calibration.md) — α_damp calibration: bending mode ω_bend=(6/√5)·√(S/ρτR²), strong α-COR tradeoff (T_contact≈T_bend)
- [project_nu_knob.md](project_nu_knob.md) — ν is the single user-level knob; full recipe: ν→q via calibration table, then TAU=√2.4, K_area=q·El_t, C=3000S(1+q), ALPHA0=2.0

## ⛔ AXED (FEM/NN path — do not use)
- [project_dtn_phase4_findings.md](project_dtn_phase4_findings.md) — AXED PATH: FEM DtN NN findings (4.58% best; 4-disk multi-contact floor ~20%; approach abandoned)
- [project_dataset_status.md](project_dataset_status.md) — AXED PATH: FEM training dataset files (data/processed/)
- [project_scaling_invariance.md](project_scaling_invariance.md) — AXED PATH: Scaling law for FEM/NN data augmentation
