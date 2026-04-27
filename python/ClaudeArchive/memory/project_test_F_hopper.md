---
name: Test F Hopper Discharge
description: Emulsion hopper simulation: HopperRegion object, Poiseuille-v preset, RSA bounding box seeding, quick-mode PASS
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
Test F implements a 2D emulsion hopper discharge simulation combining all EPD capabilities.

**Status:** Infrastructure complete, quick mode PASS (2026-04-26). Full run (80k steps) pending.

**Quick mode command:** `QUICK=1 python src/validation/test_F_hopper.py`
**Full mode command:** `python src/validation/test_F_hopper.py`

**Why:** Showcase simulation for the paper; demonstrates gravity + drag + background flow + custom geometry + polydispersity all working together.

**How to apply:** When asked about Test F or hopper simulation, refer to `src/validation/test_F_hopper.py`.

## Key design decisions

- HopperRegion (objects.py): exclusion='exterior', 4 walls (left/right funnel + vertical), open top/bottom
- Poiseuille vertical preset (tf_sim.py preset 5): U_y(x) = -U_max*(1-(x-x_c)²/H²)
- RSA bounding box: draws from container polygon bbox not full Lx×Ly
- relax_only=True: skips phi auto-expansion check in initialize()
- _accessible_area() now uses container polygon area (not Lx×Ly)

## Hopper geometry (full mode)

- W_out=5.0, W_res=12.0, theta=30°, h_res=57, Lx=12, Ly=65
- Total hopper height ≈ 59 R₀; accessible area ≈ 701 R₀²
- N=20 emulsion droplets, κ=0.02, N_nodes=36, Oh=0.15, poly_sigma=5%

## Performance

- ~55 ms/step for 10 particles (emulsion, N_nodes=36)
- ~112 ms/step for 20 particles
- Full run (20 particles, 80k steps): ~90 min
