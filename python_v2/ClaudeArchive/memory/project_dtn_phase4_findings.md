---
name: Phase 4 DtN architecture findings
description: Empirical results from Phase 4A/4B — linear model floor, Fourier architectures, contact-width conditioning, 4-disk floor analysis
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Phase 4 Architecture Search Results

## Key Finding: 5% Linear Floor (Phase 4A)

After per-sample physics normalization (u_norm = u·E/max|F|, F_norm = F/max|F|), the DtN map is well-approximated as linear: u_norm ≈ C(ν) · F_norm.

Measured linear floor (group-level Fourier fit) by ν bucket on 5256-sample dataset:
- ν ∈ [0.18, 0.30): 5.36%
- ν ∈ [0.30, 0.44): 9.82%  (includes ν=0.45 P1/P2 transition)
- ν ∈ [0.44, 0.50): 7.72%

**All trained Fourier architectures plateau near 5.4–5.9% — near the linear model limit.**

FourierDtNv2 achieves 5.39% with continuous per-sample ν conditioning (better than group-level floor because it uses exact ν per sample).

## Architecture Results (300 epochs, 5256-sample combined dataset)
- FourierDtN (68K): 5.42%
- FourierDtNv2 (160K): **5.39%** (best) — default for all further work
- FourierDtNTrunc k≤40 (26K): 5.91%
- FourierResidual (195K): 5.55% — worse than FourierDtNv2

## Phase 4B Training Results (FINAL)

Full training progression:
- Phase 4A baseline (nu_only, 5256 samples, 300ep): **5.39% test**
- Phase 4B naive (nu_cw, 9096 corrupted samples, 500ep): **9.58% test** ← REGRESSION (corrupted data)
- Phase 4B clean (nu_cw, 5256 2+3-disk, 300ep): **4.83% test**
- Phase 4B clean warm restart (+300ep, LR=2e-4): **4.76% test**
- Phase 4B clean 600ep fresh: **4.58% test / 4.80% val** ← BEST MODEL
- Phase 4B 9K clean (2+3+4-disk, 600ep): **9.59% test** ← 4-disk data HURTS

**Best model**: `results/dtn_fourier_v2_23disk_600ep.pt` — 4.58% test rel L2
Best val epoch: 332 (9K run) / 460 (5K run). Plateau at ~4.56% train / 4.78% val.
**Phase 4B <4% gate NOT achieved.** Model at linear floor for 2+3-disk geometry.

## CRITICAL: 4-Disk Per-ν-Bucket Linear Floor

Even with CLEAN 4-disk data (solver fixed), training is HARMFUL due to fundamental geometry:

**Per-ν-bucket linear floor:**
- 2-disk: 5-6% (good)
- 4-disk: **15-29%** (4-5× higher!)

**Root cause:** Two simultaneous orthogonal contacts create a CONSTRAINED force distribution
(contact angles always 90° apart). The least-squares fit learns a biased A_k matrix that
doesn't generalize. cw conditioning CANNOT distinguish:
- "One wide contact" (e.g., single contact at large δ with cw=10%)  
- "Two perpendicular contacts" (4-disk at δ=0.015 with cw≈10% total)

Adding 4-disk data: clean 5K → 4.58%, 9K (2+3+4-disk) → 9.59% — REGRESSION.

**Lesson**: Multi-contact data requires contact multiplicity/direction conditioning to be useful.
The cw scalar is insufficient to disambiguate contact configurations.

## Phase 4B Code Changes

1. `src/nn/dataset.py`: DtNDataset now returns 5-tuple (F, logE, nu, u, cw).
   cw = z-scored contact fraction, computed on-the-fly.
2. `src/nn/model_spectral.py`: FourierDtNv2 accepts `cw` as optional conditioning input.
   Default `cond_inputs='nu_cw'` (was 'nu_only').
3. `src/nn/train.py`: Updated to unpack 5-tuple batches; _model_forward() passes cw if model accepts it.
4. `data/processed/dtn_23disk_clean.h5`: clean 5256-sample 2+3-disk dataset (NO 4-disk).
5. `src/simulation/four_disk_contact.py`: Fixed warm-start + contact_shrank detection.
   Clean 4-disk data now in `dtn_fourdisk_aug.h5` and `dtn_combined.h5` but should NOT be used for training.

## 4-Disk Solver Fix

The penalty iterative solver was fixed to handle wrong fixed points for stiff materials at large δ:
- Added `u_init` parameter to `run_four_disk_contact` (warm-start)
- Added `contact_shrank` check in `run_four_disk_contact_robust`:
  if n_penetrating drops by >20% from prev_n_pen when stepping up δ → cold-start retry
- Fixed glass_stiff d=0.040: n_pen 11→22-36, F monotone with δ

## Next Steps to Break Below 4%

1. **Accept 4.58% and move to Phase 5** (N-body DEM with NN-DtN) — recommended
2. Phase 4C: Add contact multiplicity conditioning (n_contacts=1,2,...) to FourierDtNv2
3. Phase 4C alternative: Use 2-disk DtN per contact, superpose for multi-contact
