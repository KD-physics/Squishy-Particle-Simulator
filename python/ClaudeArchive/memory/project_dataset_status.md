---
name: Dataset files and filtering
description: Current dataset files, sample counts, filtering decisions, and Phase 4B status including 4-disk analysis
type: project
originSessionId: 65866e77-0c8b-4782-bfa5-f6e8cd23633f
---
# Dataset Status (as of 2026-04-22)

## Files
| File | Samples | Status | Notes |
|------|---------|--------|-------|
| `data/processed/dtn_twodisk_aug.h5` | 2400 (2376 non-deg.) | ✅ USE | 10mat × 10δ × 2disks × 12rot |
| `data/processed/dtn_threedisk_aug.h5` | 2880 | ✅ USE | 10mat × 8δ × 3disks × 12rot |
| `data/processed/dtn_fourdisk_aug.h5` | 3840 | ⚠️ DO NOT TRAIN | 4-disk per-ν floor=15-29%; hurts model |
| `data/processed/dtn_23disk_clean.h5` | 5256 | ✅ BEST TRAINING SET | 2-disk + 3-disk; used for best model |
| `data/processed/dtn_combined.h5` | 9120 | ⚠️ DO NOT TRAIN | 2+3+4 disk; 4-disk degrades results |

## Best Trained Model
`results/dtn_fourier_v2_23disk_600ep.pt` — **4.58% test / 4.80% val rel L2**
Trained on `dtn_23disk_clean.h5` (5256 samples), 600 epochs, FourierDtNv2 nu+cw.

## Why 4-Disk Data Hurts
Two simultaneous orthogonal contacts (4-disk) → constrained force distribution (angles always 90°
apart). Per-ν-bucket linear floor = 15-29% (vs 5% for 2-disk). cw cannot distinguish single wide
contact from two perpendicular contacts. Training on 9K combined: 9.59% test (vs 4.58% on 5K).

## Key Data Quality Issues
1. **24 outlier samples filtered** (rubber_jello E=1e3, very small F): `F_linf/E < 1e-5`
2. **dataset.py filter** now uses L-inf (max-abs) norm + gain < 500: `(F_linf/E > 1e-5) & (u_linf*E/F_linf < 500)`
3. **ν=0.45 anomaly**: higher floor (~10%) near P1/P2 transition (pdms, E=3e4)

## Physics Normalization
Per-sample: F_input = F/max|F|, u_output = u·E/max|F|. Conditioning: z-scored log10(E), ν, and cw.
cw = fraction of nodes with |F|/max|F| > 0.01 (contact width, derived, rotation-invariant).

## Contact Width Distribution (5256-sample dtn_23disk_clean.h5)
- cw range: [0.013, 0.133], mean: 0.045, std: 0.029

## Delta Coverage
- 2-disk: [0.002, 0.004, 0.006, 0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.070]
- 3-disk: [0.004, 0.006, 0.010, 0.015, 0.020, 0.030, 0.040, 0.070]
- 4-disk: [0.004, 0.006, 0.010, 0.015, 0.020, 0.030, 0.040, 0.070] — DO NOT USE for DtN training

## Dataset Item Format (DtNDataset.__getitem__)
Returns 5-tuple: (F, logE, nu, u, cw) — F:(N_peri,2), logE:scalar, nu:scalar, u:(N_peri,2), cw:scalar
