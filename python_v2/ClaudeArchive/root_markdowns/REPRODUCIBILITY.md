# REPRODUCIBILITY.md — How to Run All Validations

> For results, metrics, and physics explanation see **VALIDATION.md**.  
> This file is the operator's guide: environment, commands, expected exit codes.

---

## Environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib pandas h5py tqdm meshio scikit-fem triangle scikit-learn
```

Verified with: Python 3.12, scikit-fem 12.0.1, numpy 2.x, scipy 1.12+.

---

## Full Suite (all phases, in order)

Estimated wall time: **~25 minutes** at N=120 (development tier).

```bash
source .venv/bin/activate

python src/validation/phase1_gate.py            # ~1 min  — FEM baseline, analytic check
python src/validation/nu_sweep_convergence.py   # ~1 min  — P1 locking, P2 fix
python src/validation/twodisk_validation.py     # ~10 min — 130 contact configs (V1–V6)
python src/validation/hetero_validation.py      # ~5 min  — heterogeneous pairs (H1–H4)
python src/validation/primitives_validation.py  # ~5 min  — wall, arc, V-notch (B1–B3)
python src/validation/scaling_validation.py     # ~3 min  — scaling invariance + augmentation demo
```

All scripts print `PASS` / `FAIL` per test and exit with:
- `0` if all primary assertions pass  
- (non-zero not enforced; check stdout for `FAIL` lines)

---

## Individual Tests

### Phase 1 — FEM baseline

```bash
# P1: Uniform pressure analytic comparison (all 4 guide materials)
python src/validation/phase1_gate.py

# P2: ν-sweep + N-convergence (P1 locking at ν≥0.45, P2 fix)
python src/validation/nu_sweep_convergence.py

# P3: Perimeter tier consistency (N=120/240/360)
python src/validation/dirichlet_hertz_test.py
```

Expected: all PASS, errors < 0.001%.

### Phase 2 — Two-disk homogeneous

```bash
# Full suite: 10 materials × 13 δ/R = 130 configs, ~10 min
python src/validation/twodisk_validation.py

# Quick: 5 materials × 5 δ/R = 25 configs, ~2 min
python src/validation/twodisk_validation.py --quick

# Production tier (N=240): ~3 hours
python src/validation/twodisk_validation.py --N 240
```

Expected: V1 130/130, V2 machine precision, V3 CV<1%, V4 ratio [0.77–1.24], V5 15/15, V6 exact.

### Phase 2-H — Two-disk heterogeneous

```bash
python src/validation/hetero_validation.py
```

Expected: H1–H4 all PASS. Size ratios 1:1/2:1/4:1, material contrasts glass→rubber.

### Phase 2B — Contact primitives

```bash
# All three tests
python src/validation/primitives_validation.py

# Individual tests
python src/validation/primitives_validation.py B1   # two disks in a box (12 cases)
python src/validation/primitives_validation.py B2   # arc indenter sweep (5 R_arc values)
python src/validation/primitives_validation.py B3   # V-notch corner (no double-counting)
```

Expected: B1 12/12, B2 5/5 + monotonicity, B3 PASS.

### Phase 2C — Swelling initialization

```bash
# All four tests
python src/validation/swelling_validation.py

# Individual tests
python src/validation/swelling_validation.py SW1   # Case A: big swell
python src/validation/swelling_validation.py SW2   # Case B: just-touching start
python src/validation/swelling_validation.py SW3   # rate sensitivity (N_steps 1/5/25/100)
python src/validation/swelling_validation.py SW4   # reversibility (no hysteresis)
```

Expected: ALL PASS. err=0.000% (SW1/SW2), spread=0.000000% (SW3), hysteresis=0.00e+00 (SW4).

### Scaling invariance

```bash
python src/validation/scaling_validation.py
```

Expected: ALL PASS, rel_err = 0.00e+00 (exact floating-point). Augmentation demo printed.

---

## What Each Script Checks

| Script | Key assertion | Tolerance |
|--------|---------------|-----------|
| `phase1_gate.py` | FEM vs analytic u_r | < 0.001% |
| `nu_sweep_convergence.py` | ΔF at N=512→1024 | < 0.3% (P2) |
| `twodisk_validation.py` | Convergence, Newton, E-scaling, Hertz a, large strain, symmetry | Various |
| `hetero_validation.py` | Newton balance heterogeneous; F/(E*·R_eff) consistency | < 0% exact, < 1% |
| `primitives_validation.py` | F_wall ≈ F_dd; a_sim/a_hertz; no corner double-count | < 15%, [0.5,2.0], < 1.5× |
| `scaling_validation.py` | F, u, a scale as λ under X→λX | < 1e-12 (actual: 0) |
| `swelling_validation.py` | swelling ≡ pos-sweep; rate-indep; no hysteresis | 0.000% / 0% / 0 |

---

## Troubleshooting

**SparseEfficiencyWarning from scipy:** harmless — appears when `splu` converts CSR→CSC.
Filter with `2>&1 | grep -v SparseEfficiency` if desired.

**`rubber_jello` fails at large δ/R:** fixed by `k_pen_floor=1500`. If you see F=0, check
that `k_pen = max(0.5·E, 1500)` is active (present in `two_disk_contact.py`).

**`twodisk_validation.py` slow:** use `--quick` flag. Production tier (`--N 240`) takes ~3 hours.

**Bisection in `primitives_validation.py B1` slow:** each case runs 12 bisection steps × 1 inner
FEM solve per step. Normal — takes ~20 s per material/δ case.

---

*Updated: 2026-04-16*
