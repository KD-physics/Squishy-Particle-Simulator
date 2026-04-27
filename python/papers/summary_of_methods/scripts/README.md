# Scripts Archive — EPD Model Calibration

**Frozen branch:** `frozen/epdm-calibration-v1`
**Commit:** `54aa57a`
**Date archived:** 2026-04-19

These are the exact scripts that were run to produce the data in `../data/`.
They are archived as-is — they require the simulation code from the repo and
must currently be run from the repository root (or with the repo root on
`sys.path`).  If the simulation API changes on a future branch, diff against
`frozen/epdm-calibration-v1` to identify what needs updating.

**Runtime environment:** Python 3.x, `.venv/` in repo root.
All scripts are run as:
```
source .venv/bin/activate
python papers/summary_of_methods/scripts/<script>.py
```
or from the repo root with the script copied there.

**Code dependency:** All scripts import from `src/validation/twodisk_capsule.py`
which in turn depends on `src/simulation/capsule_shell.py` and
`src/simulation/contact_primitives.py`.  These are all frozen at commit
`54aa57a` on branch `frozen/epdm-calibration-v1`.

---

## Scripts (run in this order to fully reproduce)

### 1. `fixed_b_q_sweep.py`
**Purpose:** Initial squishiness-axis exploration at fixed b=0.2, N=32.
**Output:** `results/fixed_b_q_sweep/` (created relative to working directory)
- `fixed_b_results.json` → archive copy in `../data/fixed_b_results.json`
- `fixed_b_sweep.png` → `../figures/fixed_b_sweep.png`
- `perimeter_overlay.png` → `../figures/perimeter_overlay.png`
- `comparison_q_sweep.gif`, `q*.gif` → `../movies/S1*–S6*.gif`

**Runtime:** ~5 minutes (N=32, 10 q values)
**Paper figures:** Fig. 1 (both panels), Supplemental Movies S1–S6

---

### 2. `calibration_sweep.py`
**Purpose:** Full calibration sweep — N × q × ε_ref grid.
**Output:** `results/calibration_sweep/` (relative to working directory)
- `calibration_data.json` → `../data/calibration_data.json`
- `calibration_table.csv` → `../data/calibration_table.csv`
- `calib_N_sensitivity.png` → `../figures/calib_N_sensitivity.png`
- `calib_eps_sensitivity.png` → `../figures/calib_eps_sensitivity.png`
- `calib_dA_nu_trajectory.png` → `../figures/calib_dA_nu_trajectory.png`
- `calib_N_drift.png` → `../figures/calib_N_drift.png`
- `calib_full_grid.png` → `../figures/calib_full_grid.png`

**Runtime:** ~65 minutes (3 N values × 18 q values; N=72 is O(N²) per step)
**Parameters hardcoded in script:** N_VALS=[32,48,72], 18 q values, 5 ε_ref values,
b=0.2, τ=√2.4, S=1, R0=1, C₀=3000, α₀=2.0, SR=0.01, dt_factor=0.4, eps_max=0.14

**Paper figures:** Fig. 2, Fig. 3 (both panels), Fig. 4, Table I, Table II

---

### 3. `R0_spot_check.py`
**Purpose:** Verify R0 invariance at three points on the squishiness axis.
**Depends on:** `results/calibration_sweep/calibration_data.json` (from script 2)
**Output:** `results/calibration_sweep/`
- `R0_spot_check.json` → `../data/R0_spot_check.json`
- `R0_spot_check.png` → `../figures/R0_spot_check.png`

**Runtime:** ~5 minutes (3 ν_targets × 3 R0 values × 1 simulation each)
**Paper figures:** Fig. 5, Table III

---

### 4. `insert_calib_results.py`
**Purpose:** Read `calibration_data.json` and print formatted tables for
insertion into the paper (§7.1 calibration table, N-sensitivity table,
ε_ref-sensitivity table).  Does not write files — prints to stdout.
**Depends on:** `results/calibration_sweep/calibration_data.json`
**Runtime:** <1 second

---

## Notes on path conventions

Scripts write outputs to `results/<name>/` relative to the working directory.
When run from the repo root, this matches the original data location.
The archived copies in `../data/` and `../figures/` were copied manually after
each run.  The `cp` commands used are recorded in the session log (`LOG.md`).

## Tethering to a future code branch

When simulation code is updated (e.g., `capsule_shell.py` API changes):
1. Create a new branch from `master` with the updated code
2. Run these scripts against the new code
3. Compare outputs to `../data/` files — CV should be <0.01% for all observables
4. If the results match, the new branch inherits this calibration
5. If not, diagnose the discrepancy before proceeding
