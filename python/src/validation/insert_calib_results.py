"""
insert_calib_results.py — Read calibration_data.json and print the key tables
for insertion into docs/model_and_numerics.md.

Produces:
  1. N-sensitivity table: ν(q) at ε_ref=0.08 for N ∈ {32,48,72}
  2. ε_ref sensitivity table: ν(q) at N=48 for ε_ref ∈ {0.04,...,0.12}
  3. N-drift summary: max|ν(N) - ν(N=72)| across q, at ε_ref=0.08
  4. ε_ref drift summary: max|ν(ε) - ν(ε=0.08)| across q, at N=48
  5. Updated calibration table (q → ν at N=48, ε_ref=0.08) for §7.1
"""

import json, sys, pathlib
import numpy as np

CALIB_F = pathlib.Path("results/calibration_sweep/calibration_data.json")
if not CALIB_F.exists():
    print("calibration_data.json not found — run calibration_sweep.py first.")
    sys.exit(1)

raw = json.loads(CALIB_F.read_text())

# Build lookup: (N, q, eps_ref) -> {nu, dA}
data = {}
for r in raw:
    for e_str, m in r['metrics'].items():
        key = (r['N'], r['q'], float(e_str))
        data[key] = m

N_vals   = sorted(set(r['N']  for r in raw))
Q_vals   = sorted(set(r['q']  for r in raw))
eps_vals = sorted(set(float(e) for r in raw for e in r['metrics']))

print("=" * 70)
print("CALIBRATION RESULTS SUMMARY")
print("=" * 70)

# ── 1. N-sensitivity at ε_ref=0.08 ────────────────────────────────────────────
e_ref = 0.08
print(f"\n── N-sensitivity  (ε_ref={e_ref}) ──────────────────────────")
header = f"{'q':>7}  {'q·b':>5}"
for N in N_vals:
    header += f"  {'ν(N='+str(N)+')':>10}  {'ΔA':>7}"
print(header)
print("-" * (20 + 20 * len(N_vals)))

B = 0.2
for q in Q_vals:
    line = f"{q:7.3f}  {q*B:5.3f}"
    for N in N_vals:
        m = data.get((N, q, e_ref), {'nu': float('nan'), 'dA': float('nan')})
        line += f"  {m['nu']:+10.4f}  {m['dA']:7.4f}"
    print(line)

# N-drift: max spread across N at each q
nu_by_N = {}
for q in Q_vals:
    vals = [data.get((N, q, e_ref), {}).get('nu', float('nan')) for N in N_vals]
    valid = [v for v in vals if np.isfinite(v)]
    if len(valid) >= 2:
        nu_by_N[q] = max(valid) - min(valid)

max_q = max(nu_by_N, key=nu_by_N.get)
print(f"\nMax N-drift (ε_ref=0.08): Δν = {nu_by_N[max_q]:.4f} at q={max_q}")
print(f"Typical N-drift: Δν = {np.median(list(nu_by_N.values())):.4f} (median across q)")

# ── 2. ε_ref-sensitivity at N=48 ──────────────────────────────────────────────
N_ref = 48
print(f"\n── ε_ref-sensitivity  (N={N_ref}) ─────────────────────────────")
header = f"{'q':>7}  {'q·b':>5}"
for e in eps_vals:
    header += f"  {'ν@'+str(e):>9}"
print(header)
print("-" * (14 + 12 * len(eps_vals)))

for q in Q_vals:
    line = f"{q:7.3f}  {q*B:5.3f}"
    for e in eps_vals:
        m = data.get((N_ref, q, e), {'nu': float('nan')})
        line += f"  {m['nu']:+9.4f}"
    print(line)

# ε_ref drift
drifts = []
for q in Q_vals:
    vals = [data.get((N_ref, q, e), {}).get('nu', float('nan')) for e in eps_vals]
    valid = [v for v in vals if np.isfinite(v)]
    if len(valid) >= 2:
        drifts.append(max(valid) - min(valid))

print(f"\nMax ε_ref drift (N={N_ref}): Δν = {max(drifts):.4f}")
print(f"Typical ε_ref drift: Δν = {np.median(drifts):.4f} (median)")

# ── 3. Final calibration table (N=48, ε_ref=0.08) ─────────────────────────────
print(f"\n── Final calibration table  (N={N_ref}, ε_ref={e_ref}) ────────────")
print(f"{'q':>7}  {'q·b':>6}  {'ν_meas':>8}  {'ΔA%':>7}  character")
print("-" * 55)
chars = {0.05: 'area-compressible', 0.10: 'area-compressible', 0.15: 'area-compressible',
         0.20: 'area-compressible', 0.30: 'compressible', 0.50: 'compressible',
         0.75: 'compressible', 1.0: 'intermediate', 1.5: 'intermediate',
         2.0: 'intermediate', 3.0: 'shape-deformable', 5.0: 'midpoint (q·b=1)',
         7.5: 'shape-deformable', 10.0: 'shape-deformable', 15.0: 'shape-deformable',
         20.0: 'nearly incompressible', 30.0: 'nearly incompressible',
         50.0: 'nearly incompressible'}
for q in Q_vals:
    m = data.get((N_ref, q, e_ref), {'nu': float('nan'), 'dA': float('nan')})
    c = chars.get(q, '')
    print(f"{q:7.3f}  {q*B:6.3f}  {m['nu']:+8.4f}  {m['dA']:7.4f}  {c}")

print("\n" + "=" * 70)
print("Copy the above tables into docs/model_and_numerics.md §7.1 and §7.2")
