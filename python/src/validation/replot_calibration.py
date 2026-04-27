"""
replot_calibration.py — Regenerate all calibration plots from existing JSON data.

Reads results/calibration_sweep/calibration_data.json and regenerates all figures.
Safe to run at any time; uses whatever N values and q values are present in the JSON.

Usage:
    cd /root/workspace/projects/polyfem
    source .venv/bin/activate
    python src/validation/replot_calibration.py
"""

import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

OUTDIR = Path("results/calibration_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)

CAL_JSON = OUTDIR / "calibration_data.json"

# ── Load data ──────────────────────────────────────────────────────────────────
data = json.load(open(CAL_JSON))
rows = []
for entry in data:
    for eps_key, mets in entry.get('metrics', {}).items():
        rows.append({
            'N': entry['N'],
            'q': entry['q'],
            'eps_ref': float(eps_key),
            'nu': mets['nu'],
            'dA': mets['dA'],
        })
df = pd.DataFrame(rows)

# Determine which N and q values are present
N_VALS_PRESENT = sorted(df['N'].unique())
Q_VALS_ALL = sorted(df['q'].unique())
EPS_REF_VALS = [0.04, 0.06, 0.08, 0.10, 0.12]

B_TARGET = 0.2
TAU      = np.sqrt(12.0 * B_TARGET)
B_ACTUAL = TAU**2 / 12.0

N_COLORS  = {32: '#1f77b4', 48: '#ff7f0e', 72: '#2ca02c', 120: '#d62728', 240: '#9467bd'}
EPS_STYLES = {0.04: (':', 'x'), 0.06: ('--', 's'), 0.08: ('-', 'o'),
              0.10: ('-.', '^'), 0.12: ((0,(3,1,1,1)), 'D')}

print(f"Loaded {len(df)} entries from {CAL_JSON}")
print(f"N values present: {N_VALS_PRESENT}")
print(f"q values: {len(Q_VALS_ALL)} points")


def get_curve(df, N_val, eps_val, col):
    sub = df[(df['N'] == N_val) & (np.isclose(df['eps_ref'], eps_val))].sort_values('q')
    return sub['q'].values, sub[col].values


# ── Figure 1: N-sensitivity at ε_ref = 0.08 ───────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax_nu, ax_dA = axes

for N in N_VALS_PRESENT:
    q_v, nu_v = get_curve(df, N, 0.08, 'nu')
    q_v, dA_v = get_curve(df, N, 0.08, 'dA')
    if len(q_v) == 0:
        continue
    c = N_COLORS.get(N, 'k')
    ax_nu.semilogx(q_v, nu_v, 'o-', color=c, lw=2, ms=6, label=f'N={N}')
    ax_dA.semilogx(q_v, dA_v, 'o-', color=c, lw=2, ms=6, label=f'N={N}')

for ax in axes:
    ax.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2, label='q·b=1')
    ax.set_xlabel('q  (area stiffness ratio)', fontsize=12)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

ax_nu.set_ylabel('ν_meas  (2D Poisson ratio)', fontsize=12)
ax_nu.set_title(f'Calibration: ν vs q  |  b={B_ACTUAL:.2f}, ε_ref=0.08\n'
                f'N ∈ {N_VALS_PRESENT}', fontsize=10)
ax_dA.set_ylabel('ΔA (%)  (area loss)', fontsize=12)
ax_dA.set_title(f'Calibration: ΔA vs q  |  b={B_ACTUAL:.2f}, ε_ref=0.08\n'
                f'N ∈ {N_VALS_PRESENT}', fontsize=10)

fig.suptitle(f'N-sensitivity  (b={B_ACTUAL:.2f}, τ={TAU:.4f})  —  N ∈ {N_VALS_PRESENT}',
             fontsize=12)
plt.tight_layout()
p = OUTDIR / 'calib_N_sensitivity.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 2: ε_ref sensitivity at highest complete N ─────────────────────────
N_ref = 72  # use N=72 (fully converged, all 18 q values)
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax_nu, ax_dA = axes

for e_ref in EPS_REF_VALS:
    ls, mk = EPS_STYLES[e_ref]
    q_v, nu_v = get_curve(df, N_ref, e_ref, 'nu')
    q_v, dA_v = get_curve(df, N_ref, e_ref, 'dA')
    if len(q_v) == 0:
        continue
    ax_nu.semilogx(q_v, nu_v, marker=mk, linestyle=ls, lw=2, ms=6, label=f'ε_ref={e_ref:.2f}')
    ax_dA.semilogx(q_v, dA_v, marker=mk, linestyle=ls, lw=2, ms=6, label=f'ε_ref={e_ref:.2f}')

for ax in axes:
    ax.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2, label='q·b=1')
    ax.set_xlabel('q  (area stiffness ratio)', fontsize=12)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

ax_nu.set_ylabel('ν_meas', fontsize=12)
ax_nu.set_title(f'ε_ref sensitivity: ν vs q  |  N={N_ref}, b={B_ACTUAL:.2f}', fontsize=10)
ax_dA.set_ylabel('ΔA (%)', fontsize=12)
ax_dA.set_title(f'ε_ref sensitivity: ΔA vs q  |  N={N_ref}, b={B_ACTUAL:.2f}', fontsize=10)

fig.suptitle(f'Sensitivity to ε_ref  (N={N_ref}, b={B_ACTUAL:.2f})', fontsize=12)
plt.tight_layout()
p = OUTDIR / 'calib_eps_sensitivity.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 3: ΔA–ν trajectory ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 6))
for N in N_VALS_PRESENT:
    q_v, nu_v = get_curve(df, N, 0.08, 'nu')
    q_v, dA_v = get_curve(df, N, 0.08, 'dA')
    if len(q_v) == 0:
        continue
    ax.plot(dA_v, nu_v, 'o-', color=N_COLORS.get(N,'k'), lw=2, ms=6, label=f'N={N}')
    idx = np.argmin(np.abs(q_v - 1/B_ACTUAL))
    ax.scatter([dA_v[idx]], [nu_v[idx]], s=150, facecolors='none',
               edgecolors=N_COLORS.get(N,'k'), lw=2.0, zorder=10)

ax.axhline(0, color='lightgray', ls=':'); ax.axhline(1, color='lightgray', ls=':')
ax.set_xlabel('ΔA (%)  at ε_ref=0.08', fontsize=12)
ax.set_ylabel('ν_meas  at ε_ref=0.08', fontsize=12)
ax.set_title(f'ΔA – ν trajectory  b={B_ACTUAL:.2f}, hollow = q·b=1', fontsize=10)
ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
plt.tight_layout()
p = OUTDIR / 'calib_dA_nu_trajectory.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 4: Full grid ────────────────────────────────────────────────────────
N_plot = [n for n in N_VALS_PRESENT if n in [32,48,72,120,240]]
fig, axes = plt.subplots(len(N_plot), len(EPS_REF_VALS),
                         figsize=(4.5 * len(EPS_REF_VALS), 3.8 * len(N_plot)),
                         sharex=True, sharey=True)
if len(N_plot) == 1:
    axes = axes[np.newaxis, :]

for i, N in enumerate(N_plot):
    for j, e_ref in enumerate(EPS_REF_VALS):
        ax = axes[i][j]
        q_v, nu_v = get_curve(df, N, e_ref, 'nu')
        q_v, dA_v = get_curve(df, N, e_ref, 'dA')
        ax2 = ax.twinx()
        if len(q_v):
            ax.semilogx(q_v, nu_v, 'o-', color='tomato',   lw=1.8, ms=5)
            ax2.semilogx(q_v, dA_v, 's-', color='steelblue', lw=1.8, ms=5)
        ax.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.0)
        ax.set_ylim(-0.1, 1.05); ax2.set_ylim(-0.2, 8)
        if i == 0: ax.set_title(f'ε_ref = {e_ref:.2f}', fontsize=10)
        if j == 0: ax.set_ylabel(f'N={N}\nν_meas', fontsize=9, color='tomato')
        if j == len(EPS_REF_VALS)-1: ax2.set_ylabel('ΔA (%)', fontsize=9, color='steelblue')
        if i == len(N_plot)-1: ax.set_xlabel('q', fontsize=9)
        ax.grid(True, alpha=0.2)

fig.suptitle(f'Calibration grid  b={B_ACTUAL:.2f}  |  rows=N, cols=ε_ref  |  dashed=q·b=1',
             fontsize=12, y=1.01)
plt.tight_layout()
p = OUTDIR / 'calib_full_grid.png'
fig.savefig(p, dpi=130, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 5: N-drift quantification ──────────────────────────────────────────
e_ref = 0.08
nu_matrix = np.full((len(N_VALS_PRESENT), len(Q_VALS_ALL)), np.nan)
dA_matrix = np.full((len(N_VALS_PRESENT), len(Q_VALS_ALL)), np.nan)

for i, N in enumerate(N_VALS_PRESENT):
    q_v, nu_v = get_curve(df, N, e_ref, 'nu')
    q_v, dA_v = get_curve(df, N, e_ref, 'dA')
    for j, q in enumerate(Q_VALS_ALL):
        idx = np.where(np.isclose(q_v, q))[0]
        if len(idx):
            nu_matrix[i, j] = nu_v[idx[0]]
            dA_matrix[i, j] = dA_v[idx[0]]

# Drift = max - min across all available N
n_available = np.sum(np.isfinite(nu_matrix), axis=0)
nu_range = np.where(n_available >= 2,
                    np.nanmax(nu_matrix, axis=0) - np.nanmin(nu_matrix, axis=0),
                    np.nan)
dA_range = np.where(n_available >= 2,
                    np.nanmax(dA_matrix, axis=0) - np.nanmin(dA_matrix, axis=0),
                    np.nan)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax_spread, ax_range = axes

q_arr = np.array(Q_VALS_ALL)
valid = np.isfinite(nu_range)
ax_spread.semilogx(q_arr[valid], nu_range[valid] * 100, 'ko-', lw=2, ms=7)
ax_spread.axhline(1.0, color='red', ls='--', lw=1.2, label='1% threshold')
ax_spread.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2)
ax_spread.set_xlabel('q', fontsize=12)
ax_spread.set_ylabel('max(ν) − min(ν) across N  (pp)', fontsize=11)
ax_spread.set_title(f'N-sensitivity of ν calibration  (ε_ref={e_ref})\n'
                    f'N ∈ {N_VALS_PRESENT}', fontsize=10)
ax_spread.legend(); ax_spread.grid(True, alpha=0.3)

valid2 = np.isfinite(dA_range)
ax_range.semilogx(q_arr[valid2], dA_range[valid2], 'bs-', lw=2, ms=7)
ax_range.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2)
ax_range.set_xlabel('q', fontsize=12)
ax_range.set_ylabel('max(ΔA) − min(ΔA) across N  (%)', fontsize=11)
ax_range.set_title(f'N-sensitivity of ΔA calibration  (ε_ref={e_ref})', fontsize=10)
ax_range.grid(True, alpha=0.3)

plt.tight_layout()
p = OUTDIR / 'calib_N_drift.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


print(f"\nAll figures regenerated. Output: {OUTDIR}/")
