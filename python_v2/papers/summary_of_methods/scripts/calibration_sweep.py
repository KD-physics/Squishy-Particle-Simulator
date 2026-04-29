"""
calibration_sweep.py — Canonical squishiness calibration for the Elastic Perimeter Disk Model.

Fixed b = EI/El_t = 0.2 (τ = sqrt(2.4) ≈ 1.5492).
Sweep:
  N     ∈ {32, 48, 72, 120}          perimeter node counts
  q     ∈ 18-point log grid          area stiffness ratio K_area / El_t
  ε_ref ∈ {0.04, 0.06, 0.08, 0.10, 0.12}  reference strains for reporting

Produces:
  calibration_data.json          — raw results indexed by (N, q)
  calib_N_sensitivity.png        — ν(q) at ε_ref=0.08, one line per N
  calib_eps_sensitivity.png      — ν(q) at N=72, one line per ε_ref
  calib_dA_sensitivity.png       — ΔA(q) at ε_ref=0.08, one line per N
  calib_full_grid.png            — 4×5 grid: N rows × ε_ref columns, each panel ν(q)
  calibration_table.csv          — lookup table: (N, q, eps_ref) → (nu, dA)
"""

import sys, pathlib, json, itertools
sys.path.insert(0, str(pathlib.Path(__file__).parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from twodisk_capsule import run_squeeze_once

OUTDIR = pathlib.Path("results/calibration_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Fixed parameters ───────────────────────────────────────────────────────────
S        = 1.0
R0       = 1.0
RHO_F    = 1.0
ALPHA0   = 2.0
SR       = 0.01          # strain rate ratio (wall speed / c_wave)
N_FRAMES = 80
EPS_MAX  = 0.14          # run to 14% so all ε_ref values are reachable
C_0      = 3000.0

B_TARGET = 0.2
TAU      = np.sqrt(12.0 * B_TARGET)    # ≈ 1.5492
B_ACTUAL = TAU**2 / 12.0              # = 0.2 exactly

# ── Sweep axes ─────────────────────────────────────────────────────────────────
N_VALS = [32, 48, 72]

# 18-point log grid from q=0.05 to q=50; q·b=1 at q=5 falls near centre
Q_VALS = [0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75,
          1.00, 1.50, 2.00, 3.00, 5.00, 7.50,
          10.0, 15.0, 20.0, 30.0, 50.0]

EPS_REF_VALS = [0.04, 0.06, 0.08, 0.10, 0.12]

print(f"Elastic Perimeter Disk Model — Squishiness Calibration")
print(f"b = {B_ACTUAL:.4f}  (τ = {TAU:.4f})")
print(f"q·b = 1  at  q = {1/B_ACTUAL:.1f}")
print(f"N values : {N_VALS}")
print(f"q values : {len(Q_VALS)} points, {Q_VALS[0]}–{Q_VALS[-1]}")
print(f"ε_ref    : {EPS_REF_VALS}")
print(f"Total runs: {len(N_VALS)*len(Q_VALS)}")
print()


def make_params(q):
    El_t   = 12.0 * S / TAU**2
    K_area = q * El_t                   # K_area = q × El_t_base (constant, no R0 factor)
    C      = C_0 * S * (1.0 + q)       # τ-independent contact hardness
    EI     = S * R0**3                  # bending stiffness (R0-corrected)
    return dict(El_t=El_t, K_area=K_area, C=C, EI=EI)


def interp_at_eps(eps_arr, val_arr, eps_target):
    eps_arr = np.asarray(eps_arr, dtype=float)
    val_arr = np.asarray(val_arr, dtype=float)
    valid   = np.isfinite(val_arr)
    if not valid.any():
        return float('nan')
    e, v = eps_arr[valid], val_arr[valid]
    if eps_target < e[0] or eps_target > e[-1]:
        return float('nan')
    idx = int(np.clip(np.searchsorted(e, eps_target), 1, len(e) - 1))
    t   = (eps_target - e[idx-1]) / (e[idx] - e[idx-1] + 1e-30)
    return float(v[idx-1] + t * (v[idx] - v[idx-1]))


# ── Main sweep ─────────────────────────────────────────────────────────────────
all_results = []   # list of dicts, one per (N, q)

for N in N_VALS:
    print(f"\n── N = {N} ──────────────────────────────────────")
    print(f"{'q':>6}  {'q·b':>6}  ", end="")
    for e in EPS_REF_VALS:
        print(f"  ν@{e:.2f}  ΔA@{e:.2f}", end="")
    print(f"  cc%")
    print("-" * (16 + 18 * len(EPS_REF_VALS)))

    for q in Q_VALS:
        p = make_params(q)
        try:
            frames, dt, T_wave, v_wall, p1, p2, worst_pen, worst_cc_gap = run_squeeze_once(
                R0=R0, N=N, tau=TAU, S=S, C=p['C'], rho_f=RHO_F,
                strain_rate_ratio=SR, alpha_damp=1.0, alpha0=ALPHA0,
                n_frames=N_FRAMES, delta_max=1.5, dt_factor=0.4,
                verbose=False, K_area=p['K_area'], side_walls=False,
                eps_max=EPS_MAX,
            )
        except Exception as e:
            print(f"{q:6.2f}  CRASHED: {e}")
            continue

        A0     = frames[0]['A0']
        r_c    = p1.r_c
        cc_pct = abs(worst_cc_gap) / r_c * 100.0

        eps_arr = np.array([0.5 * (fr['eps1'] + fr['eps2']) for fr in frames])
        dA_arr  = np.array([(1 - 0.5 * (fr['A1'] + fr['A2']) / A0) * 100 for fr in frames])
        nu_arr  = np.array([0.5 * (fr.get('nu_meas1', np.nan) + fr.get('nu_meas2', np.nan))
                            for fr in frames])

        # Evaluate at each ε_ref
        metrics = {}
        for e_ref in EPS_REF_VALS:
            metrics[e_ref] = dict(
                nu = interp_at_eps(eps_arr, nu_arr, e_ref),
                dA = interp_at_eps(eps_arr, dA_arr, e_ref),
            )

        # Perimeter at max strain (centred)
        x1_c = frames[-1]['x1'] - frames[-1]['x1'].mean(axis=0)

        row = dict(
            N=N, q=q, q_b=q * B_ACTUAL,
            cc_pct=cc_pct,
            eps_arr=eps_arr.tolist(),
            dA_arr=dA_arr.tolist(),
            nu_arr=nu_arr.tolist(),
            x1_c=x1_c.tolist(),
            metrics={str(e): metrics[e] for e in EPS_REF_VALS},
        )
        all_results.append(row)

        print(f"{q:6.2f}  {q*B_ACTUAL:6.3f}  ", end="")
        for e in EPS_REF_VALS:
            nu_v = metrics[e]['nu']; dA_v = metrics[e]['dA']
            print(f"  {nu_v:+.4f}  {dA_v:6.3f}", end="")
        print(f"  {cc_pct:.1f}")

# ── Save raw JSON ──────────────────────────────────────────────────────────────
json_path = OUTDIR / 'calibration_data.json'
json_path.write_text(json.dumps(all_results, indent=2))
print(f"\nSaved: {json_path}")


# ── Build DataFrame for easy querying ─────────────────────────────────────────
rows = []
for r in all_results:
    for e_str, m in r['metrics'].items():
        rows.append(dict(N=r['N'], q=r['q'], q_b=r['q_b'],
                         eps_ref=float(e_str), nu=m['nu'], dA=m['dA'],
                         cc_pct=r['cc_pct']))
df = pd.DataFrame(rows)
csv_path = OUTDIR / 'calibration_table.csv'
df.to_csv(csv_path, index=False)
print(f"Saved: {csv_path}")


# ── Plotting helpers ───────────────────────────────────────────────────────────
Q_ARR = np.array(Q_VALS)
QB_ARR = Q_ARR * B_ACTUAL
N_COLORS  = {32: '#1f77b4', 48: '#ff7f0e', 72: '#2ca02c', 120: '#d62728'}
EPS_STYLES = {0.04: (':', 'x'), 0.06: ('--', 's'), 0.08: ('-', 'o'),
              0.10: ('-.', '^'), 0.12: ((0,(3,1,1,1)), 'D')}


def get_curve(df, N_val, eps_val, col):
    sub = df[(df['N'] == N_val) & (np.isclose(df['eps_ref'], eps_val))].sort_values('q')
    return sub['q'].values, sub[col].values


# ── Figure 1: N-sensitivity at ε_ref = 0.08 ───────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax_nu, ax_dA = axes

for N in N_VALS:
    q_v, nu_v = get_curve(df, N, 0.08, 'nu')
    q_v, dA_v = get_curve(df, N, 0.08, 'dA')
    c = N_COLORS[N]
    ax_nu.semilogx(q_v, nu_v, 'o-', color=c, lw=2, ms=6, label=f'N={N}')
    ax_dA.semilogx(q_v, dA_v, 'o-', color=c, lw=2, ms=6, label=f'N={N}')

for ax in axes:
    ax.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2, label='q·b=1')
    ax.set_xlabel('q  (area stiffness ratio)', fontsize=12)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

ax_nu.set_ylabel('ν_meas  (2D Poisson ratio)', fontsize=12)
ax_nu.set_title(f'Calibration: ν vs q  |  b={B_ACTUAL:.2f}, ε_ref=0.08\nLines = different perimeter resolutions N', fontsize=10)
ax_dA.set_ylabel('ΔA (%)  (area loss)', fontsize=12)
ax_dA.set_title(f'Calibration: ΔA vs q  |  b={B_ACTUAL:.2f}, ε_ref=0.08\nLines = different perimeter resolutions N', fontsize=10)

fig.suptitle(f'N-sensitivity of the squishiness calibration  (b={B_ACTUAL:.2f}, τ={TAU:.4f})',
             fontsize=12)
plt.tight_layout()
p = OUTDIR / 'calib_N_sensitivity.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 2: ε_ref sensitivity at N = 72 ──────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
ax_nu, ax_dA = axes

N_ref = 72
for e_ref in EPS_REF_VALS:
    ls, mk = EPS_STYLES[e_ref]
    q_v, nu_v = get_curve(df, N_ref, e_ref, 'nu')
    q_v, dA_v = get_curve(df, N_ref, e_ref, 'dA')
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

fig.suptitle(f'Sensitivity of calibration to ε_ref choice  (N={N_ref}, b={B_ACTUAL:.2f})',
             fontsize=12)
plt.tight_layout()
p = OUTDIR / 'calib_eps_sensitivity.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 3: ΔA–ν trajectory per N at ε_ref = 0.08 (phase plot) ─────────────
fig, ax = plt.subplots(figsize=(7, 6))
for N in N_VALS:
    q_v, nu_v = get_curve(df, N, 0.08, 'nu')
    q_v, dA_v = get_curve(df, N, 0.08, 'dA')
    ax.plot(dA_v, nu_v, 'o-', color=N_COLORS[N], lw=2, ms=6, label=f'N={N}')
    # annotate q·b=1 point
    idx = np.argmin(np.abs(q_v - 1/B_ACTUAL))
    ax.scatter([dA_v[idx]], [nu_v[idx]], s=150, facecolors='none',
               edgecolors=N_COLORS[N], lw=2.0, zorder=10)

ax.axhline(0, color='lightgray', ls=':'); ax.axhline(1, color='lightgray', ls=':')
ax.set_xlabel('ΔA (%)  at ε_ref=0.08', fontsize=12)
ax.set_ylabel('ν_meas  at ε_ref=0.08', fontsize=12)
ax.set_title(f'ΔA – ν trajectory (squishiness axis)\nb={B_ACTUAL:.2f}, hollow circles = q·b=1 midpoint',
             fontsize=10)
ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
plt.tight_layout()
p = OUTDIR / 'calib_dA_nu_trajectory.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 4: Full grid — ν(q) for all (N, ε_ref) ────────────────────────────
fig, axes = plt.subplots(len(N_VALS), len(EPS_REF_VALS),
                         figsize=(4.5 * len(EPS_REF_VALS), 3.8 * len(N_VALS)),
                         sharex=True, sharey=True)

for i, N in enumerate(N_VALS):
    for j, e_ref in enumerate(EPS_REF_VALS):
        ax = axes[i][j]
        q_v, nu_v = get_curve(df, N, e_ref, 'nu')
        q_v, dA_v = get_curve(df, N, e_ref, 'dA')

        ax2 = ax.twinx()
        ax.semilogx(q_v, nu_v, 'o-', color='tomato',  lw=1.8, ms=5, label='ν')
        ax2.semilogx(q_v, dA_v, 's-', color='steelblue', lw=1.8, ms=5, label='ΔA')
        ax.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.0)
        ax.set_ylim(-0.1, 1.05)
        ax2.set_ylim(-0.2, 8)

        if i == 0:
            ax.set_title(f'ε_ref = {e_ref:.2f}', fontsize=10)
        if j == 0:
            ax.set_ylabel(f'N={N}\nν_meas', fontsize=9, color='tomato')
        if j == len(EPS_REF_VALS) - 1:
            ax2.set_ylabel('ΔA (%)', fontsize=9, color='steelblue')
        if i == len(N_VALS) - 1:
            ax.set_xlabel('q', fontsize=9)

        ax.grid(True, alpha=0.2)

fig.suptitle(f'Calibration grid: ν (red) and ΔA (blue) vs q\n'
             f'b={B_ACTUAL:.2f}, τ={TAU:.4f}  |  rows=N, cols=ε_ref  |  dashed=q·b=1',
             fontsize=12, y=1.01)
plt.tight_layout()
p = OUTDIR / 'calib_full_grid.png'
fig.savefig(p, dpi=130, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Figure 5: N-drift quantification ──────────────────────────────────────────
# At ε_ref=0.08, for each q, show spread of ν across N values
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
ax_spread, ax_range = axes

e_ref = 0.08
nu_matrix = np.full((len(N_VALS), len(Q_VALS)), np.nan)
dA_matrix = np.full((len(N_VALS), len(Q_VALS)), np.nan)

for i, N in enumerate(N_VALS):
    q_v, nu_v = get_curve(df, N, e_ref, 'nu')
    q_v, dA_v = get_curve(df, N, e_ref, 'dA')
    for j, q in enumerate(Q_VALS):
        idx = np.where(np.isclose(q_v, q))[0]
        if len(idx):
            nu_matrix[i, j] = nu_v[idx[0]]
            dA_matrix[i, j] = dA_v[idx[0]]

nu_range = np.nanmax(nu_matrix, axis=0) - np.nanmin(nu_matrix, axis=0)
dA_range = np.nanmax(dA_matrix, axis=0) - np.nanmin(dA_matrix, axis=0)

ax_spread.semilogx(Q_VALS, nu_range * 100, 'ko-', lw=2, ms=7)
ax_spread.axhline(1.0, color='red', ls='--', lw=1.2, label='1% threshold')
ax_spread.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2)
ax_spread.set_xlabel('q', fontsize=12)
ax_spread.set_ylabel('max(ν) − min(ν) across N  (pp)', fontsize=11)
ax_spread.set_title(f'N-sensitivity of ν calibration  (ε_ref={e_ref})\n'
                    f'N ∈ {N_VALS}', fontsize=10)
ax_spread.legend(); ax_spread.grid(True, alpha=0.3)

ax_range.semilogx(Q_VALS, dA_range, 'bs-', lw=2, ms=7)
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


print(f"\nAll done. Files in: {OUTDIR}/")
print("\nFiles:")
for f in sorted(OUTDIR.glob("*.*")):
    print(f"  {f.name:45s}  {f.stat().st_size//1024:4d} KB")
