"""
R0_spot_check.py — Verify that the b=0.2 calibration is R0-invariant (v1 API).

For three target ν values spanning the squishiness axis, reads the q value
from calibration_data.json (N=72, ε_ref=0.08), then runs simulations at
R0 ∈ {0.5, 1.0, 2.0} and checks that ν and ΔA match to within tolerance.

Expected result: CV(ν) < 1%, CV(ΔA) < 1% across R0 for each ν_target.
This validates the dimensional analysis claim that R0 is a pure length scale
and the calibration curve transfers without modification.
"""

import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

from src.validation.twodisk_squeeze import run_squeeze_raw, interpolate_at_strain

OUTDIR   = pathlib.Path("results/calibration_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)
CALIB_F  = OUTDIR / 'calibration_data.json'

# ── Fixed parameters ──────────────────────────────────────────────────────────
ALPHA0    = 2.0
SR        = 0.01
N_FRAMES  = 80
EPS_MAX   = 0.12
EPS_REF   = 0.08
C_FACTOR  = 3000.0
B_TARGET  = 0.2
TAU       = np.sqrt(12.0 * B_TARGET)
B_ACTUAL  = TAU**2 / 12.0
N_CALIB   = 72
R0_VALS   = [0.5, 1.0, 2.0]
NU_TARGETS = [0.30, 0.60, 0.85]


# ── Load calibration and build q(ν_target) interpolant ───────────────────────
if not CALIB_F.exists():
    raise FileNotFoundError(f"Run calibration_sweep.py first — {CALIB_F} not found.")

raw = json.loads(CALIB_F.read_text())
calib_rows = [r for r in raw
              if r['N'] == N_CALIB
              and str(EPS_REF) in r['metrics']
              and np.isfinite(r['metrics'][str(EPS_REF)]['nu'])]
calib_rows.sort(key=lambda r: r['q'])

q_cal  = np.array([r['q']                           for r in calib_rows])
nu_cal = np.array([r['metrics'][str(EPS_REF)]['nu'] for r in calib_rows])

nu_to_logq = interp1d(nu_cal, np.log(q_cal), kind='linear', fill_value='extrapolate')

def nu_to_q(nu_target):
    return float(np.exp(nu_to_logq(nu_target)))


print(f"R0 Spot Check — Elastic Perimeter Disk Model (v1 API)")
print(f"b={B_ACTUAL:.4f}, τ={TAU:.4f}, ε_ref={EPS_REF}")
print(f"Calibration source: N={N_CALIB}, {len(calib_rows)} q-points")
print(f"ν_targets: {NU_TARGETS}")
print(f"R0 values: {R0_VALS}")
print()

results = []

for nu_tgt in NU_TARGETS:
    q_tgt = nu_to_q(nu_tgt)
    El_t  = 12.0 / TAU**2
    print(f"ν_target = {nu_tgt:.2f}  →  q = {q_tgt:.4f}")
    print(f"  El_t={El_t:.3f}  K_area={q_tgt*El_t:.3f}  C={C_FACTOR*(1+q_tgt):.1f}")
    print(f"  {'R0':>5}  {'ν_meas':>8}  {'ΔA%':>7}  {'cc%':>5}")

    row_group = []
    for R0 in R0_VALS:
        try:
            frames = run_squeeze_raw(
                q=q_tgt, tau=TAU, N=N_CALIB, R0=R0,
                delta_max=EPS_MAX, n_frames=N_FRAMES,
                C_factor=C_FACTOR, alpha_damp=ALPHA0, strain_rate_ratio=SR,
                verbose=False,
            )
        except Exception as exc:
            print(f"  {R0:5.2f}  CRASHED: {exc}")
            continue

        nu_meas, dA_frac = interpolate_at_strain(frames, EPS_REF)
        dA_meas = float(dA_frac) * 100 if np.isfinite(dA_frac) else float('nan')
        cc_pct  = abs(min(f['cc_gap_min'] for f in frames)) / (2 * np.pi * R0 / N_CALIB) * 100.0

        print(f"  {R0:5.2f}  {nu_meas:8.4f}  {dA_meas:7.4f}  {cc_pct:.1f}")

        row_group.append(dict(
            nu_target=nu_tgt, R0=R0, q=q_tgt,
            nu_meas=float(nu_meas), dA_meas=float(dA_meas), cc_pct=float(cc_pct),
            eps_arr=[f['wall_strain'] for f in frames],
            nu_arr=[f['nu_meas']     for f in frames],
            dA_arr=[f['dA_frac'] * 100 for f in frames],
        ))
        results.append(row_group[-1])

    nu_v = [r['nu_meas'] for r in row_group if np.isfinite(r['nu_meas'])]
    dA_v = [r['dA_meas'] for r in row_group if np.isfinite(r['dA_meas'])]
    if nu_v:
        cv_nu = np.std(nu_v) / (np.mean(nu_v) + 1e-10) * 100
        cv_dA = np.std(dA_v) / (np.mean(dA_v) + 1e-10) * 100
        print(f"  CV(ν) = {cv_nu:.2f}%   CV(ΔA) = {cv_dA:.2f}%")
    print()


# ── Save results ──────────────────────────────────────────────────────────────
out_json = OUTDIR / 'R0_spot_check.json'
out_json.write_text(json.dumps(results, indent=2))
print(f"Saved: {out_json}")


# ── Figure: R0 spot check ─────────────────────────────────────────────────────
n_tgt = len(NU_TARGETS)
fig, axes = plt.subplots(1, n_tgt, figsize=(6 * n_tgt, 5.5), sharey=False)

R0_COLORS = {0.5: '#1f77b4', 1.0: '#2ca02c', 2.0: '#d62728'}
R0_MARKS  = {0.5: 's', 1.0: 'o', 2.0: '^'}

for ax, nu_tgt in zip(axes, NU_TARGETS):
    group = [r for r in results if r['nu_target'] == nu_tgt]
    for r in group:
        c  = R0_COLORS.get(r['R0'], 'k')
        mk = R0_MARKS.get(r['R0'], 'x')
        ax.plot(r['eps_arr'], r['nu_arr'], color=c, marker=mk, markevery=10,
                lw=2.0, ms=7, label=f"R0={r['R0']}")
    ax.axhline(nu_tgt, color='black', ls='--', lw=1.2, label=f'ν_target={nu_tgt}')
    ax.axvline(EPS_REF, color='gray', ls=':', lw=1.0)
    ax.set_xlabel('wall_strain  (ε_wall)', fontsize=11)
    ax.set_ylabel('ν_meas', fontsize=11)
    ax.set_title(f'ν_target = {nu_tgt:.2f}  (q = {nu_to_q(nu_tgt):.3f})\n'
                 f'R0 ∈ {R0_VALS}', fontsize=10)
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    ax.set_ylim(-0.1, 1.05)

fig.suptitle(f'R0 invariance check  |  b={B_ACTUAL:.2f}, N={N_CALIB}, ε_ref={EPS_REF}\n'
             f'Curves for different R0 should overlap — validates R0 as pure length scale',
             fontsize=11)
plt.tight_layout()
p = OUTDIR / 'R0_spot_check.png'
fig.savefig(p, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {p}")


# ── Summary table ─────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("R0 SPOT CHECK SUMMARY")
print(f"{'ν_target':>9}  {'R0':>5}  {'q':>8}  {'ν_meas':>8}  {'ΔA%':>7}  {'cc%':>5}")
print("-" * 55)
for nu_tgt in NU_TARGETS:
    group  = [r for r in results if r['nu_target'] == nu_tgt]
    q_used = nu_to_q(nu_tgt)
    for r in group:
        print(f"  {nu_tgt:7.2f}  {r['R0']:5.2f}  {q_used:8.4f}  "
              f"{r['nu_meas']:8.4f}  {r['dA_meas']:7.4f}  {r['cc_pct']:.1f}")
    nu_v = [r['nu_meas'] for r in group if np.isfinite(r['nu_meas'])]
    dA_v = [r['dA_meas'] for r in group if np.isfinite(r['dA_meas'])]
    if nu_v:
        print(f"  {'CV':>7}                           "
              f"{np.std(nu_v)/np.mean(nu_v)*100:7.2f}%  "
              f"{np.std(dA_v)/np.mean(dA_v)*100:6.2f}%")
    print()

print(f"All done. Files in: {OUTDIR}/")
