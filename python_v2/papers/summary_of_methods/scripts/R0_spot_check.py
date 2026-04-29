"""
R0_spot_check.py — Verify that the b=0.2 calibration is R0-invariant.

For three target ν values spanning the squishiness axis, reads the q value
from calibration_data.json (N=72, ε_ref=0.08), then runs simulations at
R0 ∈ {0.5, 1.0, 2.0} and checks that ν and ΔA match to within tolerance.

Expected result: CV(ν) < 1%, CV(ΔA) < 1% across R0 for each ν_target.
This validates the dimensional analysis claim that R0 is a pure length scale
and the calibration curve transfers without modification.
"""

import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

from twodisk_capsule import run_squeeze_once

OUTDIR   = pathlib.Path("results/calibration_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)
CALIB_F  = OUTDIR / 'calibration_data.json'

# ── Fixed parameters ───────────────────────────────────────────────────────────
S        = 1.0
RHO_F    = 1.0
ALPHA0   = 2.0
SR       = 0.01
N_FRAMES = 80
EPS_MAX  = 0.12
EPS_REF  = 0.08
C_0      = 3000.0
B_TARGET = 0.2
TAU      = np.sqrt(12.0 * B_TARGET)
B_ACTUAL = TAU**2 / 12.0
N_CALIB  = 72        # which N to read from calibration
R0_VALS  = [0.5, 1.0, 2.0]
NU_TARGETS = [0.30, 0.60, 0.85]   # three points spanning the squishiness axis


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


# ── Load calibration and build q(ν_target) interpolant ────────────────────────
if not CALIB_F.exists():
    raise FileNotFoundError(f"Run calibration_sweep.py first — {CALIB_F} not found.")

raw = json.loads(CALIB_F.read_text())
calib_rows = [r for r in raw
              if r['N'] == N_CALIB
              and str(EPS_REF) in r['metrics']
              and np.isfinite(r['metrics'][str(EPS_REF)]['nu'])]
calib_rows.sort(key=lambda r: r['q'])

q_cal  = np.array([r['q']                             for r in calib_rows])
nu_cal = np.array([r['metrics'][str(EPS_REF)]['nu']   for r in calib_rows])
dA_cal = np.array([r['metrics'][str(EPS_REF)]['dA']   for r in calib_rows])

# Build monotone interpolant nu → q  (nu is monotone increasing with q)
# Use log(q) for interpolation so we span the range smoothly
nu_to_logq = interp1d(nu_cal, np.log(q_cal), kind='linear', fill_value='extrapolate')

def nu_to_q(nu_target):
    return float(np.exp(nu_to_logq(nu_target)))

def make_params(q, R0_val):
    El_t   = 12.0 * S / TAU**2             # base (R0=1 units), constant
    K_area = q * El_t                       # constant — does NOT scale with R0
    C      = C_0 * S * (1.0 + q)           # constant — does NOT scale with R0
    EI     = S * R0_val**3                  # scales as R0³
    return dict(El_t=El_t, K_area=K_area, C=C, EI=EI)


print(f"R0 Spot Check — Elastic Perimeter Disk Model")
print(f"b={B_ACTUAL:.4f}, τ={TAU:.4f}, ε_ref={EPS_REF}")
print(f"Calibration source: N={N_CALIB}, {len(calib_rows)} q-points")
print(f"ν_targets: {NU_TARGETS}")
print(f"R0 values: {R0_VALS}")
print()

results = []   # list of {nu_target, R0, q, nu_meas, dA_meas}

for nu_tgt in NU_TARGETS:
    q_tgt = nu_to_q(nu_tgt)
    p_R1  = make_params(q_tgt, 1.0)
    print(f"ν_target = {nu_tgt:.2f}  →  q = {q_tgt:.4f}")
    print(f"  El_t={p_R1['El_t']:.3f}  K_area={p_R1['K_area']:.3f}  "
          f"C={p_R1['C']:.1f}  EI(R0=1)={p_R1['EI']:.3f}")
    print(f"  {'R0':>5}  {'EI':>8}  {'ν_meas':>8}  {'ΔA%':>7}  {'cc%':>5}")

    row_group = []
    for R0 in R0_VALS:
        p = make_params(q_tgt, R0)
        try:
            frames, dt, T_wave, v_wall, part1, part2, worst_pen, worst_cc_gap = run_squeeze_once(
                R0=R0, N=N_CALIB, tau=TAU, S=S, C=p['C'], rho_f=RHO_F,
                strain_rate_ratio=SR, alpha_damp=1.0, alpha0=ALPHA0,
                n_frames=N_FRAMES, delta_max=1.5 * R0, dt_factor=0.4,
                verbose=False, K_area=p['K_area'], side_walls=False,
                eps_max=EPS_MAX,
            )
        except Exception as e:
            print(f"  {R0:5.2f}  CRASHED: {e}")
            continue

        A0     = frames[0]['A0']
        r_c    = part1.r_c
        cc_pct = abs(worst_cc_gap) / r_c * 100.0

        eps_arr = np.array([0.5 * (fr['eps1'] + fr['eps2']) for fr in frames])
        dA_arr  = np.array([(1 - 0.5 * (fr['A1'] + fr['A2']) / A0) * 100 for fr in frames])
        nu_arr  = np.array([0.5 * (fr.get('nu_meas1', np.nan) + fr.get('nu_meas2', np.nan))
                            for fr in frames])

        nu_meas = interp_at_eps(eps_arr, nu_arr, EPS_REF)
        dA_meas = interp_at_eps(eps_arr, dA_arr, EPS_REF)

        print(f"  {R0:5.2f}  {p['EI']:8.4f}  {nu_meas:8.4f}  {dA_meas:7.4f}  {cc_pct:.1f}")

        row_group.append(dict(nu_target=nu_tgt, R0=R0, q=q_tgt,
                              nu_meas=nu_meas, dA_meas=dA_meas, cc_pct=cc_pct,
                              eps_arr=eps_arr.tolist(), nu_arr=nu_arr.tolist(),
                              dA_arr=dA_arr.tolist()))
        results.append(row_group[-1])

    # Compute CV across R0 for this ν_target
    nu_vals = [r['nu_meas'] for r in row_group if np.isfinite(r['nu_meas'])]
    dA_vals = [r['dA_meas'] for r in row_group if np.isfinite(r['dA_meas'])]
    if nu_vals:
        cv_nu = np.std(nu_vals) / (np.mean(nu_vals) + 1e-10) * 100
        cv_dA = np.std(dA_vals) / (np.mean(dA_vals) + 1e-10) * 100
        print(f"  CV(ν) = {cv_nu:.2f}%   CV(ΔA) = {cv_dA:.2f}%")
    print()

# ── Save results ───────────────────────────────────────────────────────────────
out_json = OUTDIR / 'R0_spot_check.json'
out_json.write_text(json.dumps(results, indent=2))
print(f"Saved: {out_json}")


# ── Figure: R0 spot check ──────────────────────────────────────────────────────
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
    ax.set_xlabel('ε_p  (engineering strain)', fontsize=11)
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


# ── Summary table ──────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("R0 SPOT CHECK SUMMARY")
print(f"{'ν_target':>9}  {'R0':>5}  {'q':>8}  {'ν_meas':>8}  {'ΔA%':>7}  {'cc%':>5}")
print("-"*55)
for nu_tgt in NU_TARGETS:
    group = [r for r in results if r['nu_target'] == nu_tgt]
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
