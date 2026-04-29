"""
contact_law_extended.py — F(δ) N-convergence with calibrated N=120 and N=240 (v1 API).

Assumes calibration_sweep.py has already been run and calibration_data.json
contains entries for the requested N values.  Run this AFTER calibration_sweep.py.

Outputs
-------
    results/contact_law_fd/fd_N_convergence_extended.png
    results/contact_law_fd/fd_N_convergence_extended.json
    results/contact_law_fd/fd_curves_by_N.png
"""

import sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.validation.twodisk_squeeze import run_squeeze_raw, interpolate_at_strain
from src.validation.contact_law_fd import (
    _fd_from_frames, fit_power_law, plot_exponent_vs_nu,
    EPS_NU_REF, FIT_FRAC, R0, TAU, C_FACTOR, ALPHA0,
    OUTDIR,
)

OUTDIR_EXT = Path("results/contact_law_fd")
OUTDIR_EXT.mkdir(parents=True, exist_ok=True)

N_EXTENDED   = [32, 48, 72, 120, 240]
NU_PROBES    = [0.33, 0.56, 0.83]
CAL_PATH     = Path("results/calibration_sweep/calibration_data.json")
N_FRAMES     = 100
EPS_MAX      = 0.10
SR_BASE      = 0.001   # at N=32; scaled ∝ sqrt(32/N) to keep eq_err ≈ const


def load_calibration_by_N(cal_path=CAL_PATH, eps_ref=EPS_NU_REF):
    data = json.load(open(cal_path))
    cal  = {}
    for entry in data:
        n_val   = entry['N']
        eps_key = str(eps_ref)
        if eps_key not in entry.get('metrics', {}):
            continue
        cal.setdefault(n_val, []).append(
            (entry['q'], entry['metrics'][eps_key]['nu']))
    for n_val in cal:
        cal[n_val].sort()
        q_arr  = np.array([p[0] for p in cal[n_val]])
        nu_arr = np.array([p[1] for p in cal[n_val]])
        cal[n_val] = (q_arr, nu_arr)
    return cal


def run_probe(q, N, nu_probe):
    """Run one quasi-static contact law sim; return (nu_meas, n_exp, eq_err, delta, F_dd)."""
    sr_N  = SR_BASE * (32.0 / N) ** 0.5
    frames = run_squeeze_raw(
        q=q, tau=TAU, N=N, R0=R0, delta_max=EPS_MAX,
        n_frames=N_FRAMES, C_factor=C_FACTOR, alpha_damp=ALPHA0,
        strain_rate_ratio=sr_N, verbose=False)
    res     = _fd_from_frames(frames, R0=R0, N=N)
    nu_meas, _ = interpolate_at_strain(frames, EPS_NU_REF)
    mid     = (res['eps_p'] >= 0.05) & (res['eps_p'] <= 0.10)
    eq_err  = float(np.nanmean(res['equil_err'][mid]) if mid.any()
                    else np.nanmean(res['equil_err'][-5:]))
    A, n_exp = fit_power_law(res['delta'], res['F_dd'])
    return float(nu_meas), n_exp, eq_err, res['delta'].tolist(), res['F_dd'].tolist()


def main():
    print("=" * 65)
    print("Contact law extended N-convergence  (v1 API)")
    print(f"  N ∈ {N_EXTENDED},  ν probes ≈ {NU_PROBES}")
    print("=" * 65)

    cal = load_calibration_by_N()
    available_N = sorted(cal.keys())
    print(f"Calibration available for N ∈ {available_N}")

    missing = [n for n in N_EXTENDED if n not in cal]
    if missing:
        print(f"WARNING: no calibration for N ∈ {missing}; these will be skipped.")

    all_rows = {}

    for N in N_EXTENDED:
        if N not in cal:
            print(f"\n  N={N}: skipped (no calibration data)")
            continue
        q_arr, nu_arr = cal[N]
        sr_N = SR_BASE * (32.0 / N) ** 0.5
        print(f"\n── N={N} (SR={sr_N:.5f}) ────────────────────────")
        rows = []
        for nu_probe in NU_PROBES:
            if nu_probe < nu_arr.min() or nu_probe > nu_arr.max():
                print(f"  ν={nu_probe}: out of range, skipping")
                continue
            q_probe = float(np.exp(np.interp(nu_probe, nu_arr, np.log(q_arr))))
            print(f"  ν_probe={nu_probe:.2f}  q={q_probe:.4f}  ", end="", flush=True)
            nu_meas, n_exp, eq_err, delta, F_dd = run_probe(q_probe, N, nu_probe)
            print(f"ν_meas={nu_meas:.3f}  n={n_exp:.3f}  eq_err={eq_err:.1%}")
            rows.append({
                'N': N, 'nu_probe': nu_probe, 'nu_meas': nu_meas,
                'fit_n': n_exp, 'equil_err': eq_err, 'SR': float(sr_N),
                'delta': delta, 'F_dd': F_dd,
            })
        all_rows[N] = rows

    out_json = OUTDIR_EXT / "fd_N_convergence_extended.json"
    with open(out_json, 'w') as fh:
        json.dump(all_rows, fh, indent=2, default=float)
    print(f"\nSaved: {out_json}")

    print("\n  N    ν_probe   ν_meas    n      eq_err")
    print("  " + "-" * 45)
    for N, rows in all_rows.items():
        for r in rows:
            n_s = f"{r['fit_n']:.3f}" if r['fit_n'] is not None else "  N/A"
            print(f"  {N:3d}  {r['nu_probe']:.2f}      "
                  f"{r['nu_meas']:.3f}   {n_s}   {r['equil_err']:.1%}")

    # 1/N extrapolation
    print("\n=== 1/N extrapolation: n(N) = n_inf + a/N ===")
    n_inf_vals = {}
    for nu_probe in NU_PROBES:
        N_pts, n_pts = [], []
        for N, rows in all_rows.items():
            match = [r for r in rows if abs(r['nu_probe'] - nu_probe) < 0.01 and r['fit_n']]
            if match:
                N_pts.append(N); n_pts.append(match[0]['fit_n'])
        if len(N_pts) >= 2:
            c = np.polyfit(1.0 / np.array(N_pts, float), n_pts, 1)
            n_inf_vals[nu_probe] = c[1]
            print(f"  ν≈{nu_probe:.2f}: n_inf={c[1]:.3f}  a={c[0]:.1f}  "
                  f"n(240)≈{c[1]+c[0]/240:.3f}  n(360)≈{c[1]+c[0]/360:.3f}")

    make_plots(all_rows, n_inf_vals)


def make_plots(all_rows, n_inf_vals):
    probe_colors = ['C0', 'C1', 'C2']
    N_colors = {32: '#1f77b4', 48: '#ff7f0e', 72: '#2ca02c',
                120: '#d62728', 240: '#9467bd'}

    # Plot A: n vs 1/N with extrapolation
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax_ninv, ax_nu = axes

    for j, nu_probe in enumerate(NU_PROBES):
        N_pts, n_pts = [], []
        for N, rows in all_rows.items():
            match = [r for r in rows if abs(r['nu_probe'] - nu_probe) < 0.01 and r['fit_n']]
            if match:
                N_pts.append(N); n_pts.append(match[0]['fit_n'])
        if not N_pts:
            continue
        inv = 1.0 / np.array(N_pts, float)
        ax_ninv.plot(inv, n_pts, 'o-', color=probe_colors[j],
                     lw=1.8, ms=7, label=f"ν≈{nu_probe}")
        for N_, n_ in zip(N_pts, n_pts):
            ax_ninv.annotate(f"N={N_}", (1/N_, n_),
                             textcoords="offset points", xytext=(4, 3),
                             fontsize=7, color=probe_colors[j])
        if len(N_pts) >= 2:
            c = np.polyfit(inv, n_pts, 1)
            x_ext = np.linspace(0, inv.max() * 1.1, 60)
            ax_ninv.plot(x_ext, np.polyval(c, x_ext),
                         '--', color=probe_colors[j], lw=1, alpha=0.6)
            ax_ninv.axhline(c[1], color=probe_colors[j], ls=':', lw=0.8, alpha=0.4,
                            label=f"n∞≈{c[1]:.2f}")

    ax_ninv.axhline(1.0, ls='--', color='gray', lw=1, label='n=1')
    ax_ninv.axhline(0.5, ls=':', color='gray', lw=1, label='n=0.5 (2D Hertz)')
    ax_ninv.set_xlabel("1/N", fontsize=12); ax_ninv.set_ylabel("n", fontsize=12)
    ax_ninv.set_title("n vs 1/N  (dashed = linear extrapolation)", fontsize=10)
    ax_ninv.legend(fontsize=8, loc='lower right')
    ax_ninv.grid(True, alpha=0.3); ax_ninv.set_xlim(left=-0.001)
    ax_ninv.set_ylim(0.4, 1.8)

    for N, rows in all_rows.items():
        nu_p = [r['nu_probe'] for r in rows if r['fit_n']]
        n_p  = [r['fit_n']    for r in rows if r['fit_n']]
        if nu_p:
            ax_nu.plot(nu_p, n_p, 'o-', color=N_colors.get(N, 'k'),
                       lw=1.8, ms=7, label=f"N={N}")
    ax_nu.axhline(1.0, ls='--', color='gray', lw=1)
    ax_nu.axhline(0.5, ls=':', color='gray', lw=1)
    ax_nu.set_xlabel("ν", fontsize=12); ax_nu.set_ylabel("n", fontsize=12)
    ax_nu.set_title("n(ν) for each N", fontsize=10)
    ax_nu.legend(fontsize=9); ax_nu.grid(True, alpha=0.3); ax_nu.set_ylim(0.4, 1.8)

    fig.tight_layout()
    out = OUTDIR_EXT / "fd_N_convergence_extended.png"
    fig.savefig(out, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")

    # Plot B: F(δ) per ν probe, colored by N
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    for j, nu_probe in enumerate(NU_PROBES):
        ax = axes[j]
        for N, rows in all_rows.items():
            match = [r for r in rows if abs(r['nu_probe'] - nu_probe) < 0.01]
            if not match:
                continue
            r     = match[0]
            delta = np.array(r['delta'])
            F     = np.array(r['F_dd'])
            mask  = (F > FIT_FRAC * F.max()) & (delta > 0)
            ax.loglog(delta[mask] / R0, F[mask],
                      color=N_colors.get(N, 'k'), lw=1.5, label=f"N={N}")
        ax.set_xlabel("δ / R₀", fontsize=11); ax.set_ylabel("F", fontsize=11)
        ax.set_title(f"F(δ) log–log,  ν≈{nu_probe}", fontsize=10)
        ax.legend(fontsize=9); ax.grid(True, alpha=0.3, which='both')

    fig.tight_layout()
    out2 = OUTDIR_EXT / "fd_curves_by_N.png"
    fig.savefig(out2, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out2}")


if __name__ == '__main__':
    main()
