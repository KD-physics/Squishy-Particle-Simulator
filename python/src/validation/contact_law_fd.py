"""
contact_law_fd.py — EPD contact law: F(δ) for disk–disk contact (v1 API).

Sweeps the full ν axis using canonical b=0.2 parameters. For each ν:
  1. Looks up q by interpolation from calibration_data.json
  2. Runs two-disk squeeze to ε_wall = 10%
  3. Per frame records:
       δ    = 2*(R0+L0) − d_cm    center-to-centre contact overlap
       F_dd = disk–disk force (y-force on bottom half of top disk)
       F_wall = top-wall force (y-force on top half of top disk)
       ν_meas = lateral/axial strain ratio
  4. Plots F(δ): linear and log–log with power-law fits
  5. Prints summary table: ν_target | ν_meas | power-law exponent | equilibrium error

Force decomposition strategy:
  In quasi-static squeezing the top disk's top nodes contact the wall and
  its bottom nodes contact the lower disk.  We separate by comparing each
  node's y-coordinate to the disk CM.

Usage
-----
    cd /root/workspace/projects/polyfem
    source .venv/bin/activate
    python src/validation/contact_law_fd.py

Outputs
-------
    results/contact_law_fd/fd_curves.png
    results/contact_law_fd/fd_exponent.png
    results/contact_law_fd/fd_summary.json
"""

import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.validation.twodisk_squeeze import run_squeeze_raw, interpolate_at_strain

OUTDIR = Path("results/contact_law_fd")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Canonical parameters ──────────────────────────────────────────────────────
R0         = 1.0
N          = 32
B_TARGET   = 0.2
TAU        = np.sqrt(12.0 * B_TARGET)   # ≈ 1.5492
ALPHA0     = 2.0
SR         = 0.001          # slow → quasi-static → clean F(δ)
EPS_MAX    = 0.10
N_FRAMES   = 100
EPS_NU_REF = 0.08
C_FACTOR   = 3000.0
FIT_FRAC   = 0.02


# ── Load calibration ──────────────────────────────────────────────────────────

def load_calibration(eps_ref=0.08, N_cal=32):
    cal_path = Path("results/calibration_sweep/calibration_data.json")
    if not cal_path.exists():
        raise FileNotFoundError(f"Run calibration_sweep.py first — {cal_path}")
    data = json.load(open(cal_path))
    rows = [(r['q'], r['metrics'][str(eps_ref)]['nu'])
            for r in data
            if r['N'] == N_cal and str(eps_ref) in r.get('metrics', {})]
    rows.sort(key=lambda x: x[0])
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def nu_to_q(nu_target, q_arr, nu_arr):
    return float(np.exp(np.interp(nu_target, nu_arr, np.log(q_arr))))


# ── Force decomposition helpers ───────────────────────────────────────────────

def _fd_from_frames(frames, R0=R0, N=N):
    """
    Extract δ, F_dd, F_wall, nu_meas, eps_wall from frame list.

    Disk indices: 0 = lower disk, 1 = upper disk.
    For the upper disk:
      - nodes with y > CM_y → contact top wall → F_wall
      - nodes with y ≤ CM_y → contact lower disk → F_dd
    """
    L0    = 2.0 * np.pi * R0 / N
    d_ref = 2.0 * (R0 + L0)   # initial CM-to-CM distance (just-touching)

    delta_arr  = []
    F_dd_arr   = []
    F_wall_arr = []
    nu_arr     = []
    eps_arr    = []

    for fr in frames:
        x_all = fr['x_all']   # (2, N, 2)
        x_cm  = fr['x_cm']    # (2, 2)
        fc    = fr['f_contact']  # (2, N, 2)

        # δ = initial CM separation − current CM separation
        d_cm  = x_cm[1, 1] - x_cm[0, 1]
        delta = d_ref - d_cm

        # Split upper disk (index 1) into top/bottom node sets
        cy_top     = x_cm[1, 1]
        bot_nodes  = x_all[1, :, 1] <= cy_top  # → touching lower disk
        top_nodes  = x_all[1, :, 1] >  cy_top  # → touching upper wall

        # F_dd: upward y-force on upper disk from lower disk (positive = repulsive)
        F_dd  = float(fc[1, bot_nodes, 1].sum())
        # F_wall: downward y-force on upper disk from top wall (take abs value)
        F_wall = abs(float(fc[1, top_nodes, 1].sum()))

        delta_arr.append(float(delta))
        F_dd_arr.append(abs(F_dd))
        F_wall_arr.append(F_wall)
        nu_arr.append(fr['nu_meas'])
        eps_arr.append(fr['wall_strain'])

    return {
        'delta':   np.array(delta_arr),
        'F_dd':    np.array(F_dd_arr),
        'F_wall':  np.array(F_wall_arr),
        'nu_meas': np.array(nu_arr),
        'eps_p':   np.array(eps_arr),
        'equil_err': np.abs(np.array(F_dd_arr) - np.array(F_wall_arr))
                     / (np.maximum(np.array(F_wall_arr), 1e-10)),
    }


# ── Power-law fit ─────────────────────────────────────────────────────────────

def fit_power_law(delta, F, fit_frac=FIT_FRAC):
    mask = (F > fit_frac * F.max()) & (delta > 0)
    if mask.sum() < 4:
        return None, None
    n, lA = np.polyfit(np.log(delta[mask]), np.log(F[mask]), 1)
    return float(np.exp(lA)), float(n)


# ── Exponent-vs-ν plot helper ─────────────────────────────────────────────────

def plot_exponent_vs_nu(summary_rows, outpath=None):
    fig, ax = plt.subplots(figsize=(7, 5))
    nu_plot = [r['nu_target'] for r in summary_rows if r['fit_n'] is not None]
    n_plot  = [r['fit_n']     for r in summary_rows if r['fit_n'] is not None]
    ax.plot(nu_plot, n_plot, 'o-', lw=2, ms=6, label=f"N={N}")
    ax.axhline(1.0, ls='--', color='gray', lw=1, label='n = 1 (linear)')
    ax.axhline(0.5, ls=':',  color='gray', lw=1, label='n = 0.5 (2D Hertz)')
    ax.set_xlabel("ν  (target Poisson ratio)", fontsize=12)
    ax.set_ylabel("Power-law exponent  n", fontsize=12)
    ax.set_title(f"F(δ) contact law exponent vs ν  [N={N}]", fontsize=11)
    ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.6)
    if outpath:
        fig.tight_layout()
        fig.savefig(outpath, dpi=130, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {outpath}")
    return ax


# ── N-convergence sweep ───────────────────────────────────────────────────────

N_CONV_VALUES  = [32, 48, 72, 120]
NU_CONV_PROBES = [0.33, 0.56, 0.83]


def run_N_convergence():
    print("\n" + "=" * 65)
    print("N-convergence: F(δ) exponent vs N")
    print(f"  N ∈ {N_CONV_VALUES},  ν probes ≈ {NU_CONV_PROBES}")
    print("=" * 65)

    cal_path = Path("results/calibration_sweep/calibration_data.json")
    cal_data = json.load(open(cal_path))

    cal_by_N = {}
    for entry in cal_data:
        n_val   = entry['N']
        eps_key = str(EPS_NU_REF)
        if eps_key not in entry.get('metrics', {}):
            continue
        cal_by_N.setdefault(n_val, []).append(
            (entry['q'], entry['metrics'][eps_key]['nu']))
    for n_val in cal_by_N:
        cal_by_N[n_val].sort()

    cal_avail = sorted(cal_by_N.keys())
    def get_cal(n_val):
        if n_val in cal_by_N:
            return cal_by_N[n_val], n_val
        nearest = min(cal_avail, key=lambda x: abs(x - n_val))
        print(f"  [N={n_val}: no calibration, using N={nearest} as proxy]")
        return cal_by_N[nearest], nearest

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    ax_exp, ax_lin, ax_log = axes
    ref_colors = plt.cm.tab10(np.linspace(0, 0.5, len(N_CONV_VALUES)))
    all_conv_rows = {}

    for n_idx, n_val in enumerate(N_CONV_VALUES):
        pairs, _ = get_cal(n_val)
        q_arr_n  = np.array([p[0] for p in pairs])
        nu_arr_n = np.array([p[1] for p in pairs])

        rows_n = []
        for nu_probe in NU_CONV_PROBES:
            if nu_probe < nu_arr_n.min() or nu_probe > nu_arr_n.max():
                continue
            q_probe = float(np.exp(np.interp(nu_probe, nu_arr_n, np.log(q_arr_n))))
            print(f"  N={n_val:3d}  ν_probe={nu_probe:.2f}  q={q_probe:.4f}  ",
                  end="", flush=True)

            frames = run_squeeze_raw(
                q=q_probe, tau=TAU, N=n_val, R0=R0, delta_max=EPS_MAX,
                n_frames=N_FRAMES, C_factor=C_FACTOR, alpha_damp=ALPHA0,
                strain_rate_ratio=SR, verbose=False)
            res = _fd_from_frames(frames, R0=R0, N=n_val)

            nu_v, _ = interpolate_at_strain(frames, EPS_NU_REF)
            mid  = (res['eps_p'] >= 0.05) & (res['eps_p'] <= 0.10)
            eq_err = float(np.nanmean(res['equil_err'][mid]) if mid.any()
                           else np.nanmean(res['equil_err'][-5:]))
            A, n_exp = fit_power_law(res['delta'], res['F_dd'])
            print(f"ν_meas={nu_v:.3f}  n={n_exp:.3f}  eq_err={eq_err:.1%}")

            rows_n.append({'N': n_val, 'nu_probe': nu_probe, 'nu_target': nu_probe,
                           'nu_meas': float(nu_v), 'fit_n': n_exp, 'fit_A': A,
                           'equil_err': eq_err,
                           'delta': res['delta'].tolist(),
                           'F_dd': res['F_dd'].tolist()})

        all_conv_rows[n_val] = rows_n
        col = ref_colors[n_idx]
        nu_p = [r['nu_probe'] for r in rows_n]
        n_p  = [r['fit_n']    for r in rows_n]
        ax_exp.plot(nu_p, n_p, 'o-', color=col, lw=1.8, ms=7, label=f"N={n_val}")

        ls_cycle = ['-', '--', ':']
        for j, r in enumerate(rows_n):
            delta = np.array(r['delta']); F = np.array(r['F_dd'])
            ax_lin.plot(delta / R0, F, ls_cycle[j % 3], color=col, lw=1.4)
            mask = (F > FIT_FRAC * F.max()) & (delta > 0)
            if mask.sum() > 3:
                ax_log.loglog(delta[mask] / R0, F[mask], ls_cycle[j % 3],
                              color=col, lw=1.4)

    ax_exp.axhline(1.0, ls='--', color='gray', lw=1, label='n=1')
    ax_exp.axhline(0.5, ls=':', color='dimgray', lw=1, label='n=0.5')
    ax_exp.set_xlabel("ν"); ax_exp.set_ylabel("Exponent n")
    ax_exp.set_title("n vs ν for different N"); ax_exp.legend(fontsize=8)
    ax_exp.grid(True, alpha=0.3)
    ax_lin.set_xlabel("δ/R₀"); ax_lin.set_ylabel("F"); ax_lin.grid(True, alpha=0.3)
    ax_lin.set_title("F(δ) linear")
    ax_log.set_xlabel("δ/R₀"); ax_log.set_ylabel("F"); ax_log.grid(True, alpha=0.3)
    ax_log.set_title("F(δ) log–log")
    for n_idx, n_val in enumerate(N_CONV_VALUES):
        ax_lin.plot([], [], color=ref_colors[n_idx], lw=2, label=f"N={n_val}")
    ax_lin.legend(fontsize=8)
    fig.tight_layout()
    out_conv = OUTDIR / "fd_N_convergence.png"
    fig.savefig(out_conv, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {out_conv}")

    conv_path = OUTDIR / "fd_N_convergence.json"
    with open(conv_path, 'w') as fh:
        json.dump(all_conv_rows, fh, indent=2, default=float)
    print(f"Saved: {conv_path}")


# ── Main sweep ────────────────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("EPD Contact Law: F(δ) sweep over full ν axis  (v1 API)")
    print(f"  b={B_TARGET}  τ={TAU:.4f}  N={N}  ε_max={EPS_MAX*100:.0f}%  "
          f"SR={SR}  n_frames={N_FRAMES}")
    print("=" * 65)

    q_arr, nu_arr = load_calibration(eps_ref=EPS_NU_REF, N_cal=N)
    n_nu = len(q_arr)
    print(f"\nCalibration loaded: {n_nu} (q,ν) pairs, "
          f"ν ∈ [{nu_arr.min():.3f}, {nu_arr.max():.3f}]")

    cmap   = plt.cm.plasma
    colors = cmap(np.linspace(0.05, 0.92, n_nu))

    all_results  = []
    summary_rows = []

    for i, (q_target, nu_cal) in enumerate(zip(q_arr, nu_arr)):
        print(f"\n[{i+1:2d}/{n_nu}]  q={q_target:.4f}  ν_cal={nu_cal:.4f}  ",
              end="", flush=True)

        frames = run_squeeze_raw(
            q=q_target, tau=TAU, N=N, R0=R0, delta_max=EPS_MAX,
            n_frames=N_FRAMES, C_factor=C_FACTOR, alpha_damp=ALPHA0,
            strain_rate_ratio=SR, verbose=False)
        res = _fd_from_frames(frames, R0=R0, N=N)

        nu_final, _ = interpolate_at_strain(frames, EPS_NU_REF)
        mid_mask    = (res['eps_p'] >= 0.05) & (res['eps_p'] <= 0.10)
        eq_err      = float(np.nanmean(res['equil_err'][mid_mask]) if mid_mask.any()
                            else np.nanmean(res['equil_err'][-5:]))
        delta_max   = float(res['delta'].max())
        F_max       = float(res['F_dd'].max())
        A, n_exp    = fit_power_law(res['delta'], res['F_dd'])

        print(f"ν_meas={nu_final:.4f}  F_max={F_max:.3f}  "
              f"δ_max={delta_max:.4f}  exp={n_exp:.3f}  eq_err={eq_err:.2%}")

        all_results.append({
            'q':        float(q_target),
            'nu_target': float(nu_cal),
            'data':     {k: v.tolist() for k, v in res.items()},
            'fit_A':    float(A) if A is not None else None,
            'fit_n':    float(n_exp) if n_exp is not None else None,
        })
        summary_rows.append({
            'q':         q_target,  'nu_target': nu_cal,
            'nu_meas':   float(nu_final) if np.isfinite(nu_final) else None,
            'nu_err':    abs(float(nu_final) - nu_cal) if np.isfinite(nu_final) else None,
            'fit_n':     n_exp,     'fit_A': A,
            'equil_err': eq_err,
        })

    # ── Plot 1: F(δ) linear + log–log ─────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax_lin, ax_log = axes
    norm = plt.Normalize(nu_arr.min(), nu_arr.max())

    for i, (res_dict, row) in enumerate(zip(all_results, summary_rows)):
        data  = res_dict['data']
        delta = np.array(data['delta'])
        F     = np.array(data['F_dd'])
        col   = colors[i]
        ax_lin.plot(delta / R0, F, '-', color=col, lw=1.5)
        mask = (F > FIT_FRAC * F.max()) & (delta > 0)
        if mask.sum() > 3:
            ax_log.loglog(delta[mask] / R0, F[mask], '-', color=col, lw=1.5)
            if res_dict['fit_n'] is not None:
                d_f = np.linspace(delta[mask].min(), delta[mask].max(), 50)
                ax_log.loglog(d_f / R0,
                              res_dict['fit_A'] * d_f ** res_dict['fit_n'],
                              '--', color=col, lw=0.8, alpha=0.5)

    for ax, xsc, ysc, title in [
        (ax_lin, 'linear', 'linear', f"F(δ) — linear  [N={N}, b={B_TARGET}, SR={SR}]"),
        (ax_log, 'log',    'log',    "F(δ) — log–log  (dashed = power-law fit)"),
    ]:
        ax.set_xscale(xsc); ax.set_yscale(ysc)
        ax.set_xlabel("δ / R₀", fontsize=11)
        ax.set_ylabel("F  [model units]", fontsize=11)
        ax.set_title(title, fontsize=10)
        ax.grid(True, alpha=0.3, which='both')

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, ax=axes, label='ν (Poisson ratio)', fraction=0.02, pad=0.01)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fd_curves.png", dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved: {OUTDIR}/fd_curves.png")

    # ── Plot 2: power-law exponent n vs ν ─────────────────────────────────────
    plot_exponent_vs_nu(summary_rows, outpath=OUTDIR / "fd_exponent.png")

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\n" + "=" * 75)
    print(f"{'ν_target':>10}  {'ν_meas':>8}  {'Δν':>8}  "
          f"{'exp n':>7}  {'A':>10}  {'eq_err':>8}")
    print("-" * 75)
    for r in summary_rows:
        n_s = f"{r['fit_n']:.3f}" if r['fit_n'] is not None else "  N/A "
        A_s = f"{r['fit_A']:.3e}" if r['fit_A'] is not None else "    N/A   "
        nu_m = f"{r['nu_meas']:.4f}" if r['nu_meas'] is not None else "  N/A  "
        nu_e = f"{r['nu_err']:.4f}"  if r['nu_err']  is not None else "  N/A  "
        print(f"  {r['nu_target']:>8.4f}  {nu_m:>8}  {nu_e:>8}  "
              f"{n_s:>7}  {A_s:>10}  {r['equil_err']:>7.1%}")

    nu_errs = [r['nu_err']    for r in summary_rows if r['nu_err']    is not None]
    eq_errs = [r['equil_err'] for r in summary_rows]
    n_exps  = [r['fit_n']     for r in summary_rows if r['fit_n']     is not None]
    pass_nu = all(e < 0.04 for e in nu_errs)
    pass_eq = all(e < 0.10 for e in eq_errs)
    print("\n" + "=" * 75)
    print(f"ν verification (|Δν| < 0.04):  "
          f"{'PASS' if pass_nu else 'FAIL'}  (max = {max(nu_errs):.4f})")
    print(f"Equilibrium check (err < 10%):  "
          f"{'PASS' if pass_eq else 'FAIL'}  (max = {max(eq_errs):.1%})")
    if n_exps:
        print(f"Exponents: min={min(n_exps):.3f}  max={max(n_exps):.3f}  "
              f"mean={np.mean(n_exps):.3f}")

    out_path = OUTDIR / "fd_summary.json"
    with open(out_path, 'w') as fh:
        json.dump({'summary': summary_rows,
                   'curves':  all_results,
                   'params':  {'N': N, 'R0': R0, 'TAU': float(TAU),
                               'eps_max': EPS_MAX, 'SR': SR, 'b': B_TARGET}},
                  fh, indent=2, default=float)
    print(f"Saved: {out_path}")
    return pass_nu and pass_eq


if __name__ == '__main__':
    ok = main()
    run_N_convergence()
    sys.exit(0 if ok else 1)
