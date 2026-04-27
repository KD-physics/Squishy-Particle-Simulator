"""
Plot F(δ) curves with power-law fits for the paper figure.
3-column layout (one per ν probe) × 2 rows (log-log + residuals).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

DATA = "results/contact_law_fd/fd_N_convergence_extended.json"
OUT  = "results/contact_law_fd/fd_fits_paper.png"

N_COLORS = {32:"#1f77b4", 48:"#ff7f0e", 72:"#2ca02c", 120:"#d62728", 240:"#9467bd"}
NU_LABELS = {0.33: r"$\nu \approx 0.33$", 0.56: r"$\nu \approx 0.56$", 0.83: r"$\nu \approx 0.83$"}
NU_PROBES = [0.33, 0.56, 0.83]

def powerlaw(x, A, n):
    return A * x**n

with open(DATA) as f:
    data = json.load(f)

# flatten all records
records = []
for n_val in ["32","48","72","120","240"]:
    for rec in data[n_val]:
        records.append(rec)

fig, axes = plt.subplots(2, 3, figsize=(11, 7),
                         gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08})

for col, nu_p in enumerate(NU_PROBES):
    ax_main = axes[0, col]
    ax_res  = axes[1, col]

    ax_main.set_title(NU_LABELS[nu_p], fontsize=12)

    for N in [32, 48, 72, 120, 240]:
        rec = next((r for r in records if r["N"] == N
                    and abs(r["nu_probe"] - nu_p) < 0.05), None)
        if rec is None:
            continue

        delta = np.array(rec["delta"])
        F     = np.array(rec["F_dd"])
        n_fit = rec["fit_n"]

        # only use data where both delta>0 and F>0 (skip pre-contact zeros)
        mask = (delta > 0) & (F > 0)
        d = delta[mask]
        f = F[mask]

        # fit A using the known n (or refit both)
        try:
            (A, n_refit), _ = curve_fit(powerlaw, d, f, p0=[0.5, n_fit], maxfev=5000)
        except Exception:
            A, n_refit = np.polyfit(np.log(d), np.log(f), 1)[::-1]
            A, n_refit = np.exp(A), n_refit

        lbl = f"$N={N}$, $n={n_refit:.3f}$"
        color = N_COLORS[N]
        ax_main.loglog(d, f, color=color, lw=1.4, alpha=0.8, label=lbl)

        # fit line (extend slightly)
        d_fit = np.geomspace(d.min(), d.max(), 200)
        ax_main.loglog(d_fit, powerlaw(d_fit, A, n_refit),
                       color=color, lw=1.0, ls="--", alpha=0.95)

        # residuals: only in the "clean" regime (d > 20% of max)
        clean = d > 0.20 * d.max()
        f_pred = powerlaw(d, A, n_refit)
        resid  = 100.0 * (f - f_pred) / f_pred
        ax_res.semilogx(d[clean], resid[clean], color=color, lw=1.0, alpha=0.8)

    # reference slope lines: n=1 (linear) and n=1.5 (3D Hertz)
    d_ref = np.array([3e-3, 4e-2])
    for n_ref, ls_ref, lbl_ref in [(1.0, ":", "$n=1$"), (1.5, "-.", "$n=1.5$")]:
        A_ref = 0.35 / (d_ref[-1]**n_ref)
        ax_main.loglog(d_ref, A_ref * d_ref**n_ref, "k", lw=0.9, ls=ls_ref,
                       alpha=0.55, label=lbl_ref if col == 2 else None)

    ax_main.set_ylabel(r"$F_{\rm dd}$" if col == 0 else "")
    ax_main.tick_params(labelbottom=False)
    ax_main.grid(True, which="both", lw=0.3, alpha=0.4)

    ax_res.axhline(0, color="k", lw=0.6, ls="--")
    ax_res.set_xlabel(r"$\delta / R_0$")
    ax_res.set_ylabel("Resid. (%)" if col == 0 else "")
    ax_res.set_ylim(-25, 25)
    ax_res.grid(True, which="both", lw=0.3, alpha=0.4)
    ax_res.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(10))

    if col == 0:
        ax_main.legend(fontsize=7.5, loc="upper left", framealpha=0.8)
    if col == 2:
        ax_main.legend(fontsize=7.5, loc="upper left", framealpha=0.8)

# overall labels
fig.text(0.5, 0.02, r"Overlap $\delta/R_0$", ha="center", fontsize=11)
fig.text(0.04, 0.6, r"Contact force $F_{\rm dd}$", va="center",
         rotation="vertical", fontsize=11)

fig.suptitle(r"$F_{\rm dd}(\delta)$ with power-law fits $A\delta^n$ — dashed lines",
             fontsize=12, y=0.99)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
plt.savefig(OUT, dpi=180, bbox_inches="tight")
print(f"Saved: {OUT}")
plt.close()
