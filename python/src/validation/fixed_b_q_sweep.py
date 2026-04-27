"""
fixed_b_q_sweep.py — Fixed b=0.2 squishiness sweep (v1 API).

b = EI/El_t = τ²/12 = 0.2  →  τ = sqrt(2.4) ≈ 1.549

Sweeps q (area stiffness ratio) at fixed b, producing:
  fixed_b_results.json     — sweep data
  fixed_b_sweep.png        — ΔA and ν vs q  (main science figure)
  perimeter_overlay.png    — disk shapes at max strain
  comparison_q_sweep.gif   — side-by-side animation (5 q values)
"""

import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.animation as animation
from pathlib import Path

from src.validation.twodisk_squeeze import run_squeeze_raw, interpolate_at_strain

OUTDIR = pathlib.Path("results/fixed_b_q_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)

N        = 32
SR       = 0.01
N_FRAMES = 80
EPS_MAX  = 0.12
C_FACTOR = 3000.0
ALPHA0   = 2.0
EPS_REF  = 0.08
FPS      = 15

B_TARGET = 0.2
TAU      = np.sqrt(12.0 * B_TARGET)     # ≈ 1.5492
B_ACTUAL = TAU**2 / 12.0               # = 0.2 exactly

Q_VALS   = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]

print(f"Fixed b = {B_ACTUAL:.4f}  (τ = {TAU:.4f})")
print(f"q·b = 1  at  q = {1/B_ACTUAL:.1f}")
print()


# ── Run sweep ─────────────────────────────────────────────────────────────────
print(f"{'q':>6}  {'q·b':>6}  {'ΔA@8%':>8}  {'ν@8%':>8}  {'cc%':>5}")
print("-"*40)

results = []
for q in Q_VALS:
    try:
        frames = run_squeeze_raw(
            q=q, tau=TAU, N=N, delta_max=EPS_MAX, n_frames=N_FRAMES,
            C_factor=C_FACTOR, alpha_damp=ALPHA0, strain_rate_ratio=SR,
            verbose=False,
        )
    except Exception as exc:
        print(f"{q:6.2f}  CRASHED: {exc}")
        continue

    nu_ref, dA_ref = interpolate_at_strain(frames, EPS_REF)

    # Worst capsule-capsule gap as fraction of r_c
    L0   = 2.0 * np.pi / N        # R0=1
    r_c  = L0
    cc_pct = abs(min(f['cc_gap_min'] for f in frames)) / r_c * 100.0

    # Centred perimeter of disk 0 at max strain (for shape overlay)
    x_last = frames[-1]['x_all'][0]
    x0_c   = x_last - x_last.mean(axis=0)

    print(f"{q:6.2f}  {q*B_ACTUAL:6.3f}  {dA_ref*100:8.4f}  {nu_ref:8.4f}  {cc_pct:5.1f}")

    results.append(dict(
        q=q, q_b=q * B_ACTUAL,
        nu_ref=float(nu_ref), dA_ref=float(dA_ref),
        cc_pct=float(cc_pct),
        eps_arr=[f['wall_strain'] for f in frames],
        dA_arr=[f['dA_frac']     for f in frames],
        nu_arr=[f['nu_meas']     for f in frames],
        x0_c=x0_c.tolist(),
        frames=frames,
    ))

# Save JSON (without bulky frame arrays)
json_rows = [{k: v for k, v in r.items() if k != 'frames'} for r in results]
(OUTDIR / 'fixed_b_results.json').write_text(json.dumps(json_rows, indent=2))
print(f"\nSaved: {OUTDIR}/fixed_b_results.json")


# ── Figure 1: ΔA and ν vs q ───────────────────────────────────────────────────
q_vals  = np.array([r['q']      for r in results])
dA_vals = np.array([r['dA_ref'] for r in results]) * 100   # → percent
nu_vals = np.array([r['nu_ref'] for r in results])
qb_vals = q_vals * B_ACTUAL

norm_q = plt.Normalize(vmin=np.log10(q_vals.min()), vmax=np.log10(q_vals.max()))

fig1, axes1 = plt.subplots(1, 3, figsize=(16, 5))
ax_q, ax_qb, ax_phase = axes1

ax2 = ax_q.twinx()
l1, = ax_q.semilogx(q_vals, dA_vals, 'b^-', ms=9, lw=2, label='ΔA (%)')
l2, = ax2.semilogx(q_vals, nu_vals,  'ro-', ms=9, lw=2, label='ν_meas')
ax_q.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2, label=f'q·b=1')
ax_q.set_xlabel('q', fontsize=12)
ax_q.set_ylabel('ΔA (%) at ε=0.08', color='blue', fontsize=11)
ax2.set_ylabel('ν at ε=0.08', color='red', fontsize=11)
ax_q.set_title(f'ΔA and ν vs q\n(b={B_ACTUAL:.2f}, τ={TAU:.3f})', fontsize=10)
ax_q.legend(handles=[l1, l2], fontsize=9); ax_q.grid(True, alpha=0.3)

ax2b = ax_qb.twinx()
l3, = ax_qb.semilogx(qb_vals, dA_vals, 'b^-', ms=9, lw=2, label='ΔA (%)')
l4, = ax2b.semilogx(qb_vals, nu_vals,  'ro-', ms=9, lw=2, label='ν_meas')
ax_qb.axvline(1.0, color='gray', ls='--', lw=1.5, label='q·b = 1')
ax_qb.set_xlabel('q · b', fontsize=12)
ax_qb.set_ylabel('ΔA (%)', color='blue', fontsize=11)
ax2b.set_ylabel('ν', color='red', fontsize=11)
ax_qb.set_title('Same axes rescaled to q·b\n(natural midpoint at q·b=1)', fontsize=10)
ax_qb.legend(handles=[l3, l4], fontsize=9); ax_qb.grid(True, alpha=0.3)

sc = ax_phase.scatter(dA_vals, nu_vals,
                      c=np.log10(q_vals), cmap='RdYlGn',
                      s=120, zorder=5, edgecolors='k', linewidths=0.5)
ax_phase.plot(dA_vals, nu_vals, 'k-', lw=1.0, alpha=0.4, zorder=3)
plt.colorbar(sc, ax=ax_phase, label='log₁₀(q)')
for r in results:
    if r['q'] in [0.1, 1.0, 5.0, 20.0]:
        ax_phase.annotate(f"q={r['q']}", (r['dA_ref'] * 100, r['nu_ref']),
                          textcoords='offset points', xytext=(6, 4), fontsize=8)
qb1_pts = [r for r in results if abs(r['q'] * B_ACTUAL - 1.0) < 0.3]
if qb1_pts:
    r_ = min(qb1_pts, key=lambda r: abs(r['q'] * B_ACTUAL - 1.0))
    ax_phase.scatter([r_['dA_ref'] * 100], [r_['nu_ref']], s=250,
                     facecolors='none', edgecolors='red', lw=2.5, zorder=10, label='q·b≈1')
ax_phase.axhline(0, color='lightgray', ls=':', lw=0.8)
ax_phase.axhline(1, color='gray',      ls=':', lw=0.8, label='ν=1')
ax_phase.set_xlabel('ΔA (%) at ε=0.08', fontsize=12)
ax_phase.set_ylabel('ν_meas at ε=0.08', fontsize=12)
ax_phase.set_title('ΔA – ν trajectory  (b=0.2 fixed, q varies)', fontsize=10)
ax_phase.legend(fontsize=8); ax_phase.grid(True, alpha=0.3)

fig1.suptitle(f'Fixed b = {B_ACTUAL:.2f}  (τ = {TAU:.3f})  |  q sweeps the squishiness axis\n'
              f'q·b=1 at q={1/B_ACTUAL:.0f}  |  left=area-squishy, right=shape-squishy',
              fontsize=11)
plt.tight_layout()
p1 = OUTDIR / 'fixed_b_sweep.png'
fig1.savefig(p1, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {p1}')


# ── Figure 2: Perimeter overlay at max strain ─────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(17, 5.5))
ax_full, ax_top, ax_eq = axes2

cmap_q = plt.cm.RdYlGn

def ref_circle(ax):
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(th), np.sin(th), ':', color='#cccccc', lw=1.0, zorder=0)

for ax in axes2:
    ref_circle(ax); ax.set_aspect('equal'); ax.grid(True, alpha=0.15)
    ax.set_xlabel('x / R₀', fontsize=10)

for r in results:
    c   = cmap_q(norm_q(np.log10(r['q'])))
    lbl = f"q={r['q']}  ΔA={r['dA_ref']*100:.2f}%  ν={r['nu_ref']:.2f}"
    x   = np.array(r['x0_c'])
    pts = np.vstack([x, x[:1]])
    ax_full.plot(pts[:, 0], pts[:, 1], color=c, lw=2.0, label=lbl)

    mask_top = x[:, 1] > 0.3
    ix = np.where(mask_top)[0]
    if len(ix) > 1:
        ax_top.plot(x[ix, 0], x[ix, 1], color=c, lw=2.2)

    mask_eq = np.abs(x[:, 1]) < 0.25
    ix = np.where(mask_eq)[0]
    if len(ix) > 1:
        ax_eq.plot(x[ix, 0], x[ix, 1], color=c, lw=2.2, marker='o', ms=3)

ax_full.set_xlim(-1.3, 1.3); ax_full.set_ylim(-1.3, 1.3)
ax_full.set_title(f'Full perimeter at max strain\n(b={B_ACTUAL:.2f}, green=small q → red=large q)', fontsize=9)
ax_full.legend(fontsize=6, loc='lower right')
ax_full.set_ylabel('y / R₀', fontsize=10)

ax_top.set_xlim(-1.2, 1.2); ax_top.set_ylim(0.3, 1.25)
ax_top.set_title('Wall-contact side (top half)', fontsize=9)

ax_eq.set_xlim(0.7, 1.15); ax_eq.set_ylim(-0.25, 0.25)
ax_eq.set_title('Equator zoom\n(shape-squishy bulges more here)', fontsize=9)

sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=norm_q)
sm.set_array([])
fig2.colorbar(sm, ax=axes2, label='log₁₀(q)', fraction=0.012, pad=0.01)
fig2.suptitle(f'Perimeter shape at max strain  (b={B_ACTUAL:.2f}, τ={TAU:.3f})\n'
              f'Green=area-squishy (small q)  →  Red=shape-squishy (large q)',
              fontsize=10)
plt.tight_layout()
p2 = OUTDIR / 'perimeter_overlay.png'
fig2.savefig(p2, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {p2}')


# ── Movies: 5 representative q values ────────────────────────────────────────
MOVIE_Q = [0.1, 0.5, 2.0, 10.0, 50.0]
movie_results = [r for r in results if r['q'] in MOVIE_Q]

L0_vis    = 2.0 * np.pi / N       # R0=1
r_c_vis   = L0_vis
half_w_vis = 1.0 + 2.0 * L0_vis   # same as in run_squeeze_raw geometry


def save_comparison(cases, title, outpath):
    n = len(cases)
    n_frames_movie = min(len(r['frames']) for r in cases)
    fig, axes = plt.subplots(1, n, figsize=(5.0 * n, 6.0))
    panels = []
    for ax, r in zip(axes, cases):
        frames   = r['frames']
        lim_y    = 2.0 + 4.0 * r_c_vis + 0.15
        lim_x    = half_w_vis + r_c_vis + 0.15
        ax.set_xlim(-lim_x, lim_x); ax.set_ylim(-lim_y, lim_y)
        ax.set_aspect('equal')
        qb = r['q'] * B_ACTUAL
        ax.set_title(f"q={r['q']}  q·b={qb:.2f}\n"
                     f"ΔA@8%={r['dA_ref']*100:.2f}%  ν@8%={r['nu_ref']:.2f}", fontsize=9)
        ax.set_xlabel('x / R₀', fontsize=9)
        if ax is axes[0]:
            ax.set_ylabel('y / R₀', fontsize=9)

        # Wall patches (positions updated each frame)
        wall_h = 2.0 * r_c_vis
        top_p = mpatches.Rectangle((-half_w_vis - r_c_vis, 0),
                                   2 * (half_w_vis + r_c_vis), wall_h,
                                   fc='#888', ec='#555', lw=0.5, alpha=0.75, zorder=1)
        bot_p = mpatches.Rectangle((-half_w_vis - r_c_vis, 0),
                                   2 * (half_w_vis + r_c_vis), wall_h,
                                   fc='#888', ec='#555', lw=0.5, alpha=0.75, zorder=1)
        ax.add_patch(top_p); ax.add_patch(bot_p)

        # Disk polygons
        f1 = MplPolygon(np.zeros((4, 2)), closed=True,
                        fc='cornflowerblue', ec='navy', lw=1.2, alpha=0.60, zorder=3)
        f2 = MplPolygon(np.zeros((4, 2)), closed=True,
                        fc='salmon', ec='darkred', lw=1.2, alpha=0.60, zorder=3)
        ax.add_patch(f1); ax.add_patch(f2)

        # Contact force scatter (colored by |f_contact| per node)
        fmag_max = max(
            np.linalg.norm(fr['f_contact'], axis=-1).max() for fr in frames
        )
        sc = ax.scatter([], [], c=[], cmap='hot_r', vmin=0, vmax=max(fmag_max, 1e-10),
                        s=22, zorder=5, linewidths=0)
        tx = ax.text(0.03, 0.98, '', transform=ax.transAxes,
                     fontsize=7.5, va='top', family='monospace', zorder=10)
        panels.append(dict(f1=f1, f2=f2, top=top_p, bot=bot_p, sc=sc, tx=tx, frames=frames))

    fig.suptitle(title, fontsize=10)
    plt.tight_layout()

    def init():
        out = []
        for p in panels:
            p['f1'].set_xy(np.zeros((4, 2))); p['f2'].set_xy(np.zeros((4, 2)))
            p['sc'].set_offsets(np.zeros((0, 2))); p['sc'].set_array(np.zeros(0))
            p['tx'].set_text('')
            out += [p['f1'], p['f2'], p['sc'], p['tx'], p['top'], p['bot']]
        return out

    def animate(i):
        out = []
        for p in panels:
            fi  = min(i, len(p['frames']) - 1)
            fr  = p['frames'][fi]
            x_all = fr['x_all']            # (2, N, 2)
            p['f1'].set_xy(np.vstack([x_all[0], x_all[0][:1]]))
            p['f2'].set_xy(np.vstack([x_all[1], x_all[1][:1]]))
            p['top'].set_xy((-half_w_vis - r_c_vis, fr['y_top'] - r_c_vis))
            p['bot'].set_xy((-half_w_vis - r_c_vis, fr['y_bot']))
            # Node contact force magnitudes
            fmags = np.linalg.norm(fr['f_contact'], axis=-1)  # (2, N)
            all_pts = x_all.reshape(-1, 2)
            p['sc'].set_offsets(all_pts)
            p['sc'].set_array(fmags.ravel())
            p['tx'].set_text(
                f"ε_w={fr['wall_strain']:.3f}\n"
                f"ε_v={fr['eps_v']:.3f}\n"
                f"ν ={fr['nu_meas']:.3f}\n"
                f"ΔA={fr['dA_frac']*100:.2f}%"
            )
            out += [p['f1'], p['f2'], p['sc'], p['tx'], p['top'], p['bot']]
        return out

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=n_frames_movie, interval=1000 / FPS, blit=False)
    anim.save(str(outpath), writer=animation.PillowWriter(fps=FPS))
    plt.close(fig)
    print(f'    → {outpath.name}  ({outpath.stat().st_size // 1024} KB)')


print('\nRendering comparison movie...')
save_comparison(
    movie_results,
    title=f'Fixed b={B_ACTUAL:.2f} (τ={TAU:.3f}): q sweeps the squishiness axis',
    outpath=OUTDIR / 'comparison_q_sweep.gif',
)

print(f"\nAll done. Files in: {OUTDIR}/")
for f in sorted(OUTDIR.glob("*.*")):
    print(f"  {f.name:45s}  {f.stat().st_size//1024:4d} KB")
