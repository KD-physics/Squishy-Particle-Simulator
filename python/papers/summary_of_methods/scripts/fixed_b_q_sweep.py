"""
Fixed b=0.2 sweep: vary q only.

b = EI/El_t = τ²/12 = 0.2  →  τ = sqrt(12 * 0.2) = sqrt(2.4) ≈ 1.549

At this b, q sweeps the full squishiness axis:
  q small  →  ΔA large,  ν → 0     (area-compressible: area shrinks, shape barely changes)
  q large  →  ΔA → 0,    ν → ~0.9  (shape-compressible: area conserved, shape deforms)

This is the canonical 1D squishiness axis with b=0.2 as the fixed bending context.

q·b = 1 occurs at q = 1/0.2 = 5  (natural midpoint of the axis)
"""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.animation as animation
from pathlib import Path

from twodisk_capsule import run_squeeze_once, _capsule_outline_polygon

OUTDIR = pathlib.Path("results/fixed_b_q_sweep")
OUTDIR.mkdir(parents=True, exist_ok=True)

S        = 1.0
R0       = 1.0
N        = 32
RHO_F    = 1.0
ALPHA0   = 2.0
SR       = 0.01
N_FRAMES = 80
DELTA_MAX = 1.5
EPS_MAX  = 0.12
C_0      = 3000.0
EPS_REF  = 0.08
FPS      = 15

B_TARGET = 0.2
TAU      = np.sqrt(12.0 * B_TARGET)        # ≈ 1.5492
B_ACTUAL = TAU**2 / 12.0                   # = 0.2 exactly

Q_VALS = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
# q·b=1 at q=5.0

print(f"Fixed b = {B_ACTUAL:.4f}  (τ = {TAU:.4f})")
print(f"q·b = 1  at  q = {1/B_ACTUAL:.1f}")
print()


def make_params(q):
    El_t   = 12.0 * S / TAU**2
    K_area = q * El_t
    C      = C_0 * S * (1.0 + q)
    EI     = S * R0**3
    return dict(El_t=El_t, K_area=K_area, C=C, EI=EI, b=B_ACTUAL, q_b=q*B_ACTUAL)


def interp_at_eps(eps_arr, val_arr, eps_target):
    eps_arr = np.asarray(eps_arr); val_arr = np.asarray(val_arr)
    valid = np.isfinite(val_arr)
    if not valid.any(): return float('nan')
    e, v = eps_arr[valid], val_arr[valid]
    if eps_target < e[0] or eps_target > e[-1]: return float('nan')
    idx = np.clip(np.searchsorted(e, eps_target), 1, len(e)-1)
    t = (eps_target - e[idx-1]) / (e[idx] - e[idx-1] + 1e-30)
    return float(v[idx-1] + t*(v[idx] - v[idx-1]))


# ── Run sweep ─────────────────────────────────────────────────────────────────
print(f"{'q':>6}  {'q·b':>6}  {'ΔA@8%':>8}  {'ν@8%':>8}  {'θ_w@8%':>9}  {'cc%':>5}")
print("-"*52)

results = []
for q in Q_VALS:
    p = make_params(q)
    try:
        frames, dt, T_wave, v_wall, p1, p2, worst_pen, worst_cc_gap = run_squeeze_once(
            R0=R0, N=N, tau=TAU, S=S, C=p['C'], rho_f=RHO_F,
            strain_rate_ratio=SR, alpha_damp=1.0, alpha0=ALPHA0,
            n_frames=N_FRAMES, delta_max=DELTA_MAX, dt_factor=0.4,
            verbose=False, K_area=p['K_area'], side_walls=False, eps_max=EPS_MAX,
        )
    except Exception as e:
        print(f"{q:6.2f}  CRASHED: {e}"); continue

    A0     = frames[0]['A0']
    r_c    = p1.r_c
    cc_pct = abs(worst_cc_gap) / r_c * 100.0

    eps_arr = np.array([0.5*(fr['eps1']+fr['eps2']) for fr in frames])
    dA_arr  = np.array([(1-0.5*(fr['A1']+fr['A2'])/A0)*100 for fr in frames])
    nu_arr  = np.array([0.5*(fr.get('nu_meas1',np.nan)+fr.get('nu_meas2',np.nan)) for fr in frames])
    tw_arr  = np.array([fr.get('turn_angle_wall_deg', np.nan) for fr in frames])

    dA_ref = interp_at_eps(eps_arr, dA_arr, EPS_REF)
    nu_ref = interp_at_eps(eps_arr, nu_arr, EPS_REF)
    tw_ref = interp_at_eps(eps_arr, tw_arr, EPS_REF)

    x1_c = frames[-1]['x1'] - frames[-1]['x1'].mean(axis=0)

    print(f"{q:6.2f}  {p['q_b']:6.3f}  {dA_ref:8.4f}  {nu_ref:8.4f}  {tw_ref:9.3f}  {cc_pct:5.1f}")

    results.append(dict(q=q, **p,
                        eps_arr=eps_arr, dA_arr=dA_arr, nu_arr=nu_arr, tw_arr=tw_arr,
                        dA_ref=dA_ref, nu_ref=nu_ref, tw_ref=tw_ref, cc_pct=cc_pct,
                        frames=frames, x1_c=x1_c))

(OUTDIR / 'fixed_b_results.json').write_text(
    json.dumps([{k: (v.tolist() if isinstance(v, np.ndarray) else
                     (v if not isinstance(v, list) or k not in ('frames','x1_c') else None))
                 for k, v in r.items() if k not in ('frames',)} for r in results], indent=2))


# ── Figure 1: ΔA and ν vs q  ──────────────────────────────────────────────────
q_vals  = np.array([r['q']      for r in results])
dA_vals = np.array([r['dA_ref'] for r in results])
nu_vals = np.array([r['nu_ref'] for r in results])
qb_vals = q_vals * B_ACTUAL

# Color by q·b position relative to 1
cmap_q = plt.cm.RdYlGn
norm_q = plt.Normalize(vmin=np.log10(q_vals.min()), vmax=np.log10(q_vals.max()))

fig1, axes1 = plt.subplots(1, 3, figsize=(16, 5))
ax_q, ax_qb, ax_phase = axes1

# ΔA and ν vs q (linear)
ax2 = ax_q.twinx()
l1, = ax_q.semilogx(q_vals, dA_vals, 'b^-', ms=9, lw=2, label='ΔA (%)')
l2, = ax2.semilogx(q_vals, nu_vals,  'ro-', ms=9, lw=2, label='ν_meas')
ax_q.axvline(1/B_ACTUAL, color='gray', ls='--', lw=1.2, label=f'q·b=1 (q={1/B_ACTUAL:.0f})')
ax_q.set_xlabel('q', fontsize=12); ax_q.set_ylabel('ΔA (%) at ε=0.08', color='blue', fontsize=11)
ax2.set_ylabel('ν at ε=0.08', color='red', fontsize=11)
ax_q.set_title(f'ΔA and ν vs q\n(b={B_ACTUAL:.2f}, τ={TAU:.3f})', fontsize=10)
ax_q.legend(handles=[l1,l2], fontsize=9); ax_q.grid(True, alpha=0.3)

# ΔA and ν vs q·b (shows where q·b=1 falls naturally)
ax2b = ax_qb.twinx()
l3, = ax_qb.semilogx(qb_vals, dA_vals, 'b^-', ms=9, lw=2, label='ΔA (%)')
l4, = ax2b.semilogx(qb_vals, nu_vals,  'ro-', ms=9, lw=2, label='ν_meas')
ax_qb.axvline(1.0, color='gray', ls='--', lw=1.5, label='q·b = 1')
ax_qb.set_xlabel('q · b', fontsize=12); ax_qb.set_ylabel('ΔA (%)', color='blue', fontsize=11)
ax2b.set_ylabel('ν', color='red', fontsize=11)
ax_qb.set_title('Same axes rescaled to q·b\n(natural midpoint at q·b=1)', fontsize=10)
ax_qb.legend(handles=[l3,l4], fontsize=9); ax_qb.grid(True, alpha=0.3)

# ΔA–ν trajectory as q varies
sc = ax_phase.scatter(dA_vals, nu_vals,
                      c=np.log10(q_vals), cmap='RdYlGn',
                      s=120, zorder=5, edgecolors='k', linewidths=0.5)
ax_phase.plot(dA_vals, nu_vals, 'k-', lw=1.0, alpha=0.4, zorder=3)
plt.colorbar(sc, ax=ax_phase, label='log₁₀(q)')
# Annotate key points
for r in results:
    if r['q'] in [0.1, 1.0, 5.0, 20.0]:
        ax_phase.annotate(f"q={r['q']}", (r['dA_ref'], r['nu_ref']),
                          textcoords='offset points', xytext=(6,4), fontsize=8)
# Mark q·b=1
qb1_pts = [r for r in results if abs(r['q']*B_ACTUAL - 1.0) < 0.3]
if qb1_pts:
    r = min(qb1_pts, key=lambda r: abs(r['q']*B_ACTUAL - 1.0))
    ax_phase.scatter([r['dA_ref']], [r['nu_ref']], s=250, facecolors='none',
                     edgecolors='red', lw=2.5, zorder=10, label='q·b≈1')
ax_phase.axhline(0, color='lightgray', ls=':', lw=0.8)
ax_phase.axhline(1, color='gray',      ls=':', lw=0.8, label='ν=1')
ax_phase.set_xlabel('ΔA (%) at ε=0.08', fontsize=12)
ax_phase.set_ylabel('ν_meas at ε=0.08', fontsize=12)
ax_phase.set_title('ΔA – ν trajectory  (b=0.2 fixed, q varies)\nq·b=1 circled in red', fontsize=10)
ax_phase.legend(fontsize=8); ax_phase.grid(True, alpha=0.3)

fig1.suptitle(f'Fixed b = {B_ACTUAL:.2f}  (τ = {TAU:.3f})  |  q sweeps the squishiness axis\n'
              f'q·b = 1  at  q = {1/B_ACTUAL:.0f}  |  left=area-squishy, right=shape-squishy',
              fontsize=11)
plt.tight_layout()
p1 = OUTDIR / 'fixed_b_sweep.png'
fig1.savefig(p1, dpi=150, bbox_inches='tight')
plt.close()
print(f'\nSaved: {p1}')


# ── Figure 2: Perimeter overlay at max strain ─────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(17, 5.5))
ax_full, ax_top, ax_eq = axes2

def ref_circle(ax):
    th = np.linspace(0, 2*np.pi, 300)
    ax.plot(np.cos(th), np.sin(th), ':', color='#cccccc', lw=1.0, zorder=0)

for ax in axes2:
    ref_circle(ax); ax.set_aspect('equal'); ax.grid(True, alpha=0.15)
    ax.set_xlabel('x / R₀', fontsize=10)

for r in results:
    c   = cmap_q(norm_q(np.log10(r['q'])))
    lbl = f"q={r['q']}  ΔA={r['dA_arr'][-1]:.2f}%  ν={r['nu_arr'][np.isfinite(r['nu_arr'])][-1]:.2f}"
    x   = r['x1_c']
    pts = np.vstack([x, x[:1]])
    ax_full.plot(pts[:,0], pts[:,1], color=c, lw=2.0, label=lbl)

    mask_top = x[:,1] > 0.3
    ix = np.where(mask_top)[0]
    if len(ix)>1: ax_top.plot(x[ix,0], x[ix,1], color=c, lw=2.2)

    mask_eq = np.abs(x[:,1]) < 0.25
    ix = np.where(mask_eq)[0]
    if len(ix)>1: ax_eq.plot(x[ix,0], x[ix,1], color=c, lw=2.2, marker='o', ms=3)

ax_full.set_xlim(-1.3,1.3); ax_full.set_ylim(-1.3,1.3)
ax_full.set_title(f'Full perimeter at max strain\n(b={B_ACTUAL:.2f}, q: green=small → red=large)', fontsize=9)
ax_full.legend(fontsize=6, loc='lower right'); ax_full.set_ylabel('y / R₀', fontsize=10)

ax_top.set_xlim(-1.2,1.2); ax_top.set_ylim(0.3,1.25)
ax_top.set_title('Wall-contact side (top half)', fontsize=9)

ax_eq.set_xlim(0.7,1.15); ax_eq.set_ylim(-0.25,0.25)
ax_eq.set_title('Equator zoom\n(shape-squishy bulges more here)', fontsize=9)

sm = plt.cm.ScalarMappable(cmap='RdYlGn', norm=norm_q)
sm.set_array([])
fig2.colorbar(sm, ax=axes2, label='log₁₀(q)', fraction=0.012, pad=0.01)
fig2.suptitle(f'Perimeter shape at max strain  (b={B_ACTUAL:.2f}, τ={TAU:.3f})\n'
              f'Green=area-squishy (small q, large ΔA)  →  Red=shape-squishy (large q, small ΔA)',
              fontsize=10)
plt.tight_layout()
p2 = OUTDIR / 'perimeter_overlay.png'
fig2.savefig(p2, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {p2}')


# ── Movies: select 5 representative q values ─────────────────────────────────
MOVIE_Q = [0.1, 0.5, 2.0, 10.0, 50.0]   # spans the full axis

def save_movie(r, outpath):
    frames = r['frames']
    p      = make_params(r['q'])
    r_c      = frames[0]['r_c'];  r_c_wall = frames[0]['r_c_wall']
    half_w   = frames[0]['half_w'];  A0 = frames[0]['A0']
    dR_all   = [0.5*(fr.get('eps1',0)+fr.get('eps2',0)) for fr in frames]
    F_all    = [np.linalg.norm(fr['F_top']) for fr in frames]

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 5.5))
    lim_y = 2.2*R0+4.0*r_c_wall; lim_x = half_w+2.0*r_c_wall+0.15
    ax_l.set_xlim(-lim_x,lim_x); ax_l.set_ylim(-lim_y,lim_y); ax_l.set_aspect('equal')
    ax_l.set_xlabel('x / R₀', fontsize=10); ax_l.set_ylabel('y / R₀', fontsize=10)
    ax_r.set_xlim(-0.002, max(dR_all)*1.12)
    ax_r.set_ylim(-0.02, max(F_all)*1.25+0.02)
    ax_r.set_xlabel('ε_p', fontsize=10); ax_r.set_ylabel('|F_top|', fontsize=10)
    ax_r.set_title('Force vs strain', fontsize=10); ax_r.grid(True, alpha=0.3)

    qb = r['q'] * B_ACTUAL
    fig.suptitle(f'b={B_ACTUAL:.2f}  |  q={r["q"]},  q·b={qb:.2f}\n'
                 f'El_t={p["El_t"]:.2f},  K_area={p["K_area"]:.1f},  τ={TAU:.3f}', fontsize=10)

    top_p = mpatches.Rectangle((-half_w-r_c_wall,0),2*(half_w+r_c_wall),2*r_c_wall,
                                fc='#888',ec='#555',lw=0.5,alpha=0.75,zorder=1)
    bot_p = mpatches.Rectangle((-half_w-r_c_wall,0),2*(half_w+r_c_wall),2*r_c_wall,
                                fc='#888',ec='#555',lw=0.5,alpha=0.75,zorder=1)
    ax_l.add_patch(top_p); ax_l.add_patch(bot_p)
    f1 = MplPolygon(np.zeros((4,2)),closed=True,fc='cornflowerblue',ec='navy',lw=1.2,alpha=0.60,zorder=3)
    f2 = MplPolygon(np.zeros((4,2)),closed=True,fc='salmon',ec='darkred',lw=1.2,alpha=0.60,zorder=3)
    ax_l.add_patch(f1); ax_l.add_patch(f2)
    fmag_max = max((fr['contact_fmags'].max() if len(fr['contact_fmags'])>0 else 0
                    for fr in frames), default=1.0)
    scat = ax_l.scatter([],[],c=[],cmap='hot_r',vmin=0,vmax=max(fmag_max,1e-10),
                        s=25,zorder=5,linewidths=0)
    fig.colorbar(scat,ax=ax_l,fraction=0.03,pad=0.01).set_label('contact force',fontsize=8)
    txt  = ax_l.text(0.02,0.98,'',transform=ax_l.transAxes,
                     fontsize=8.5,va='top',family='monospace',zorder=10)
    fline,  = ax_r.plot([],[],'k-',lw=1.5)
    cur_pt, = ax_r.plot([],[],'ro',ms=7)
    plt.tight_layout()

    def init():
        f1.set_xy(np.zeros((4,2))); f2.set_xy(np.zeros((4,2)))
        scat.set_offsets(np.zeros((0,2))); scat.set_array(np.zeros(0))
        fline.set_data([],[]); cur_pt.set_data([],[]); txt.set_text('')
        return f1,f2,scat,fline,cur_pt,txt,top_p,bot_p

    def animate(i):
        fr = frames[i]
        f1.set_xy(_capsule_outline_polygon(fr['x1'],r_c))
        f2.set_xy(_capsule_outline_polygon(fr['x2'],r_c))
        top_p.set_xy((-half_w-r_c_wall,fr['y_top']-r_c_wall))
        bot_p.set_xy((-half_w-r_c_wall,fr['y_bot']-r_c_wall))
        cpts=fr['contact_pts']; cfmag=fr['contact_fmags']
        if len(cpts)>0: scat.set_offsets(cpts); scat.set_array(cfmag)
        else: scat.set_offsets(np.zeros((0,2))); scat.set_array(np.zeros(0))
        eps_p  = 0.5*(fr.get('eps1',0)+fr.get('eps2',0))
        dA_pct = (1-0.5*(fr['A1']+fr['A2'])/A0)*100
        nu1    = fr.get('nu_meas1',np.nan); nu2 = fr.get('nu_meas2',np.nan)
        nu_m   = 0.5*(nu1+nu2) if np.isfinite(nu1) else np.nan
        tw_deg = fr.get('turn_angle_wall_deg',np.nan)
        txt.set_text(f'ε_p  = {eps_p:.4f}\n|F|  = {np.linalg.norm(fr["F_top"]):.4f}\n'
                     f'ΔA   = {dA_pct:.3f}%\nν    = {nu_m:.3f}\nθ_w  = {tw_deg:.2f}°')
        fline.set_data(dR_all[:i+1],F_all[:i+1]); cur_pt.set_data([dR_all[i]],[F_all[i]])
        return f1,f2,scat,fline,cur_pt,txt,top_p,bot_p

    anim = animation.FuncAnimation(fig,animate,init_func=init,
                                   frames=len(frames),interval=1000/FPS,blit=False)
    anim.save(str(outpath),writer=animation.PillowWriter(fps=FPS))
    plt.close(fig)
    print(f'    → {outpath.name}  ({outpath.stat().st_size//1024} KB)')


print('\nRendering movies...')
movie_results = [r for r in results if r['q'] in MOVIE_Q]
for r in movie_results:
    out = OUTDIR / f'q{str(r["q"]).replace(".","p")}.gif'
    save_movie(r, out)


# ── 5-panel comparison movie ───────────────────────────────────────────────────
def save_comparison(cases, title, outpath):
    n = len(cases); n_frames = min(len(r['frames']) for r in cases)
    fig, axes = plt.subplots(1, n, figsize=(5.0*n, 6.0))
    panels = []
    for ax, r in zip(axes, cases):
        frames   = r['frames']
        r_c      = frames[0]['r_c']; r_c_wall = frames[0]['r_c_wall']
        half_w   = frames[0]['half_w']; A0 = frames[0]['A0']
        lim_y = 2.2*R0+4.0*r_c_wall; lim_x = half_w+2.0*r_c_wall+0.15
        ax.set_xlim(-lim_x,lim_x); ax.set_ylim(-lim_y,lim_y); ax.set_aspect('equal')
        qb = r['q'] * B_ACTUAL
        ax.set_title(f'q={r["q"]}  q·b={qb:.2f}\n'
                     f'ΔA@8%={r["dA_ref"]:.2f}%  ν@8%={r["nu_ref"]:.2f}', fontsize=9)
        ax.set_xlabel('x / R₀', fontsize=9)
        if ax is axes[0]: ax.set_ylabel('y / R₀', fontsize=9)
        fmag_max = max((fr['contact_fmags'].max() if len(fr['contact_fmags'])>0 else 0
                        for fr in frames), default=1.0)
        top_p = mpatches.Rectangle((-half_w-r_c_wall,0),2*(half_w+r_c_wall),2*r_c_wall,
                                    fc='#888',ec='#555',lw=0.5,alpha=0.75,zorder=1)
        bot_p = mpatches.Rectangle((-half_w-r_c_wall,0),2*(half_w+r_c_wall),2*r_c_wall,
                                    fc='#888',ec='#555',lw=0.5,alpha=0.75,zorder=1)
        ax.add_patch(top_p); ax.add_patch(bot_p)
        f1 = MplPolygon(np.zeros((4,2)),closed=True,fc='cornflowerblue',ec='navy',lw=1.2,alpha=0.60,zorder=3)
        f2 = MplPolygon(np.zeros((4,2)),closed=True,fc='salmon',ec='darkred',lw=1.2,alpha=0.60,zorder=3)
        ax.add_patch(f1); ax.add_patch(f2)
        sc = ax.scatter([],[],c=[],cmap='hot_r',vmin=0,vmax=max(fmag_max,1e-10),
                        s=22,zorder=5,linewidths=0)
        tx = ax.text(0.03,0.98,'',transform=ax.transAxes,
                     fontsize=7.5,va='top',family='monospace',zorder=10)
        panels.append(dict(f1=f1,f2=f2,top=top_p,bot=bot_p,sc=sc,tx=tx,
                           frames=frames,r_c=r_c,r_c_wall=r_c_wall,half_w=half_w,A0=A0))

    fig.suptitle(title, fontsize=10)
    plt.tight_layout()

    def init():
        out = []
        for p in panels:
            p['f1'].set_xy(np.zeros((4,2))); p['f2'].set_xy(np.zeros((4,2)))
            p['sc'].set_offsets(np.zeros((0,2))); p['sc'].set_array(np.zeros(0))
            p['tx'].set_text('')
            out += [p['f1'],p['f2'],p['sc'],p['tx'],p['top'],p['bot']]
        return out

    def animate(i):
        out = []
        for p in panels:
            fi=min(i,len(p['frames'])-1); fr=p['frames'][fi]
            r_c=p['r_c']; r_c_wall=p['r_c_wall']; hw=p['half_w']; A0=p['A0']
            p['f1'].set_xy(_capsule_outline_polygon(fr['x1'],r_c))
            p['f2'].set_xy(_capsule_outline_polygon(fr['x2'],r_c))
            p['top'].set_xy((-hw-r_c_wall,fr['y_top']-r_c_wall))
            p['bot'].set_xy((-hw-r_c_wall,fr['y_bot']-r_c_wall))
            cpts=fr['contact_pts']; cfmag=fr['contact_fmags']
            if len(cpts)>0: p['sc'].set_offsets(cpts); p['sc'].set_array(cfmag)
            else: p['sc'].set_offsets(np.zeros((0,2))); p['sc'].set_array(np.zeros(0))
            eps_p  = 0.5*(fr.get('eps1',0)+fr.get('eps2',0))
            dA_pct = (1-0.5*(fr['A1']+fr['A2'])/A0)*100
            nu1=fr.get('nu_meas1',np.nan); nu_m=0.5*(nu1+fr.get('nu_meas2',np.nan))
            tw_deg = fr.get('turn_angle_wall_deg',np.nan)
            p['tx'].set_text(f'ε_p={eps_p:.4f}\nΔA={dA_pct:.3f}%\nν={nu_m:.3f}\nθ_w={tw_deg:.2f}°')
            out += [p['f1'],p['f2'],p['sc'],p['tx'],p['top'],p['bot']]
        return out

    anim = animation.FuncAnimation(fig,animate,init_func=init,
                                   frames=n_frames,interval=1000/FPS,blit=False)
    anim.save(str(outpath),writer=animation.PillowWriter(fps=FPS))
    plt.close(fig)
    print(f'    → {outpath.name}  ({outpath.stat().st_size//1024} KB)')


save_comparison(
    movie_results,
    title=f'Fixed b={B_ACTUAL:.2f} (τ={TAU:.3f}): q sweeps the squishiness axis  '
          f'|  q·b=1 at q={1/B_ACTUAL:.0f}',
    outpath=OUTDIR / 'comparison_q_sweep.gif',
)

print(f'\nAll done. Files in: {OUTDIR}/')
print('\nFiles:')
for f in sorted(OUTDIR.glob('*.*')):
    print(f'  {f.name:45s}  {f.stat().st_size//1024:4d} KB')
