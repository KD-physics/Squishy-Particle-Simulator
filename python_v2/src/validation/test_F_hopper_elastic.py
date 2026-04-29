"""
test_F_hopper_elastic.py — Elastic-shell hopper discharge, spot-check (25k steps).

50 elastic capsules (ν=0.5, N=32, 5% polydispersity) falling through the same
30°-walled hopper as test_F_hopper (W_out=4, g=0.05).

Run modes
---------
    QUICK=1 python src/validation/test_F_hopper_elastic.py   # 2 outer batches
    python src/validation/test_F_hopper_elastic.py           # full 25k steps
"""

import sys, os, time
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import io

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import HopperRegion
from src.epd.system import System

# ── Paths ─────────────────────────────────────────────────────────────────────

OUT_DIR   = "results"
os.makedirs(OUT_DIR, exist_ok=True)
GIF_PATH  = os.path.join(OUT_DIR, "F1_hopper_elastic.gif")
SNAP_PATH = os.path.join(OUT_DIR, "F1_hopper_elastic_init.png")

# ── Mode ──────────────────────────────────────────────────────────────────────

QUICK     = os.environ.get("QUICK",     "0") == "1"
INIT_ONLY = os.environ.get("INIT_ONLY", "0") == "1"

if QUICK:
    N_PARTICLES  = 10
    H_RES        = 15.0
    LY           = 25.0
    N_RELAX      = 20
    OUTER_BATCH  = 600
    SAMPLE_EVERY = 300
    N_OUTER      = 2
    MIN_DELTA_Y  = 0.0
else:
    N_PARTICLES  = 50
    H_RES        = 50.0
    LY           = 66.0
    N_RELAX      = 0
    OUTER_BATCH  = 5_000
    SAMPLE_EVERY = 500
    N_OUTER      = 5        # 5 × 5000 = 25 000 steps
    MIN_DELTA_Y  = 0.5

# ── Fixed physics parameters ──────────────────────────────────────────────────

LX        = 12.0
X_C       = LX / 2.0
W_OUT     = 4.0
W_RES     = 12.0
THETA_DEG = 30.0
Y_BOT     = 0.0
G_GRAV    = 0.05
NU        = 0.5     # effective Poisson ratio (moderate squishiness)
N_NODES   = 32
POLY_SIG  = 0.05

# ── Build system ──────────────────────────────────────────────────────────────

print("=" * 60)
print(f"Test F1-elastic — Gravity-driven hopper  ({'QUICK' if QUICK else 'FULL'} mode)")
print(f"  N={N_PARTICLES}  N_nodes={N_NODES}  ν={NU}  g={G_GRAV}")
print(f"  W_out={W_OUT}  W_res={W_RES}  θ={THETA_DEG}°  h_res={H_RES}")
print(f"  outer_batch={OUTER_BATCH}  sample_every={SAMPLE_EVERY}  "
      f"n_outer={N_OUTER}  (total={N_OUTER*OUTER_BATCH} steps)")
print("=" * 60)

sys_ = System(LX, LY, periodic_x=False, periodic_y=False, g=G_GRAV)
hopper = HopperRegion(X_C, W_OUT, W_RES, THETA_DEG, H_RES, Y_BOT)
sys_.add_object(hopper)
spec = ParticleSpec(count=N_PARTICLES, type='elastic',
                    nu=NU, N_nodes=N_NODES, poly_dist=POLY_SIG)
sys_.add_particles(spec)

print("\nInitialising ...")
sys_.initialize(phi_target=0.35, seed=42, verbose=True,
                relax_only=True, n_relax_init=(0 if INIT_ONLY else N_RELAX))

mean_y_init = float(sys_.state['x_cm'].numpy()[:, 1].mean())
print(f"  φ_outer={sys_.phi_outer:.4f}  mean_y={mean_y_init:.2f}")

_r_c_arr = sys_._params['r_c_per_p'].numpy()   # (P,)

# ── Render helper ─────────────────────────────────────────────────────────────

def _draw_frame(ax, x_all, x_cm, t):
    sys_._render_objects(ax, t=t)
    cmap = plt.cm.tab20
    for i in range(len(x_cm)):
        outer = sys_._outer_contour(x_all[i], x_cm[i], _r_c_arr[i])
        ax.fill(outer[:, 0], outer[:, 1],
                color=cmap(i % 20), alpha=0.80, zorder=2)
        xs = np.append(outer[:, 0], outer[0, 0])
        ys = np.append(outer[:, 1], outer[0, 1])
        ax.plot(xs, ys, 'k-', lw=0.4, zorder=3)

def _render(snap, label=""):
    x_all = snap['x_all']
    x_cm  = snap['x_cm']
    t     = snap['t']
    if np.isnan(x_cm).any():
        return None
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.set_xlim(-0.5, LX + 0.5);  ax.set_ylim(-2.0, 10.0);  ax.set_aspect('equal')
    _draw_frame(ax, x_all, x_cm, t)
    n_out = int(np.sum(x_cm[:, 1] < Y_BOT - 0.5))
    ax.set_title(f't={t:.0f}  n_exited={n_out}  {label}', fontsize=8)
    ax.axis('off');  fig.tight_layout(pad=0.2)
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=72)
    plt.close(fig);  buf.seek(0)
    return imageio.imread(buf)

# ── Initial snapshot ──────────────────────────────────────────────────────────

x_all_np = sys_.state['x_all'].numpy()
x_cm_np  = sys_.state['x_cm'].numpy()

margin = 2.0
x_lo_s = -margin;          x_hi_s = LX + margin
y_lo_s = Y_BOT - margin;   y_hi_s = hopper.h_total + Y_BOT + margin
aspect = (y_hi_s - y_lo_s) / (x_hi_s - x_lo_s)
fig, ax = plt.subplots(figsize=(5.0, 5.0 * aspect))
ax.set_xlim(x_lo_s, x_hi_s);  ax.set_ylim(y_lo_s, y_hi_s);  ax.set_aspect('equal')
_draw_frame(ax, x_all_np, x_cm_np, 0.0)
ax.set_title(f"t=0  N={N_PARTICLES}  φ_outer={sys_.phi_outer:.3f}  ν={NU}", fontsize=9)
ax.axis('off');  plt.tight_layout(pad=0.2)
plt.savefig(SNAP_PATH, dpi=120);  plt.close()
print(f"  Init snapshot → {SNAP_PATH}")

if INIT_ONLY:
    sys.exit(0)

# ── Diagnostics ───────────────────────────────────────────────────────────────

def _print_diag(batch_i, n_new_chunks, elapsed_batch):
    recent = sys_.diag[-n_new_chunks:]
    cand_checks  = sum(d.get('n_cand_checks',  0) for d in recent)
    cand_updates = sum(d.get('n_cand_updates', 0) for d in recent)
    retraces     = sum(d.get('n_retraces_step_full', 0) for d in recent)
    steps_done   = sum(d.get('n_steps', SAMPLE_EVERY) for d in recent)
    ms_per_step  = elapsed_batch / max(steps_done, 1) * 1000
    n_out = int(np.sum(sys_.state['x_cm'].numpy()[:, 1] < Y_BOT - 0.5))
    print(f"  [batch {batch_i+1:>2}/{N_OUTER}]  "
          f"t={sys_.t:.1f}  n_exited={n_out}  "
          f"cand_checks={cand_checks}  cand_updates={cand_updates}  "
          f"retraces={retraces}  {ms_per_step:.1f} ms/step  "
          f"({elapsed_batch:.0f}s)")

# ── Main loop ─────────────────────────────────────────────────────────────────

print(f"\nRunning {N_OUTER} outer batches of {OUTER_BATCH} steps "
      f"(sample_every={SAMPLE_EVERY}) ...")
print("GIF updated every outer batch →", GIF_PATH)
print()

n_chunks_per_batch = OUTER_BATCH // SAMPLE_EVERY
t_run_start = time.time()

for outer_i in range(N_OUTER):
    t_batch = time.time()
    sys_.run(OUTER_BATCH, sample_every=SAMPLE_EVERY,
             record_initial=(outer_i == 0), verbose=False)
    elapsed_batch = time.time() - t_batch
    _print_diag(outer_i, n_chunks_per_batch, elapsed_batch)

    imgs = [_render(fr, f"batch {outer_i+1}") for fr in sys_.frames]
    imgs = [im for im in imgs if im is not None]
    if imgs:
        imageio.mimwrite(GIF_PATH, imgs, duration=0.15, loop=0)

# ── Final checks ──────────────────────────────────────────────────────────────

total_elapsed = time.time() - t_run_start
mean_y_final = float(sys_.state['x_cm'].numpy()[:, 1].mean())
n_exited     = int(np.sum(sys_.state['x_cm'].numpy()[:, 1] < Y_BOT - 0.5))
nan_ok       = not np.isnan(sys_.state['x_cm'].numpy()).any()
delta_y      = mean_y_init - mean_y_final

print()
print("=" * 60)
print("Test F1-elastic — Final results")
print("=" * 60)
print(f"  total_steps  = {N_OUTER * OUTER_BATCH}")
print(f"  total_time   = {total_elapsed:.0f}s  ({total_elapsed/60:.1f} min)")
print(f"  mean_y       = {mean_y_init:.2f} → {mean_y_final:.2f}  (Δ={delta_y:.2f})")
print(f"  n_exited     = {n_exited}")
print(f"  NaN-free     = {nan_ok}")
print(f"  frames in GIF= {len(sys_.frames)}")
print(f"  GIF saved to {GIF_PATH}")

assert nan_ok, "NaN positions detected — FAIL"
assert delta_y >= MIN_DELTA_Y, \
    f"Particles moved up or not enough (Δy={delta_y:.3f} < {MIN_DELTA_Y}) — FAIL"

print("\n  ALL PASS")
