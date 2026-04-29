"""Quick wrap-check: seed only + 300 steps. Should finish in ~2 min (mostly XLA compile)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np
import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Channel, SquareObstacle
from src.epd.motion import MotionSpec
from src.epd.system import System

Lx, Ly   = 20.0, 20.0
omega_box = 0.2
vx_box    = 0.5 / 3.0
vy_box    = 0.3 / 3.0

sys_q = System(Lx, Ly, periodic_x=True)
channel = Channel(width=Lx, height=Ly, x0=Lx/2, y0=Ly/2, exclusion='exterior')
box = SquareObstacle(side=5.0, x0=10.0, y0=10.0, exclusion='interior')
box.set_motion(MotionSpec(omega=omega_box, vx=vx_box, vy=vy_box, r_ref=(10.0, 10.0)))
sys_q.add_object(channel)
sys_q.add_object(box)

spec = ParticleSpec(count=25, nu=0.2, N_nodes=48, poly_dist=0.05)
sys_q.add_particles(spec)

sys_q.initialize(phi_target=0.78, seed=7, verbose=True, n_relax=50,
                 swell_alpha=10.0)

# 300 steps total, 10 frames → each chunk = 30 steps
n_steps      = 1500
sample_every = 150
print(f"dt={sys_q._dt:.5f}  running {n_steps} steps → {n_steps//sample_every} frames")

sys_q.run(n_steps, sample_every=sample_every, verbose=True)
print(f"Done. t={sys_q.t:.4f}  frames={len(sys_q.frames)}")

# Verify: all snapshots have CMs in [0, Lx)
x_cms = np.array([fr['x_cm'] for fr in sys_q.frames])
bad = np.sum((x_cms[:, :, 0] < 0) | (x_cms[:, :, 0] >= Lx))
print(f"CMs outside [0,Lx): {bad}  ({'OK' if bad==0 else 'FAIL'})")

# ── render GIF ────────────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Polygon as MplPolygon

colors  = sys_q._particle_colors
r_c_arr = sys_q._params['r_c_per_p'].numpy()
P       = len(sys_q._particles)

rLx = sys_q.Lx
rLy = sys_q.Ly

fig, ax = plt.subplots(figsize=(5, 5))
ax.set_xlim(-0.5, rLx + 0.5); ax.set_ylim(-0.5, rLy + 0.5)
ax.set_aspect('equal')
ax.plot([0, rLx, rLx, 0, 0], [0, 0, rLy, rLy, 0], 'b--', lw=0.6, alpha=0.4)

for d in channel.resolved(0.0):
    p = d['prim']
    ax.plot([p.p0[0], p.p1[0]], [p.p0[1], p.p1[1]], color='#333333', lw=2)

box_patch = MplPolygon(box.region_polygon(0.0)['vertices'], closed=True,
                       facecolor='#ffcccc', edgecolor='#aa4444', lw=2, alpha=0.85, zorder=2)
ax.add_patch(box_patch)

x_offsets = [-rLx, 0.0, rLx]
p_patches = []
for pi in range(P):
    c = colors[pi] if pi < len(colors) else 'cornflowerblue'
    row = []
    for _ in x_offsets:
        pat = MplPolygon(np.zeros((4, 2)), closed=True,
                         fc=c, ec='k', lw=0.8, alpha=0.7, zorder=3, visible=False)
        ax.add_patch(pat)
        row.append(pat)
    p_patches.append(row)

all_patches = [p for row in p_patches for p in row]
time_txt = ax.text(0.02, 0.97, '', transform=ax.transAxes, fontsize=8, va='top', family='monospace')

def _update(i):
    fr    = sys_q.frames[i]
    t_now = fr['t']
    box_patch.set_xy(box.region_polygon(t_now)['vertices'])
    for pi in range(P):
        outer = System._outer_contour(fr['x_all'][pi], fr['x_cm'][pi], float(r_c_arr[pi]))
        for k, odx in enumerate(x_offsets):
            sh  = outer + np.array([odx, 0.0])
            pat = p_patches[pi][k]
            if sh[:, 0].max() > -1.5 and sh[:, 0].min() < rLx + 1.5:
                pat.set_xy(sh); pat.set_visible(True)
            else:
                pat.set_visible(False)
    time_txt.set_text(f't={t_now:.4f}s')
    return all_patches + [box_patch, time_txt]

ani = animation.FuncAnimation(fig, _update, frames=len(sys_q.frames), interval=150, blit=True)
out = 'results/test_B2_quick.gif'
os.makedirs('results', exist_ok=True)
ani.save(out, writer='pillow', fps=6)
plt.close(fig)
print(f"GIF → {out}")
