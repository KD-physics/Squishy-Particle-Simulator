"""
test_B2_spin.py — Test B2: SquareObstacle with spin + translation

Verifies:
  1. Geometry: box.region_polygon(t) corners match analytic R(omega*t)+v*t formula
     to machine precision at several t values.
  2. to_make_prim_list() encodes omega and r_ref correctly.
  3. Physics: after running, all particle CMs remain outside the moving obstacle
     at every recorded frame.
  4. GIF movie saved showing the spinning, translating box and particles.

PASS criteria:
  - |corner_pos(t) - analytic| < 1e-10  for all test times
  - omega and r_ref encoded correctly
  - 0 particle CMs inside rotated+translated obstacle at every frame

Usage:
    python src/validation/test_B2_spin.py [--out PATH]
"""

import sys, os, argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import Channel, SquareObstacle, _point_in_polygon
from src.epd.motion import MotionSpec
from src.epd.system import System

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=str, default='results/test_B2_spin.gif')
args = parser.parse_args()

os.makedirs(os.path.dirname(args.out), exist_ok=True)

print("=" * 60)
print("Test B2: SquareObstacle spin + translation")
print("=" * 60)

# ── Parameters ────────────────────────────────────────────────────────────────
Lx, Ly    = 20.0, 20.0
box_side  = 5.0
box_cx    = Lx / 2.0       # 10.0
box_cy    = Ly / 2.0       # 10.0
omega_box = 0.2             # rad/s — 1/5 of original (full rotation in ~31 s)
vx_box    = 0.5 / 3.0      # ≈ 0.167 units/s — 1/3 of original
vy_box    = 0.3 / 3.0      # = 0.100 units/s — 1/3 of original

# ── Part 1: geometry verification (pure Python) ────────────────────────────────
print("\n--- Part 1: geometry (region_polygon vs analytic formula) ---")

box_geo = SquareObstacle(side=box_side, x0=box_cx, y0=box_cy, exclusion='interior')
box_geo.set_motion(MotionSpec(omega=omega_box, vx=vx_box, vy=vy_box,
                               r_ref=(box_cx, box_cy)))

w2 = box_side / 2.0
corners_local = np.array([
    [-w2, -w2],
    [ w2, -w2],
    [ w2,  w2],
    [-w2,  w2],
])

def analytic_corners(t):
    """Corners of a box that spins (omega_box) and translates (vx_box, vy_box)."""
    theta = omega_box * t
    c, s  = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    ctr = np.array([box_cx + vx_box * t,
                    box_cy + vy_box * t])
    return (R @ corners_local.T).T + ctr   # (4, 2)

test_times = [0.5, 2.0, 5.0, 10.0, np.pi / omega_box]
geo_ok = True
for t_val in test_times:
    verts = box_geo.region_polygon(t_val)['vertices']
    ref   = analytic_corners(t_val)
    err   = np.max(np.abs(verts - ref))
    status = "PASS" if err < 1e-10 else "FAIL"
    if err >= 1e-10:
        geo_ok = False
    print(f"  t={t_val:5.2f}  max corner err = {err:.2e}  {status}")

if geo_ok:
    print("  PASS: all corner positions match analytic formula to machine precision")
else:
    print("  FAIL: corner position mismatch")

# Check to_make_prim_list encoding
prim_list = box_geo.to_make_prim_list(t=0.0)
assert len(prim_list) == 1
_, k_pen, vel, r_c_wall, omega_enc, r_ref_enc = prim_list[0]
omega_ok = abs(omega_enc - omega_box) < 1e-14
vel_ok   = np.allclose(vel, [vx_box, vy_box], atol=1e-14)
r_ref_ok = np.allclose(r_ref_enc, [box_cx, box_cy], atol=1e-14)
print(f"\n  omega enc={omega_enc}  (exp {omega_box})  {'PASS' if omega_ok else 'FAIL'}")
print(f"  vel   enc={vel}  (exp [{vx_box},{vy_box}])  {'PASS' if vel_ok else 'FAIL'}")
print(f"  r_ref enc={r_ref_enc}  (exp [{box_cx},{box_cy}])  {'PASS' if r_ref_ok else 'FAIL'}")

# ── Part 2: physics simulation ─────────────────────────────────────────────────
print("\n--- Part 2: physics (run + frame recording) ---")

sys_b2 = System(Lx, Ly, periodic_x=True)

channel = Channel(width=Lx, height=Ly, x0=Lx/2, y0=Ly/2, exclusion='exterior')
channel.set_render(color='#333333', linewidth=2.5, alpha=0.9)

box = SquareObstacle(side=box_side, x0=box_cx, y0=box_cy, exclusion='interior')
box.set_motion(MotionSpec(omega=omega_box, vx=vx_box, vy=vy_box,
                          r_ref=(box_cx, box_cy)))
box.set_render(color='#aa4444', linewidth=2.0, alpha=0.85, fill=True)

sys_b2.add_object(channel)
sys_b2.add_object(box)

spec = ParticleSpec(count=25, nu=0.2, N_nodes=48, poly_dist=0.05)
sys_b2.add_particles(spec)

sys_b2.initialize(phi_target=0.78, seed=7, verbose=True, n_relax=200,
                  swell_alpha=10.0)

dt = sys_b2._dt
# Cap at 2000 steps (~10 min wall time at ~235 ms/step for 25 particles)
n_steps      = 2000
sample_every = max(1, n_steps // 80)   # ~80 frames
print(f"  dt={dt:.6f}, n_steps={n_steps}, sample_every={sample_every}")
print(f"  Box: ω={omega_box} rad/s, v=({vx_box},{vy_box}) units/s")

sys_b2.run(n_steps, sample_every=sample_every, verbose=True)

t_final  = sys_b2.t
print(f"  Simulation ended at t={t_final:.4f} s")
print(f"  Box angle = {np.degrees(omega_box * t_final):.1f}°")
print(f"  Box centre = ({box_cx + vx_box*t_final:.2f}, {box_cy + vy_box*t_final:.2f})")

# Verify: all CMs outside the obstacle at every recorded frame
n_fail_total = 0
for fr in sys_b2.frames:
    t_fr   = fr['t']
    x_cm   = fr['x_cm']
    verts  = box.region_polygon(t_fr)['vertices']
    for i, cm in enumerate(x_cm):
        if _point_in_polygon(cm, verts):
            n_fail_total += 1
            print(f"  FAIL: particle {i} inside box at t={t_fr:.3f}")

phys_ok = (n_fail_total == 0)
if phys_ok:
    print(f"  PASS: all CMs outside obstacle at all {len(sys_b2.frames)} recorded frames")
else:
    print(f"  FAIL: {n_fail_total} violations across frames")

# ── Movie ──────────────────────────────────────────────────────────────────────
print("\nRendering movie…")
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Polygon as MplPolygon

colors = sys_b2._particle_colors
r_c_arr = sys_b2._params['r_c_per_p'].numpy()

# Use actual (post-swell) box dimensions for all rendering
rLx = sys_b2.Lx
rLy = sys_b2.Ly

fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.set_xlim(-0.5, rLx + 0.5)
ax.set_ylim(-0.5, rLy + 0.5)
ax.set_aspect('equal')
ax.set_xlabel('x / R₀')
ax.set_ylabel('y / R₀')

# Channel walls (static)
for d in channel.resolved(0.0):
    p = d['prim']
    ax.plot([p.p0[0], p.p1[0]], [p.p0[1], p.p1[1]],
            color='#333333', lw=2.5, zorder=1)

# Spinning + translating box — filled polygon, updated each frame
init_verts = box.region_polygon(0.0)['vertices']
box_patch  = MplPolygon(init_verts, closed=True,
                        facecolor='#ffcccc', edgecolor='#aa4444', lw=2.0,
                        alpha=0.85, zorder=2)
ax.add_patch(box_patch)

# Particle patches — 3 copies per particle: primary + ±rLx ghosts (periodic x)
from src.epd.system import System as _Sys

P = len(sys_b2._particles)
x_offsets = [-rLx, 0.0, rLx]  # use actual post-swell Lx
margin = 1.5                   # only draw copies within this margin of box edges

# p_patches[pi][k] = MplPolygon for particle pi at x_offset x_offsets[k]
p_patches = []
for pi in range(P):
    c   = colors[pi] if pi < len(colors) else 'cornflowerblue'
    row = []
    for _ in x_offsets:
        pat = MplPolygon(np.zeros((4, 2)), closed=True,
                         fc=c, ec='k', lw=0.8, alpha=0.65, zorder=3,
                         visible=False)
        ax.add_patch(pat)
        row.append(pat)
    p_patches.append(row)

all_patches = [pat for row in p_patches for pat in row]

time_txt = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                   fontsize=8, va='top', family='monospace')

def _update(i):
    fr    = sys_b2.frames[i]
    t_now = fr['t']

    # Update box
    verts_now = box.region_polygon(t_now)['vertices']
    box_patch.set_xy(verts_now)

    # Update particles with periodic ghost images
    for pi in range(P):
        xy  = fr['x_all'][pi]   # (N, 2) node positions (CMs already wrapped)
        cm  = fr['x_cm'][pi]    # (2,) centre of mass
        r_c = float(r_c_arr[pi])
        outer = _Sys._outer_contour(xy, cm, r_c)   # (N, 2) outer contour

        for k, odx in enumerate(x_offsets):
            shifted = outer + np.array([odx, 0.0])
            # Visibility: only draw if shifted contour overlaps the canvas
            in_canvas = (shifted[:, 0].max() > -margin and
                         shifted[:, 0].min() < rLx + margin)
            pat = p_patches[pi][k]
            if in_canvas:
                pat.set_xy(shifted)
                pat.set_visible(True)
            else:
                pat.set_visible(False)

    theta_now = omega_box * t_now
    cx_now    = box_cx + vx_box * t_now
    cy_now    = box_cy + vy_box * t_now
    time_txt.set_text(
        f't={t_now:.2f} s\n'
        f'θ={np.degrees(theta_now):.0f}°\n'
        f'ctr=({cx_now:.1f},{cy_now:.1f})'
    )
    return all_patches + [box_patch, time_txt]

ani = animation.FuncAnimation(fig, _update,
                               frames=len(sys_b2.frames),
                               interval=80, blit=True)
ani.save(args.out, writer='pillow', fps=12)
plt.close(fig)
print(f"Movie saved → {args.out}")

# ── Summary ────────────────────────────────────────────────────────────────────
all_pass = geo_ok and omega_ok and vel_ok and r_ref_ok and phys_ok
print(f"\n{'PASS' if all_pass else 'FAIL'}: Test B2 {'passed' if all_pass else 'FAILED'}")
if not all_pass:
    sys.exit(1)
