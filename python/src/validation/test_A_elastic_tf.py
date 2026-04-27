"""
test_A_elastic_tf.py — Test A: Two-disk oscillatory squeeze via C++/TF backend

Geometry and rendering are identical to the approved test_A_twodisk_oscillate.py.
The integration loop uses step_full_tf + CandidacyManager (C++ pybind11) instead
of the NumPy per-node Euler loop.

Usage:
    python src/validation/test_A_elastic_tf.py [--cycles N] [--out PATH]
"""

import sys, os, argparse, time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Polygon as MplPolygon
import matplotlib.animation as animation
import tensorflow as tf

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

# ── TF backend setup ──────────────────────────────────────────────────────────
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)
from src.simulation.tf_sim import (make_state, make_prim_data, step_full_tf,
                                    DTYPE, NP_DTYPE)
from src.simulation.candidacy_manager import CandidacyManager

# ── particle / contact primitives ─────────────────────────────────────────────
from src.epd.particles import ParticleSpec
from src.simulation.capsule_shell import CapsuleParticle
from src.simulation.contact_primitives import LineSegment


# ── geometry helpers (identical to NumPy test) ────────────────────────────────

def _capsule_outline_polygon(x, r_c, n_arc=6):
    N  = len(x)
    xn = np.roll(x, -1, axis=0)
    e  = xn - x
    L  = np.linalg.norm(e, axis=1, keepdims=True)
    t  = e / np.where(L > 1e-15, L, 1.0)
    n_edge = np.column_stack([t[:, 1], -t[:, 0]])

    pts = []
    for i in range(N):
        i_prev = (i - 1) % N
        a1 = np.arctan2(n_edge[i_prev, 1], n_edge[i_prev, 0])
        a2 = np.arctan2(n_edge[i,      1], n_edge[i,      0])
        da = (a2 - a1) % (2 * np.pi)
        if da < 1e-10:
            da = 2 * np.pi
        if da <= np.pi:
            n_pts  = max(2, int(n_arc * da / np.pi) + 1)
            thetas = np.linspace(a1, a1 + da, n_pts)
            arc    = x[i] + r_c * np.column_stack([np.cos(thetas), np.sin(thetas)])
        else:
            arc = x[i] + r_c * np.vstack([n_edge[i_prev], n_edge[i]])
        pts.append(arc)

    outline = np.vstack(pts)
    return np.vstack([outline, outline[:1]])


class OscillatingWall:
    def __init__(self, is_top, y0, half_w, A, omega, r_c_wall):
        self.is_top   = is_top
        self.y0       = float(y0)
        self.half_w   = float(half_w)
        self.A        = float(A)
        self.omega    = float(omega)
        self.r_c_wall = float(r_c_wall)
        self.sign     = -1.0 if is_top else +1.0
        inward_n = [0, -1] if is_top else [0, +1]
        self.seg = LineSegment([-half_w, y0], [half_w, y0],
                               inward_normal=inward_n)
        self.seg.r_c = r_c_wall

    def update(self, t):
        y = self.y0 + self.sign * self.A * np.sin(self.omega * t)
        self.seg.p0[1] = y
        self.seg.p1[1] = y
        return self.seg

    @property
    def y(self):
        return float(self.seg.p0[1])

    def wall_strain(self):
        return abs(self.y - self.y0) / 1.0   # normalised to R0=1


# ── main ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument('--cycles', type=int, default=3)
parser.add_argument('--out',    type=str, default='results/test_A_elastic_tf.gif')
args = parser.parse_args()

R0    = 1.0
N     = 32
nu    = 0.5

spec = ParticleSpec(count=2, nu=nu, N_nodes=N)
d    = spec.derived
tau  = d['tau']
C    = d['C']
K_area = d['K_area']
S    = 1.0

print("=" * 60)
print(f"Test A (elastic TF): Two-disk oscillatory squeeze  (ν={nu})")
print("=" * 60)
print(f"  q={d['q']:.4f}  TAU={d['TAU']:.4f}")
print(f"  tau={tau:.4f}  El_t={d['El_t']:.4f}  K_area={K_area:.4f}")
print(f"  C={C:.2f}  alpha={d['alpha']:.4f}")

# Reference particle for dt / T_wave
p_ref = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area)
L0    = p_ref.L0
r_c   = p_ref.r_c

El_t    = d['El_t']
c_edge  = np.sqrt(El_t / (p_ref.rho_d * L0))
T_wave  = 2.0 * np.pi * R0 / c_edge
alpha_damp = 2.0 / T_wave

# Oscillation: quasi-static (T_osc = 5 × T_wave)
T_osc    = 5.0 * T_wave
omega_osc = 2.0 * np.pi / T_osc
A_wall   = 0.1 * R0

# dt from wave-speed stability
from src.simulation.capsule_shell import CapsuleSim
sim_probe = CapsuleSim([p_ref])
dt_max, _ = sim_probe.estimate_dt_max()
dt = 0.4 * dt_max

n_cycles = args.cycles
t_total  = n_cycles * T_osc
n_steps  = int(np.ceil(t_total / dt))

print(f"\n  Geometry:  R0={R0}  N={N}  L0={L0:.4f}  r_c={r_c:.4f}")
print(f"  Dynamics:  dt={dt:.5f}  T_wave={T_wave:.4f}  alpha={alpha_damp:.4f}")
print(f"  Oscillation: A={A_wall}  T_osc={T_osc:.4f}  omega={omega_osc:.4f}")
print(f"  Total: {n_cycles} cycles  {n_steps} steps  t_total={t_total:.3f}")

# ── initial positions (identical to NumPy test) ───────────────────────────────
cy_disk = R0 + L0   # = 1 + L0 ≈ 1.196
p1 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0, -cy_disk))
p2 = CapsuleParticle(N=N, R0=R0, tau=tau, S=S, C=C, K_area=K_area,
                     center=(0.0,  cy_disk))
A0 = p1.A0

# Walls (same as NumPy test)
half_w    = R0 + 2.0 * L0
y_top0    =  2.0 * R0 + 3.0 * L0
y_bot0    = -(2.0 * R0 + 3.0 * L0)
r_c_wall  = L0

left_wall  = LineSegment([-half_w, -y_top0*2], [-half_w,  y_top0*2],
                         inward_normal=[1, 0])
left_wall.r_c  = r_c_wall
right_wall = LineSegment([ half_w, -y_top0*2], [ half_w,  y_top0*2],
                         inward_normal=[-1, 0])
right_wall.r_c = r_c_wall

top_wall = OscillatingWall(is_top=True,  y0=y_top0, half_w=half_w,
                            A=A_wall, omega=omega_osc, r_c_wall=r_c_wall)
bot_wall = OscillatingWall(is_top=False, y0=y_bot0, half_w=half_w,
                            A=A_wall, omega=omega_osc, r_c_wall=r_c_wall)

print(f"  Box: half_w={half_w:.4f}  y_top0={y_top0:.4f}")
cc_gap_init = np.linalg.norm(p2.x_cm - p1.x_cm) - 2.0 * (R0 + r_c)
print(f"  Initial cc_gap = {cc_gap_init:.4e}")

# ── TF state and params ───────────────────────────────────────────────────────
state, params = make_state([p1, p2])

# Non-periodic box: large dummy Lx/Ly (only used for candidacy skin, not wrapping)
Lx_cm, Ly_cm = 20.0, 20.0

# Candidacy manager
R0_arr = np.array([p1.R0, p2.R0])
cm_mgr = CandidacyManager(
    P=2, N=N,
    R0=float(np.mean(R0_arr)),
    E=128,
    skin=1.0 * float(np.mean(R0_arr)),
    periodic=False,
    periodic_x=False,
    periodic_y=False,
    Lx=Lx_cm, Ly=Ly_cm,
    R0_arr=R0_arr,
)
cm_mgr.update(state['x_cm'].numpy(), state['theta'].numpy())

# TF scalars
dt_tf    = tf.constant(NP_DTYPE(dt),         dtype=DTYPE)
alpha_tf = tf.constant(NP_DTYPE(alpha_damp), dtype=DTYPE)
g_tf     = tf.constant(NP_DTYPE(0.0),        dtype=DTYPE)


def _make_prim(t_now):
    """Build prim_data from current wall positions."""
    prims = [left_wall, right_wall,
             top_wall.update(t_now), bot_wall.update(t_now)]
    return make_prim_data([(seg, 1.0, np.zeros(2), seg.r_c, 0.0, np.zeros(2))
                           for seg in prims])


# ── frame recording ───────────────────────────────────────────────────────────
n_frames_per_cycle = 20
n_frames   = n_cycles * n_frames_per_cycle
record_times = np.linspace(0.0, t_total, n_frames + 1)
frames       = []
record_idx   = 0


def _area_poly(x):
    """Shoelace area of (N,2) polygon."""
    xn = np.roll(x, -1, axis=0)
    return 0.5 * abs(np.sum(x[:, 0] * xn[:, 1] - xn[:, 0] * x[:, 1]))


def _vert_strain(x_nodes, R0):
    """Vertical compression strain from top/bottom node y positions."""
    N = len(x_nodes)
    return 1.0 - (x_nodes[N // 4, 1] - x_nodes[3 * N // 4, 1]) / (2.0 * R0)


def record_frame(t_now):
    x_all = state['x_all'].numpy()   # (2, N, 2)
    x1, x2 = x_all[0], x_all[1]
    eps1 = _vert_strain(x1, R0)
    eps2 = _vert_strain(x2, R0)
    A1   = _area_poly(x1)
    A2   = _area_poly(x2)
    x_cm = state['x_cm'].numpy()
    cc_gap = np.linalg.norm(x_cm[1] - x_cm[0]) - 2.0 * (R0 + r_c)
    frames.append({
        't': t_now,
        'x1': x1.copy(), 'x2': x2.copy(),
        'y_top': top_wall.y, 'y_bot': bot_wall.y,
        'half_w': half_w,
        'wall_strain': top_wall.wall_strain(),
        'eps1': eps1, 'eps2': eps2,
        'cc_gap': cc_gap,
        'A1': A1, 'A2': A2, 'A0': A0,
        'r_c': r_c, 'r_c_wall': r_c_wall,
    })


# record initial frame (walls at t=0, no step yet)
top_wall.update(0.0); bot_wall.update(0.0)
record_frame(0.0)
record_idx = 1

# ── step loop ─────────────────────────────────────────────────────────────────
t = 0.0
t_start = time.time()

for step_i in range(n_steps):
    # Update prim_data at current time (wall positions at t)
    prim_data = _make_prim(t)

    # Candidacy update
    x_cm_np  = state['x_cm'].numpy()
    theta_np = state['theta'].numpy()
    if cm_mgr.needs_update(x_cm_np, theta_np):
        cm_mgr.update(x_cm_np, theta_np)
    caps = tf.constant(cm_mgr.CapCandidates, dtype=tf.int32)

    # TF step
    t_tf   = tf.constant(NP_DTYPE(t), dtype=DTYPE)
    state, _ = step_full_tf(state, caps, dt_tf, alpha_tf, g_tf,
                             params, t=t_tf, prim_data=prim_data)

    t += dt

    if record_idx < len(record_times) and t >= record_times[record_idx]:
        top_wall.update(t); bot_wall.update(t)
        record_frame(t)
        record_idx += 1

    if (step_i + 1) % max(1, n_steps // 5) == 0:
        fr = frames[-1]
        elapsed = time.time() - t_start
        print(f"    step {step_i+1}/{n_steps}  t={t:.3f}  "
              f"wall_strain={fr['wall_strain']:.4f}  "
              f"eps_p={fr['eps1']:.4f}  "
              f"cc_gap={fr['cc_gap']:.4e}  ({elapsed:.0f}s)")

print(f"\n  Recorded {len(frames)} frames  ({time.time()-t_start:.1f}s total)")


# ── diagnostics ───────────────────────────────────────────────────────────────
print("\n  ── Diagnostics ──────────────────────────────")
max_wall_strain = max(f['wall_strain'] for f in frames)
max_eps         = max(0.5 * (f['eps1'] + f['eps2']) for f in frames)
max_dA          = max(abs(1.0 - 0.5 * (f['A1'] + f['A2']) / f['A0']) for f in frames)
min_cc_gap      = min(f['cc_gap'] for f in frames)

print(f"  peak wall_strain : {max_wall_strain:.4f}  (target ≈ {A_wall:.4f})")
print(f"  peak ε_particle  : {max_eps:.4f}")
print(f"  peak |1-A/A₀|   : {max_dA:.6f}  (area conservation)")
print(f"  min cc_gap       : {min_cc_gap:.4e}")

ok_strain  = abs(max_wall_strain - A_wall) < 0.005
ok_area    = max_dA < 0.05
ok_contact = min_cc_gap < 0.0

print(f"\n  [{'PASS' if ok_strain  else 'FAIL'}] wall_strain ≈ A_wall  ({max_wall_strain:.4f} vs {A_wall})")
print(f"  [{'PASS' if ok_area   else 'FAIL'}] area conservation < 5%  (actual {max_dA:.4f})")
print(f"  [{'PASS' if ok_contact else 'FAIL'}] disks make contact  (min_cc_gap = {min_cc_gap:.4e})")


# ── render GIF (identical to NumPy test) ─────────────────────────────────────
import pathlib
out_path = pathlib.Path(args.out)
out_path.parent.mkdir(parents=True, exist_ok=True)

r_c_wall_f = frames[0]['r_c_wall']
half_w_f   = frames[0]['half_w']
A0_f       = frames[0]['A0']

y_top_max  = max(abs(f['y_top']) for f in frames)
lim_y      = y_top_max + r_c_wall_f + 0.15
lim_x      = half_w_f  + r_c_wall_f + 0.15

fig, axes = plt.subplots(1, 2, figsize=(10, 5.5))
ax_sim = axes[0]
ax_sim.set_aspect('equal')
ax_sim.set_xlim(-lim_x, lim_x)
ax_sim.set_ylim(-lim_y, lim_y)
ax_sim.set_xlabel('x / R₀')
ax_sim.set_ylabel('y / R₀')
ax_sim.set_title('Two-disk oscillatory squeeze  (ν=0.5, 10% strain)  [TF]')

wall_kw = dict(fc='#888888', ec='#555555', lw=0.5, alpha=0.85, zorder=1)

ax_sim.add_patch(mpatches.Rectangle(
    (-lim_x, -lim_y * 1.5), (lim_x - half_w_f + r_c_wall_f), lim_y * 3.0, **wall_kw))
ax_sim.add_patch(mpatches.Rectangle(
    (half_w_f - r_c_wall_f, -lim_y * 1.5), (lim_x - half_w_f + r_c_wall_f), lim_y * 3.0, **wall_kw))

top_patch = mpatches.Rectangle((-lim_x, 0.0), 2*lim_x, lim_y, **wall_kw)
bot_patch = mpatches.Rectangle((-lim_x, -lim_y), 2*lim_x, lim_y, **wall_kw)
ax_sim.add_patch(top_patch)
ax_sim.add_patch(bot_patch)

fill1 = MplPolygon(np.zeros((4, 2)), closed=True,
                    fc='cornflowerblue', ec='navy', lw=1.0, alpha=0.6, zorder=3)
fill2 = MplPolygon(np.zeros((4, 2)), closed=True,
                    fc='salmon', ec='darkred', lw=1.0, alpha=0.6, zorder=3)
ax_sim.add_patch(fill1)
ax_sim.add_patch(fill2)

txt = ax_sim.text(0.02, 0.98, '', transform=ax_sim.transAxes,
                  fontsize=8, va='top', family='monospace')

ax_ts = axes[1]
ax_ts.set_xlabel('t / T_osc')
ax_ts.set_ylabel('Value')
ax_ts.set_title('Time series')

t_arr      = np.array([f['t'] for f in frames])
strain_arr = np.array([f['wall_strain'] for f in frames])
eps_arr    = np.array([0.5 * (f['eps1'] + f['eps2']) for f in frames])

l_strain, = ax_ts.plot([], [], 'b-',  lw=1.5, label='wall strain')
l_eps,    = ax_ts.plot([], [], 'r-',  lw=1.5, label='ε_particle (avg)')
ax_ts.legend(fontsize=7)
ax_ts.set_xlim(0, t_arr[-1] / T_osc)
y_max = max(strain_arr.max(), eps_arr.max(), 0.02) * 1.3
ax_ts.set_ylim(-0.005, y_max)
ax_ts.axhline(0, color='k', lw=0.5)
vline = ax_ts.axvline(0, color='gray', lw=0.8, ls='--')

plt.tight_layout()


def _update(i):
    fr = frames[i]
    y_t = fr['y_top']
    y_b = fr['y_bot']

    top_patch.set_xy((-lim_x, y_t - r_c_wall_f))
    top_patch.set_height(lim_y - (y_t - r_c_wall_f))
    bot_patch.set_xy((-lim_x, -lim_y))
    bot_patch.set_height((y_b + r_c_wall_f) - (-lim_y))

    fill1.set_xy(_capsule_outline_polygon(fr['x1'], fr['r_c']))
    fill2.set_xy(_capsule_outline_polygon(fr['x2'], fr['r_c']))

    t_norm = fr['t'] / T_osc
    l_strain.set_data(t_arr[:i+1] / T_osc, strain_arr[:i+1])
    l_eps.set_data(t_arr[:i+1]    / T_osc, eps_arr[:i+1])
    vline.set_xdata([t_norm])

    eps_avg = 0.5 * (fr['eps1'] + fr['eps2'])
    txt.set_text(
        f"t/T={t_norm:.2f}  ε_wall={fr['wall_strain']:.3f}\n"
        f"ε_p={eps_avg:.3f}\n"
        f"cc_gap={fr['cc_gap']:+.4f}"
    )
    return [top_patch, bot_patch, fill1, fill2, l_strain, l_eps, vline, txt]


ani = animation.FuncAnimation(fig, _update, frames=len(frames),
                              interval=1000 // 10, blit=True)
ani.save(str(out_path), writer='pillow', fps=10)
plt.close(fig)
print(f"\n  GIF saved → {out_path}")

if ok_strain and ok_area and ok_contact:
    print(f"\nPASS: Test A (elastic, TF)")
else:
    print(f"\nFAIL: Test A (elastic, TF)")
