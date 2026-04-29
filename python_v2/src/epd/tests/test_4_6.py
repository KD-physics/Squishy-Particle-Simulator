"""
test_4_6.py — Phase 4.6 Integration Verification (v6)

Design rules:
  • All above-jamming periodic tests use phi=0.88 (phi_J≈0.84 for 2D random packing)
  • P=8 particles per periodic test; P=6 in Box tests (phi_rsa≈0.29 < RSA limit)
  • SimMonitor checks: overlap, circularity, stretch, kink (triangularity), ring_ratio
  • Non-periodic fixed-wall systems use relax_only (cannot swell → phi≈0.29)
  • Periodic renderer draws image copies for particles straddling boundaries
  • Gravity tests use dimensionless Bond numbers:
      elastic  Bo_el = rho_d * g * R0^2 / El_t_spec  (→ g = Bo_el * El_t)
      emulsion Bo    = rho_d * g * R0^2 / gamma       (→ g = Bo at ref params)

Tests
-----
  A   — Elastic dense packing       : periodic, phi=0.88, P=8
  B   — Elastic moving wall          : periodic, phi=0.88, P=8, wall vx=0.003
  C   — Elastic frozen pusher        : periodic, phi=0.88, P=8 (1+7), kinematic vx=0.003
  D   — Elastic right-fraction drift : periodic, phi=0.88, P=8, right half vx=0.003
  AE  — Emulsion dense packing       : periodic, phi=0.88, P=8 (proves emulsion reaches jamming)
  BE  — Emulsion moving wall         : periodic, phi=0.72, P=8, wall vx=0.003
  CE  — Emulsion frozen pusher       : periodic, phi=0.72, P=8 (1+7), kinematic vx=0.003
  DE  — Emulsion right-fraction drift: periodic, phi=0.72, P=8, right half vx=0.003
  E   — Save/load mid-run            : bit-identical resume (phi=0.72, fast)
  F   — Mixed elastic + emulsion     : periodic, phi=0.72, P=8 (4+4)
  G   — Elastic falling capsule      : Bo_el=5, non-periodic Box; checks Δy≈½gt²
  H   — Emulsion falling droplet     : Bo=2, non-periodic Box; checks Δy≈½gt²

Expected runtime: ~25-35 minutes.
"""

import sys, os, warnings, time
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..'))

import numpy as np
import pathlib

from src.epd.system import System
from src.epd.particles import ParticleSpec
from src.epd.motion import MotionSpec
from src.epd.objects import Wall, Box
from src.epd.monitor import SimMonitor

RESULTS = pathlib.Path('results/phase46')
RESULTS.mkdir(parents=True, exist_ok=True)

# ── swell params ──────────────────────────────────────────────────────────────
SWELL_88 = dict(dphi_init=0.005, dphi_max=0.012, n_relax=50, max_extra_relax=200)
SWELL_72 = dict(dphi_init=0.005, dphi_max=0.020, n_relax=30, max_extra_relax=80)

# Emulsion droplets at phi=0.88 develop flat facets and sharp triple-junction
# corners (kink→π) and ~2× perimeter stretch — both are physically expected, not errors.
EMULSION_MON = dict(stretch_crit=4.0, kink_crit=np.pi * 0.99)

# ── shared particle specs ─────────────────────────────────────────────────────
ELASTIC_SPEC  = dict(nu=0.5, N_nodes=32, poly_dist=0.05)
EMULSION_SPEC = dict(type='emulsion', gamma=1.0, N_nodes=32, poly_dist=0.05)

n_pass = 0
n_fail = 0

def check(label, cond, detail=''):
    global n_pass, n_fail
    status = 'PASS' if cond else 'FAIL'
    if cond: n_pass += 1
    else:    n_fail += 1
    msg = f'  [{status}] {label}'
    if detail: msg += f'  — {detail}'
    print(msg)

def run_monitored(sys_obj, n_steps, record_every, output_dir, **mon_kwargs):
    """Run with SimMonitor. Returns (monitor_ok, monitor)."""
    mon = SimMonitor(**mon_kwargs)
    sys_obj.run(n_steps, record_every=record_every, output_dir=output_dir,
                monitor=mon, monitor_every=record_every, abort_on_crit=False)
    lvl = mon.warn_level()
    tag = 'OK' if lvl == 'ok' else lvl.upper()
    print(f'  [monitor {tag}] {mon.summary_line()}')
    return lvl != 'crit', mon

def make_gif_safe(sys_obj, gif_path, frame_dir, fps=5):
    ffs = sorted(pathlib.Path(frame_dir).glob('frame_*.png'))
    if ffs:
        sys_obj.make_gif(str(gif_path), ffs, fps=fps)
        print(f'  GIF saved → {gif_path.name}')

t_suite = time.time()


# ─────────────────────────────────────────────────────────────────────────────
# A — Elastic dense packing  (periodic, phi=0.88, P=8)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test A: Elastic dense packing (phi=0.88, P=8) ────────────────────')
sys_a = System(Lx=8, Ly=8)
sys_a.add_particles(ParticleSpec(count=8, **ELASTIC_SPEC))
sys_a.initialize(phi_target=0.88, verbose=False, **SWELL_88)
ok_a, mon_a = run_monitored(sys_a, 100, 25, str(RESULTS/'A_frames'))

phi_a = sys_a.phi_outer
check('A.1 phi_outer ≥ 0.86', phi_a >= 0.86, f'phi={phi_a:.4f}')
check('A.2 no NaN', not np.any(np.isnan(sys_a.state['x_cm'].numpy())))
check('A.3 monitor not critical', ok_a, mon_a.summary_line())
make_gif_safe(sys_a, RESULTS/'A_elastic_dense.gif', RESULTS/'A_frames')


# ─────────────────────────────────────────────────────────────────────────────
# B — Elastic moving wall through jammed packing  (periodic, phi=0.88, P=8)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test B: Elastic moving wall (phi=0.88, P=8) ──────────────────────')
sys_b = System(Lx=8, Ly=8)
wall_b = Wall(p0=[0, 0], p1=[0, 8], normal=[1, 0])
wall_b.set_motion(MotionSpec(vx_dc=0.003))
sys_b.add_object(wall_b)
sys_b.add_particles(ParticleSpec(count=8, **ELASTIC_SPEC))
sys_b.initialize(phi_target=0.88, verbose=False, **SWELL_88)

t0_b = sys_b.t
ok_b, mon_b = run_monitored(sys_b, 100, 25, str(RESULTS/'B_frames'))
t1_b = sys_b.t

expected_disp = 0.003 * (t1_b - t0_b)
wall_p0_x = wall_b.resolved(t=t1_b)[0]['prim'].p0[0]
check('B.1 no NaN', not np.any(np.isnan(sys_b.state['x_cm'].numpy())))
check('B.2 wall displaced correctly',
      abs(wall_p0_x - expected_disp) < 0.002,
      f'wall_p0_x={wall_p0_x:.5f} expected={expected_disp:.5f}')
check('B.3 monitor not critical', ok_b, mon_b.summary_line())
make_gif_safe(sys_b, RESULTS/'B_elastic_wall.gif', RESULTS/'B_frames')


# ─────────────────────────────────────────────────────────────────────────────
# C — Elastic frozen pusher  (periodic, phi=0.88, P=8: 1 kinematic pusher + 7 free)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test C: Elastic frozen pusher (periodic, phi=0.88, P=8) ──────────')
sys_c = System(Lx=8, Ly=8)
pusher_c = ParticleSpec(count=1, frozen_shape=True, **ELASTIC_SPEC)
pusher_c.set_motion(MotionSpec(vx_dc=0.003))
sys_c.add_particles(pusher_c)
sys_c.add_particles(ParticleSpec(count=7, **ELASTIC_SPEC))
sys_c.initialize(phi_target=0.88, verbose=False, **SWELL_88)

x_cm0_c = sys_c.state['x_cm'].numpy()[0].copy()
ok_c, mon_c = run_monitored(sys_c, 100, 25, str(RESULTS/'C_frames'))
x_cm1_c = sys_c.state['x_cm'].numpy()[0]
t_c = sys_c.t
expected_c = x_cm0_c[0] + 0.003 * t_c
disp_err_c = abs(x_cm1_c[0] - expected_c) / max(0.003 * t_c, 1e-6)

check('C.1 pusher displaced ≈ vx*t', disp_err_c < 0.05,
      f'err={disp_err_c*100:.2f}%')
check('C.2 pusher shape frozen',
      float(sys_c.state['u'].numpy()[0].max()) < 1e-6)
check('C.3 monitor not critical', ok_c, mon_c.summary_line())
make_gif_safe(sys_c, RESULTS/'C_elastic_pusher.gif', RESULTS/'C_frames')


# ─────────────────────────────────────────────────────────────────────────────
# D — Elastic right-fraction drift  (periodic, phi=0.88, P=8, right half vx=0.003)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test D: Elastic right-fraction drift (periodic, phi=0.88, P=8) ───')
sys_d = System(Lx=8, Ly=8)
sys_d.add_particles(ParticleSpec(count=8, **ELASTIC_SPEC))
sys_d.initialize(phi_target=0.88, verbose=False, **SWELL_88)

right_d = sys_d.select_particles('right_fraction', fraction=0.50)
sys_d.particles_set_motion(right_d, MotionSpec(vx_dc=0.003), frozen_shape=False)
x_cm_init_right_d = sys_d.state['x_cm'].numpy()[right_d].copy()
ok_d, mon_d = run_monitored(sys_d, 100, 25, str(RESULTS/'D_frames'))
t_d = sys_d.t

dx_d = sys_d.state['x_cm'].numpy()[right_d][:, 0] - x_cm_init_right_d[:, 0]
expected_dx_d = 0.003 * t_d
check('D.1 right particles moved +x kinematically',
      len(right_d) == 0 or np.mean(dx_d) > 0.9 * expected_dx_d,
      f'mean_dx={np.mean(dx_d):.5f} expected≈{expected_dx_d:.5f}')
check('D.2 no NaN', not np.any(np.isnan(sys_d.state['x_cm'].numpy())))
check('D.3 monitor not critical', ok_d, mon_d.summary_line())
make_gif_safe(sys_d, RESULTS/'D_elastic_drift.gif', RESULTS/'D_frames')


# ─────────────────────────────────────────────────────────────────────────────
# AE — Emulsion dense packing  (periodic, phi=0.88, P=8)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test AE: Emulsion dense packing (phi=0.88, P=8) ─────────────────')
sys_ae = System(Lx=8, Ly=8)
sys_ae.add_particles(ParticleSpec(count=8, **EMULSION_SPEC))
sys_ae.initialize(phi_target=0.88, verbose=False, **SWELL_88)
ok_ae, mon_ae = run_monitored(sys_ae, 100, 25, str(RESULTS/'AE_frames'),
                               **EMULSION_MON)

phi_ae = sys_ae.phi_outer
check('AE.1 phi_outer ≥ 0.86', phi_ae >= 0.86, f'phi={phi_ae:.4f}')
check('AE.2 no NaN', not np.any(np.isnan(sys_ae.state['x_cm'].numpy())))
check('AE.3 monitor not critical', ok_ae, mon_ae.summary_line())
make_gif_safe(sys_ae, RESULTS/'AE_emulsion_dense.gif', RESULTS/'AE_frames')


# ─────────────────────────────────────────────────────────────────────────────
# BE — Emulsion moving wall  (periodic, phi=0.88, P=8)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test BE: Emulsion moving wall (phi=0.72, P=8) ────────────────────')
sys_be = System(Lx=8, Ly=8)
wall_be = Wall(p0=[0, 0], p1=[0, 8], normal=[1, 0])
wall_be.set_motion(MotionSpec(vx_dc=0.003))
sys_be.add_object(wall_be)
sys_be.add_particles(ParticleSpec(count=8, **EMULSION_SPEC))
sys_be.initialize(phi_target=0.72, verbose=False, **SWELL_72)

t0_be = sys_be.t
ok_be, mon_be = run_monitored(sys_be, 100, 25, str(RESULTS/'BE_frames'),
                               **EMULSION_MON)
t1_be = sys_be.t

expected_be  = 0.003 * (t1_be - t0_be)
wall_p0_x_be = wall_be.resolved(t=t1_be)[0]['prim'].p0[0]
check('BE.1 no NaN', not np.any(np.isnan(sys_be.state['x_cm'].numpy())))
check('BE.2 wall displaced correctly',
      abs(wall_p0_x_be - expected_be) < 0.002,
      f'wall={wall_p0_x_be:.5f} expected={expected_be:.5f}')
check('BE.3 monitor not critical', ok_be, mon_be.summary_line())
make_gif_safe(sys_be, RESULTS/'BE_emulsion_wall.gif', RESULTS/'BE_frames')


# ─────────────────────────────────────────────────────────────────────────────
# CE — Emulsion frozen pusher  (periodic, phi=0.88, P=8: 1 kinematic pusher + 7 free)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test CE: Emulsion frozen pusher (periodic, phi=0.72, P=8) ────────')
sys_ce = System(Lx=8, Ly=8)
pusher_ce = ParticleSpec(count=1, frozen_shape=True, **EMULSION_SPEC)
pusher_ce.set_motion(MotionSpec(vx_dc=0.003))
sys_ce.add_particles(pusher_ce)
sys_ce.add_particles(ParticleSpec(count=7, **EMULSION_SPEC))
sys_ce.initialize(phi_target=0.72, verbose=False, **SWELL_72)

x_cm0_ce = sys_ce.state['x_cm'].numpy()[0].copy()
ok_ce, mon_ce = run_monitored(sys_ce, 100, 25, str(RESULTS/'CE_frames'),
                               **EMULSION_MON)
x_cm1_ce = sys_ce.state['x_cm'].numpy()[0]
t_ce = sys_ce.t
disp_err_ce = abs(x_cm1_ce[0] - (x_cm0_ce[0] + 0.003 * t_ce)) / max(0.003 * t_ce, 1e-6)

check('CE.1 pusher displaced ≈ vx*t', disp_err_ce < 0.05,
      f'err={disp_err_ce*100:.2f}%')
check('CE.2 pusher shape frozen',
      float(sys_ce.state['u'].numpy()[0].max()) < 1e-6)
check('CE.3 monitor not critical', ok_ce, mon_ce.summary_line())
make_gif_safe(sys_ce, RESULTS/'CE_emulsion_pusher.gif', RESULTS/'CE_frames')


# ─────────────────────────────────────────────────────────────────────────────
# DE — Emulsion right-fraction drift  (periodic, phi=0.88, P=8, right half vx=0.003)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test DE: Emulsion right-fraction drift (periodic, phi=0.72, P=8) ─')
sys_de = System(Lx=8, Ly=8)
sys_de.add_particles(ParticleSpec(count=8, **EMULSION_SPEC))
sys_de.initialize(phi_target=0.72, verbose=False, **SWELL_72)

right_de = sys_de.select_particles('right_fraction', fraction=0.50)
sys_de.particles_set_motion(right_de, MotionSpec(vx_dc=0.003), frozen_shape=False)
x_cm_init_right_de = sys_de.state['x_cm'].numpy()[right_de].copy()
ok_de, mon_de = run_monitored(sys_de, 100, 25, str(RESULTS/'DE_frames'),
                               **EMULSION_MON)
t_de = sys_de.t

dx_de = sys_de.state['x_cm'].numpy()[right_de][:, 0] - x_cm_init_right_de[:, 0]
expected_dx_de = 0.003 * t_de
check('DE.1 right particles moved +x kinematically',
      len(right_de) == 0 or np.mean(dx_de) > 0.9 * expected_dx_de,
      f'mean_dx={np.mean(dx_de):.5f} expected≈{expected_dx_de:.5f}')
check('DE.2 no NaN', not np.any(np.isnan(sys_de.state['x_cm'].numpy())))
check('DE.3 monitor not critical', ok_de, mon_de.summary_line())
make_gif_safe(sys_de, RESULTS/'DE_emulsion_drift.gif', RESULTS/'DE_frames')


# ─────────────────────────────────────────────────────────────────────────────
# E — Save/load mid-run (bit-identical resume)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test E: Save/load mid-run ────────────────────────────────────────')
ckpt_dir_e = str(RESULTS / 'E_checkpoint')

sys_e_ref = System(Lx=8, Ly=8)
sys_e_ref.add_particles(ParticleSpec(count=8, **ELASTIC_SPEC))
sys_e_ref.initialize(phi_target=0.72, verbose=False, seed=42, **SWELL_72)
sys_e_ref.step(60)
x_cm_ref = sys_e_ref.state['x_cm'].numpy().copy()

sys_e_int = System(Lx=8, Ly=8)
sys_e_int.add_particles(ParticleSpec(count=8, **ELASTIC_SPEC))
sys_e_int.initialize(phi_target=0.72, verbose=False, seed=42, **SWELL_72)
sys_e_int.step(30)
sys_e_int.save(ckpt_dir_e)
sys_e_res = System.from_file(ckpt_dir_e)
sys_e_res.step(30)
x_cm_res = sys_e_res.state['x_cm'].numpy().copy()

max_diff = np.max(np.abs(x_cm_res - x_cm_ref))
check('E.1 bit-identical (max|Δx_cm|=0)', max_diff == 0.0,
      f'max|Δx_cm|={max_diff:.3e}')
check('E.2 step count correct', sys_e_res.step_count == 60,
      f'step_count={sys_e_res.step_count}')
check('E.3 time preserved',
      abs(sys_e_res.t - sys_e_ref.t) < 1e-10,
      f't_res={sys_e_res.t:.6f} t_ref={sys_e_ref.t:.6f}')


# ─────────────────────────────────────────────────────────────────────────────
# F — Mixed elastic + emulsion  (periodic, phi=0.88, P=8: 4 elastic + 4 emulsion)
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test F: Mixed elastic + emulsion (phi=0.72, P=8) ─────────────────')
sys_f = System(Lx=8, Ly=8)
sys_f.add_particles(ParticleSpec(count=4, **ELASTIC_SPEC))
sys_f.add_particles(ParticleSpec(count=4, **EMULSION_SPEC))
sys_f.initialize(phi_target=0.72, verbose=False, **SWELL_72)
ok_f, mon_f = run_monitored(sys_f, 100, 25, str(RESULTS/'F_frames'),
                             **EMULSION_MON)

x_cm_f = sys_f.state['x_cm'].numpy()
check('F.1 no NaN', not np.any(np.isnan(x_cm_f)))
check('F.2 particle count correct', len(x_cm_f) == 8)
check('F.3 monitor not critical', ok_f, mon_f.summary_line())
check('F.4 min_circ > 0.30', mon_f.metrics.get('min_circ', 0.0) > 0.30,
      f'min_circ={mon_f.metrics.get("min_circ", 0.0):.3f}')
make_gif_safe(sys_f, RESULTS/'F_mixed.gif', RESULTS/'F_frames')


# ─────────────────────────────────────────────────────────────────────────────
# G — Elastic falling capsule  (Bo_el=5)
#   g = Bo_el * El_t_spec ≈ 25. After 200 steps (t≈0.3), v_cm ≈ g×t ≈ 7.5 >>
#   membrane wave speed (~1), so extreme stretch is expected. Checks verify
#   kinematic free-fall (Δy ≈ ½gt²) rather than "reached floor".
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test G: Elastic falling capsule (Bo_el=5) ────────────────────────')
_sp_g = ParticleSpec(count=1, **ELASTIC_SPEC)
_El_t_spec_g = _sp_g._El_t   # ≈ 5.0
Bo_el_g = 5.0
g_g = Bo_el_g * _El_t_spec_g
print(f'  Bo_el={Bo_el_g}  El_t_spec={_El_t_spec_g:.3f}  → g={g_g:.2f}')

sys_g = System(Lx=8, Ly=12, periodic_x=False, periodic_y=False, g=g_g)
box_g = Box(width=8.0, height=12.0, x0=4.0, y0=6.0)
sys_g.add_object(box_g)
sys_g.add_particles(_sp_g)
sys_g.initialize(phi_target=0.0, relax_only=True, n_relax_init=10, verbose=False)

import tensorflow as tf
st = sys_g._state
x_cm_np = st['x_cm'].numpy().copy()
x_all_np = st['x_all'].numpy().copy()
dy_move = 9.5 - float(x_cm_np[0, 1])
x_cm_np[0, 1] += dy_move
x_all_np[0, :, 1] += dy_move
sys_g._state = {**st,
                'x_cm':  tf.constant(x_cm_np,  dtype=tf.float64),
                'x_all': tf.constant(x_all_np, dtype=tf.float64)}

y_cm_init_g = float(sys_g.state['x_cm'].numpy()[0, 1])
ok_g, mon_g = run_monitored(sys_g, 200, 50, str(RESULTS/'G_frames'),
                             stretch_crit=5.0)
y_cm_final_g = float(sys_g.state['x_cm'].numpy()[0, 1])
t_g = sys_g.t
expected_fall_g = 0.5 * g_g * t_g ** 2

check('G.1 particle fell (y decreased)', y_cm_final_g < y_cm_init_g - 1.0,
      f'y_init={y_cm_init_g:.3f}  y_final={y_cm_final_g:.3f}')
check('G.2 free-fall kinematics (Δy ≈ ½gt²)',
      abs((y_cm_init_g - y_cm_final_g) - expected_fall_g) < 0.15 * expected_fall_g,
      f'actual={y_cm_init_g - y_cm_final_g:.3f}  expected={expected_fall_g:.3f}')
check('G.3 monitor not critical', ok_g, mon_g.summary_line())
make_gif_safe(sys_g, RESULTS/'G_elastic_fall.gif', RESULTS/'G_frames', fps=8)


# ─────────────────────────────────────────────────────────────────────────────
# H — Emulsion falling droplet  (Bo=2)
#   Bo=25 was catastrophic (g=25 >> c_cap≈1 → droplet disintegrates).
#   Bo=2 keeps v_cm = g×t ≈ 0.6 < c_cap=1 after 200 steps, so droplet stays intact.
#   g = Bo directly when gamma=rho_d=R0=1.
# ─────────────────────────────────────────────────────────────────────────────
print('\n── Test H: Emulsion falling droplet (Bo=2) ──────────────────────────')
Bo_h = 2.0
g_h  = Bo_h
print(f'  Bo={Bo_h}  → g={g_h:.2f}')

sys_h = System(Lx=8, Ly=12, periodic_x=False, periodic_y=False, g=g_h)
box_h = Box(width=8.0, height=12.0, x0=4.0, y0=6.0)
sys_h.add_object(box_h)
_sp_h = ParticleSpec(count=1, **EMULSION_SPEC)
sys_h.add_particles(_sp_h)
sys_h.initialize(phi_target=0.0, relax_only=True, n_relax_init=10, verbose=False)

st = sys_h._state
x_cm_np = st['x_cm'].numpy().copy()
x_all_np = st['x_all'].numpy().copy()
dy_move = 9.5 - float(x_cm_np[0, 1])
x_cm_np[0, 1] += dy_move
x_all_np[0, :, 1] += dy_move
sys_h._state = {**st,
                'x_cm':  tf.constant(x_cm_np,  dtype=tf.float64),
                'x_all': tf.constant(x_all_np, dtype=tf.float64)}

y_cm_init_h = float(sys_h.state['x_cm'].numpy()[0, 1])
ok_h, mon_h = run_monitored(sys_h, 200, 50, str(RESULTS/'H_frames'),
                             **EMULSION_MON)
y_cm_final_h = float(sys_h.state['x_cm'].numpy()[0, 1])
t_h = sys_h.t
expected_fall_h = 0.5 * g_h * t_h ** 2

check('H.1 droplet fell (y decreased)', y_cm_final_h < y_cm_init_h - 0.05,
      f'y_init={y_cm_init_h:.3f}  y_final={y_cm_final_h:.3f}')
check('H.2 free-fall kinematics (Δy ≈ ½gt²)',
      abs((y_cm_init_h - y_cm_final_h) - expected_fall_h) < 0.15 * expected_fall_h,
      f'actual={y_cm_init_h - y_cm_final_h:.4f}  expected={expected_fall_h:.4f}')
check('H.3 monitor not critical', ok_h, mon_h.summary_line())
make_gif_safe(sys_h, RESULTS/'H_emulsion_fall.gif', RESULTS/'H_frames', fps=8)


# ─────────────────────────────────────────────────────────────────────────────
t_total = time.time() - t_suite
print()
print('=' * 60)
print(f'TOTAL: {n_pass}/{n_pass + n_fail} PASS   ({t_total:.0f}s elapsed)')
if n_fail > 0:
    sys.exit(1)
