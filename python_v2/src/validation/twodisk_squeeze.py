"""
twodisk_squeeze.py — Two-disk squeeze benchmark using v1 System API.

Two identical elastic disks are stacked vertically (touching at t=0) inside
top and bottom moving walls.  Both walls move inward via OscillatingWall
(quarter-period ramp), giving near-constant-velocity compression.  No side
walls → free lateral expansion → direct ν measurement.

Geometry (continuous-limit approximation, each capsule node circle has r_c=L0):
  disk centres :  y = ±(R0 + L0)
  top/bot walls:  y = ±(2·R0 + 3·L0)   ← wall capsule surface just touches disk
  Lx            :  2·(R0 + 2·L0)        ← enough lateral room (no side walls)

Public API
----------
run_squeeze(nu, N, R0=1.0, delta_max=0.12, n_frames=60,
            Oh=None, strain_rate_ratio=0.02, verbose=False)
    → (frames, spec)          ← high-level: ν as input, ParticleSpec handles params

run_squeeze_raw(q, tau, N, R0=1.0, delta_max=0.14, n_frames=80,
                C_factor=3000.0, alpha_damp=2.0, strain_rate_ratio=0.01,
                verbose=False)
    → frames                  ← low-level: raw (q, τ) input; used by calibration_sweep

interpolate_at_strain(frames, eps_ref)
    → (nu_meas, dA_frac)

Frame dict keys
---------------
t             : float   simulation time
wall_strain   : float   top-wall displacement / R0  (0 → delta_max at end)
eps_v         : float   mean vertical strain of the two disks
eps_l         : float   mean lateral strain
nu_meas       : float   ν = eps_l / eps_v  (nan if eps_v < 1e-4)
dA_frac       : float   mean fractional area change (A - A0) / A0
cc_gap_min    : float   most negative capsule-capsule gap seen so far
x_all         : (2, N, 2)  node positions of both disks
x_cm          : (2, 2)     CM positions
f_elastic     : (2, N, 2)  elastic forces
f_contact     : (2, N, 2)  contact forces
"""

import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.objects import OscillatingWall
from src.epd.system import System
from src.simulation.capsule_shell import CapsuleParticle


# ── geometry helpers ──────────────────────────────────────────────────────────

def _bbox(x_all):
    return x_all[:, 0].min(), x_all[:, 0].max(), x_all[:, 1].min(), x_all[:, 1].max()


def _particle_strains(x_all, R0):
    """Vertical and lateral strain of a single disk relative to undeformed radius."""
    xmn, xmx, ymn, ymx = _bbox(x_all)
    h_v = (ymx - ymn) / 2.0
    h_l = (xmx - xmn) / 2.0
    eps_v =  (R0 - h_v) / R0   # positive = compressed
    eps_l = -(R0 - h_l) / R0   # positive = expanded laterally
    return float(eps_v), float(eps_l)


def _polygon_area(x):
    """Shoelace area of a closed polygon (N,2)."""
    xr = np.roll(x, -1, axis=0)
    return 0.5 * abs((x[:, 0] * xr[:, 1] - xr[:, 0] * x[:, 1]).sum())


# ── shared simulation core ────────────────────────────────────────────────────

def _run_squeeze_core(d, N, R0, delta_max, n_frames, strain_rate_ratio, verbose):
    """
    Internal: run the two-disk squeeze given a pre-built param dict d with keys:
      El_t, K_area, C, tau, alpha  (all scalars)
    Returns frames list.
    """
    L0  = 2.0 * np.pi * R0 / N
    r_c = L0

    cy     =   R0 + L0
    y_top0 =   2.0 * R0 + 3.0 * L0
    y_bot0 = -(2.0 * R0 + 3.0 * L0)
    half_w =   R0 + 2.0 * L0
    Lx     = 2.0 * half_w
    Ly     = 2.0 * y_top0

    A_wall     = delta_max * R0
    c_wave     = np.sqrt(d['El_t'])
    v_wall     = strain_rate_ratio * c_wave
    t_squeeze  = A_wall / v_wall
    omega_wall = np.pi / (2.0 * t_squeeze)

    if verbose:
        print(f"  L0={L0:.4f}  r_c={r_c:.4f}  cy={cy:.4f}"
              f"  y_top0={y_top0:.4f}  A_wall={A_wall:.4f}")
        print(f"  c_wave={c_wave:.4f}  v_wall={v_wall:.4f}"
              f"  t_squeeze={t_squeeze:.3f}  omega_wall={omega_wall:.4f}")

    sys_ = System(Lx, Ly, periodic_x=False, periodic_y=False, g=0.0)
    top_wall = OscillatingWall(y0=y_top0, half_w=half_w,
                               A=A_wall, omega=omega_wall, is_top=True,  r_c_wall=L0)
    bot_wall = OscillatingWall(y0=y_bot0, half_w=half_w,
                               A=A_wall, omega=omega_wall, is_top=False, r_c_wall=L0)
    sys_.add_object(top_wall)
    sys_.add_object(bot_wall)

    p_bot = CapsuleParticle(N=N, R0=R0, tau=d['tau'], S=1.0,
                            C=d['C'], K_area=d['K_area'], center=[0.0, -cy])
    p_top = CapsuleParticle(N=N, R0=R0, tau=d['tau'], S=1.0,
                            C=d['C'], K_area=d['K_area'], center=[0.0,  cy])
    for p in [p_bot, p_top]:
        p.alpha_damp = d['alpha']

    sys_.initialize_from_particles([p_bot, p_top], alpha_damp=d['alpha'])

    x0    = sys_._state['x_all'].numpy()
    A0_arr = np.array([_polygon_area(x0[p]) for p in range(2)])

    frames       = []
    cc_gap_worst = [0.0]

    def _cb(sys_obj):
        snap  = sys_obj.snapshot()
        fb    = sys_obj.eval_forces()
        t_now = snap['t']
        x_all = snap['x_all']
        x_cm  = snap['x_cm']

        y_top_now   = top_wall.y_at(t_now)
        y_bot_now   = bot_wall.y_at(t_now)
        wall_strain = abs(y_top0 - y_top_now) / R0

        eps_v_arr, eps_l_arr = [], []
        for p in range(2):
            ev, el = _particle_strains(x_all[p], R0)
            eps_v_arr.append(ev); eps_l_arr.append(el)
        eps_v   = float(np.mean(eps_v_arr))
        eps_l   = float(np.mean(eps_l_arr))
        nu_meas = eps_l / eps_v if abs(eps_v) > 1e-4 else float('nan')

        A_arr   = np.array([_polygon_area(x_all[p]) for p in range(2)])
        dA_frac = float(np.mean((A_arr - A0_arr) / A0_arr))

        d12    = np.linalg.norm(x_all[0][:, None, :] - x_all[1][None, :, :], axis=-1)
        cc_gap = float(d12.min()) - 2.0 * r_c
        if cc_gap < cc_gap_worst[0]:
            cc_gap_worst[0] = cc_gap

        frames.append({
            't':           t_now,
            'wall_strain': wall_strain,
            'y_top':       float(y_top_now),
            'y_bot':       float(y_bot_now),
            'eps_v':       eps_v,
            'eps_l':       eps_l,
            'nu_meas':     nu_meas,
            'dA_frac':     dA_frac,
            'cc_gap_min':  cc_gap_worst[0],
            'x_all':       x_all.copy(),
            'x_cm':        x_cm.copy(),
            'f_elastic':   fb['f_elastic'].copy(),
            'f_contact':   fb['f_contact'].copy(),
        })
        return {}

    dt           = sys_._dt
    n_tot        = max(1, int(t_squeeze / dt))
    sample_every = max(1, n_tot // n_frames)

    if verbose:
        print(f"  n_steps={n_tot}  sample_every={sample_every}"
              f"  dt={dt:.2e}  n_frames≈{n_tot // sample_every}")

    sys_.run(n_tot, sample_every=sample_every,
             callback=_cb, record_initial=True, verbose=verbose)

    if verbose:
        last = frames[-1]
        print(f"  done: wall_strain={last['wall_strain']:.4f}"
              f"  ν_meas={last['nu_meas']:.4f}  dA={last['dA_frac']:.4f}")

    return frames


# ── public entry points ───────────────────────────────────────────────────────

def run_squeeze(nu, N, R0=1.0, delta_max=0.12, n_frames=60,
                Oh=None, strain_rate_ratio=0.02, verbose=False):
    """
    High-level entry point: ν is the input; all EPD params derived via ParticleSpec.

    Parameters
    ----------
    nu               : float  — target Poisson ratio
    N                : int    — number of perimeter nodes
    R0               : float  — particle radius
    delta_max        : float  — max wall-strain fraction (default 0.12)
    n_frames         : int    — recorded snapshots
    Oh               : float  — Ohnesorge number (None → auto)
    strain_rate_ratio: float  — v_wall / c_wave (quasi-static loading speed)
    verbose          : bool

    Returns
    -------
    frames : list of dicts  (see module docstring)
    spec   : ParticleSpec
    """
    from src.epd.particles import ParticleSpec
    spec = ParticleSpec(count=2, type='elastic', nu=nu, N_nodes=N, Oh=Oh)
    d    = spec.derived
    frames = _run_squeeze_core(d, N, R0, delta_max, n_frames, strain_rate_ratio, verbose)
    return frames, spec


def run_squeeze_raw(q, tau, N, R0=1.0, delta_max=0.14, n_frames=80,
                    C_factor=3000.0, alpha_damp=2.0, strain_rate_ratio=0.01,
                    verbose=False):
    """
    Low-level entry point: raw mechanical parameters (q, τ) as input.
    Used by calibration_sweep.py to generate the q→ν lookup table.

    Parameters
    ----------
    q            : float — area stiffness ratio K_area / El_t
    tau          : float — bending parameter (τ = √(12b), b=0.2 → τ≈1.549)
    N            : int   — number of perimeter nodes
    R0           : float — particle radius
    delta_max    : float — max wall-strain fraction
    n_frames     : int   — recorded snapshots
    C_factor     : float — contact hardness prefactor (C = C_factor*(1+q))
    alpha_damp   : float — damping coefficient
    strain_rate_ratio : float — v_wall / c_wave
    verbose      : bool

    Returns
    -------
    frames : list of dicts  (same keys as run_squeeze output)
    """
    El_t   = 12.0 / tau**2
    K_area = q * El_t
    C      = C_factor * (1.0 + q)
    d = dict(El_t=El_t, K_area=K_area, C=C, tau=tau, alpha=alpha_damp)
    return _run_squeeze_core(d, N, R0, delta_max, n_frames, strain_rate_ratio, verbose)


def interpolate_at_strain(frames, eps_ref):
    """
    Interpolate nu_meas and dA_frac at a reference wall_strain = eps_ref.
    Returns (nu_meas, dA_frac) or (nan, nan) if eps_ref is out of range.
    """
    ws  = np.array([f['wall_strain'] for f in frames])
    nu  = np.array([f['nu_meas']     for f in frames])
    dA  = np.array([f['dA_frac']     for f in frames])

    ok = np.isfinite(nu) & (ws > 1e-6)
    if ok.sum() < 2 or eps_ref < ws[ok].min() or eps_ref > ws[ok].max():
        return float('nan'), float('nan')

    nu_at = float(np.interp(eps_ref, ws[ok], nu[ok]))
    dA_at = float(np.interp(eps_ref, ws[ok], dA[ok]))
    return nu_at, dA_at


if __name__ == '__main__':
    import time
    print("twodisk_squeeze v1 — quick test")
    t0 = time.time()
    frames, spec = run_squeeze(nu=0.5, N=32, delta_max=0.10, n_frames=20, verbose=True)
    print(f"  {len(frames)} frames in {time.time()-t0:.1f}s")

    nu_at, dA_at = interpolate_at_strain(frames, eps_ref=0.08)
    print(f"  At eps_ref=0.08:  nu_meas={nu_at:.4f}  dA={dA_at:.4f}")
    assert np.isfinite(nu_at),      "nu_meas should be finite"
    assert 0.0 < nu_at < 1.0,      f"nu_meas out of range: {nu_at}"
    assert abs(dA_at) < 0.05,      f"area not conserved: dA={dA_at:.4f}"
    print("  PASS")

    # Also test raw path (q=2.0 → roughly ν≈0.7)
    print("\nraw path test (q=2.0, N=32) ...")
    t1 = time.time()
    TAU = np.sqrt(12.0 * 0.2)
    frames_raw = run_squeeze_raw(q=2.0, tau=TAU, N=32, delta_max=0.10,
                                 n_frames=20, verbose=False)
    nu_r, dA_r = interpolate_at_strain(frames_raw, eps_ref=0.08)
    print(f"  {len(frames_raw)} frames in {time.time()-t1:.1f}s"
          f"  nu_meas@0.08={nu_r:.4f}  dA={dA_r:.4f}")
    assert np.isfinite(nu_r), "raw path: nu_meas should be finite"
    print("  PASS")
