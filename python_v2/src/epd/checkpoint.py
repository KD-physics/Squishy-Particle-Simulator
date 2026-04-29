"""
checkpoint.py — Save/load system state for bit-identical resume.

Checkpoint format (directory):
    state.npz    — all numpy arrays (x_all, x_cm, v_cm, omega, theta, u, u_dot,
                   X_ref, t, step, CapCandidates)
    config.json  — {Lx, Ly, periodic_x, periodic_y, dt, alpha_damp, skin, E_candidates}
    specs.json   — list of ParticleSpec serialized dicts
    objects.json — list of SimulationObject serialized dicts
    motion_samples.npz — pre-sampled motion arrays (if any motion uses from_samples)
"""

import os, json
import numpy as np


# ── ParticleSpec serialization ────────────────────────────────────────────────

def serialize_spec(spec):
    """Convert ParticleSpec to a JSON-safe dict."""
    return {
        'count':        spec.count,
        'type':         spec.type,
        'nu':           spec.nu,
        'Bo':           spec.Bo,
        'q':            spec.q,
        'tau_b':        spec.tau_b,
        'alpha_damp':   spec.alpha_damp,
        'C':            spec.C,
        'gamma':        getattr(spec, 'gamma', None),
        'N_nodes':      spec.N_nodes,
        'R0_mean':      spec.R0_mean,
        'poly_dist':    _serialize_poly_dist(spec.poly_dist),
        'frozen_shape': spec.frozen_shape,
        'extra_forces': spec.extra_forces,
        'motion':       serialize_motion(spec.motion) if spec.motion else None,
    }


def deserialize_spec(d):
    """Reconstruct ParticleSpec from dict."""
    from src.epd.particles import ParticleSpec
    motion = deserialize_motion(d.get('motion')) if d.get('motion') else None
    return ParticleSpec(
        count=d['count'],
        type=d.get('type', 'elastic'),
        nu=d.get('nu'),
        Bo=d.get('Bo'),
        q=d.get('q'),
        tau_b=d.get('tau_b', 0.2),
        alpha_damp=d.get('alpha_damp'),
        C=d.get('C'),
        gamma=d.get('gamma'),
        N_nodes=d.get('N_nodes', 60),
        R0_mean=d.get('R0_mean', 1.0),
        poly_dist=d.get('poly_dist'),
        frozen_shape=d.get('frozen_shape', False),
        extra_forces=d.get('extra_forces') or {},
        motion=motion,
    )


def _serialize_poly_dist(pd):
    if pd is None:
        return None
    if isinstance(pd, (int, float)):
        return float(pd)
    return dict(pd)   # dict is already JSON-safe (values are primitives)


# ── MotionSpec serialization ──────────────────────────────────────────────────

def serialize_motion(ms):
    """
    Serialize MotionSpec to JSON-safe dict.
    TF-native callables cannot be serialized — they are stored as None.
    Pre-sampled arrays are stored with a special key (loaded separately from npz).
    """
    if ms is None:
        return None

    from src.epd.motion import MotionSpec

    if ms._sampled_t is not None:
        return {
            'mode': 'sampled',
            'dt': ms._sampled_dt,
            'duration': float(ms._sampled_t[-1]),
            'r_ref': list(ms.r_ref),
        }

    if ms._vx_fn is not None or ms._vy_fn is not None or ms._omega_fn is not None:
        return {
            'mode': 'callable',
            'note': 'TF-native callable — not serializable; re-register after load',
        }

    # Parametric DC+AC
    return {
        'mode':               'parametric',
        'vx_dc':              ms.vx_dc,
        'vy_dc':              ms.vy_dc,
        'vx_ac':              ms.vx_ac,
        'vy_ac':              ms.vy_ac,
        'freq_x':             ms.freq_x,
        'freq_y':             ms.freq_y,
        'phase_x':            ms.phase_x,
        'phase_y':            ms.phase_y,
        'omega_dc':           ms.omega_dc,
        'omega_ac':           ms.omega_ac,
        'freq_omega':         ms.freq_omega,
        'phase_omega':        ms.phase_omega,
        'omega_orbit_dc':     getattr(ms, 'omega_orbit_dc',    0.0),
        'omega_orbit_ac':     getattr(ms, 'omega_orbit_ac',    0.0),
        'omega_orbit_freq':   getattr(ms, 'omega_orbit_freq',  0.0),
        'omega_orbit_phase':  getattr(ms, 'omega_orbit_phase', 0.0),
        'r_ref':              list(ms.r_ref),
    }


def deserialize_motion(d, motion_samples=None):
    """Reconstruct MotionSpec from dict. motion_samples: npz for 'sampled' mode."""
    if d is None:
        return None

    from src.epd.motion import MotionSpec

    mode = d.get('mode', 'parametric')

    if mode == 'parametric':
        return MotionSpec(
            vx_dc=d.get('vx_dc', 0.0), vy_dc=d.get('vy_dc', 0.0),
            vx_ac=d.get('vx_ac', 0.0), vy_ac=d.get('vy_ac', 0.0),
            freq_x=d.get('freq_x', 0.0), freq_y=d.get('freq_y', 0.0),
            phase_x=d.get('phase_x', 0.0), phase_y=d.get('phase_y', 0.0),
            omega_dc=d.get('omega_dc', 0.0),
            omega_ac=d.get('omega_ac', 0.0),
            freq_omega=d.get('freq_omega', 0.0),
            phase_omega=d.get('phase_omega', 0.0),
            omega_orbit_dc=d.get('omega_orbit_dc', 0.0),
            omega_orbit_ac=d.get('omega_orbit_ac', 0.0),
            omega_orbit_freq=d.get('omega_orbit_freq', 0.0),
            omega_orbit_phase=d.get('omega_orbit_phase', 0.0),
            r_ref=tuple(d.get('r_ref', [0.0, 0.0])),
        )

    elif mode == 'sampled' and motion_samples is not None:
        tag = d.get('tag', '0')
        ms  = MotionSpec(r_ref=tuple(d.get('r_ref', [0.0, 0.0])))
        ms._sampled_dt    = float(d['dt'])
        ms._sampled_vx    = motion_samples[f'{tag}_vx'].astype(np.float64)
        ms._sampled_vy    = motion_samples[f'{tag}_vy'].astype(np.float64)
        ms._sampled_omega = motion_samples[f'{tag}_omega'].astype(np.float64)
        ms._sampled_t     = np.arange(len(ms._sampled_vx)) * ms._sampled_dt
        return ms

    elif mode == 'callable':
        return None   # cannot restore; user must re-register

    return None


# ── Object serialization ──────────────────────────────────────────────────────

def serialize_object(obj):
    """Convert SimulationObject to JSON-safe dict (metadata + geometry)."""
    d = {
        'kind':      obj._kind,
        'exclusion': obj._exclusion,
        'k_pen':     obj._k_pen,
        'r_c_wall':  obj._r_c_wall,
        'motion':    serialize_motion(obj._motion),
    }
    # Geometry fields per kind
    if obj._kind == 'Wall':
        d['p0']     = obj._p0.tolist()
        d['p1']     = obj._p1.tolist()
        d['normal'] = obj._normal.tolist()
    elif obj._kind == 'ArcWall':
        d['center'] = obj._center.tolist()
        d['radius'] = obj._radius
        d['convex'] = obj._convex
    elif obj._kind in ('Box', 'Channel', 'CouetteCell', 'CircleObstacle',
                       'RegularPolygon', 'CompositeObject', 'CustomObject'):
        # Composite: save origin, theta, exclusion; children serialized recursively
        d['origin']   = list(obj._origin) if hasattr(obj, '_origin') else [0.0, 0.0]
        d['theta']    = float(obj._theta)  if hasattr(obj, '_theta')  else 0.0
        d['children'] = [serialize_object(c) for c in obj._children] \
                        if hasattr(obj, '_children') else []
        # Kind-specific constructor params
        if obj._kind == 'Box':
            d['width']  = obj._width
            d['height'] = obj._height
            d['x0']     = float(obj._origin[0]) if hasattr(obj, '_origin') else 0.0
            d['y0']     = float(obj._origin[1]) if hasattr(obj, '_origin') else 0.0
        elif obj._kind == 'CircleObstacle':
            d['radius'] = obj._radius
            d['x0']     = float(obj._x0) if hasattr(obj, '_x0') else float(obj._origin[0])
            d['y0']     = float(obj._y0) if hasattr(obj, '_y0') else float(obj._origin[1])
    return d


def deserialize_object(d):
    """Reconstruct a SimulationObject from a dict saved by serialize_object."""
    from src.epd.objects import (Wall, ArcWall, Box, CircleObstacle,
                                  CompositeObject)
    import numpy as np
    kind   = d.get('kind', 'Wall')
    motion = deserialize_motion(d.get('motion')) if d.get('motion') else None

    if kind == 'Wall':
        obj = Wall(p0=d['p0'], p1=d['p1'], normal=d['normal'])
    elif kind == 'ArcWall':
        obj = ArcWall(center=d['center'], radius=d['radius'], convex=d.get('convex', True))
    elif kind == 'Box':
        obj = Box(width=d['width'], height=d['height'],
                  x0=d.get('x0', d['width']/2), y0=d.get('y0', d['height']/2),
                  exclusion=d.get('exclusion', 'exterior'))
    elif kind == 'CircleObstacle':
        obj = CircleObstacle(radius=d['radius'],
                             x0=d.get('x0', 0.0), y0=d.get('y0', 0.0),
                             exclusion=d.get('exclusion', 'interior'))
    else:
        # Fallback: return a bare SimulationObject (cannot reconstruct geometry)
        from src.epd.objects import SimulationObject
        obj = SimulationObject(kind=kind)

    obj._k_pen    = d.get('k_pen',    1.0)
    obj._r_c_wall = d.get('r_c_wall', 0.0)
    if d.get('exclusion') is not None:
        obj._exclusion = d['exclusion']
    if motion is not None:
        obj._motion = motion
    return obj


# ── Checkpoint I/O ────────────────────────────────────────────────────────────

def save_checkpoint(path, system):
    """
    Save full system state to directory `path`.

    Parameters
    ----------
    path   : str or Path — output directory (created if absent)
    system : System — the system to save
    """
    import pathlib
    path = pathlib.Path(path)
    path.mkdir(parents=True, exist_ok=True)

    st  = system._state
    cfg = system._config

    # 1. State arrays
    state_arrays = {}
    for k, v in st.items():
        arr = v.numpy() if hasattr(v, 'numpy') else np.asarray(v)
        state_arrays[k] = arr
    # Include candidacy table
    state_arrays['CapCandidates'] = system._cm_mgr.CapCandidates
    # Save per-particle physics params for bit-identical restore (avoids re-seeding)
    p_keys = ['L0', 'A0', 'K_area', 'El_t', 'EI', 'm_node', 'M_disk', 'I_disk',
              'rho_d', 'theta0', 'r_c_per_p', 'k_c_per_p', 'gamma_lt',
              'driven_mask', 'shape_frozen', 'traj', 'box']
    for k in p_keys:
        if k in system._params:
            arr = system._params[k]
            state_arrays[f'params_{k}'] = arr.numpy() if hasattr(arr, 'numpy') else np.asarray(arr)
    np.savez(path / 'state.npz', **state_arrays)

    # 2. Config
    config_data = {
        'Lx':           float(system.Lx),
        'Ly':           float(system.Ly),
        'periodic_x':   cfg['periodic_x'],
        'periodic_y':   cfg['periodic_y'],
        'dt':           float(system._dt),
        'alpha_damp':   float(system._alpha_damp),
        'skin':         float(system._skin),
        'E_candidates': int(system._E_candidates),
        'g':            float(getattr(system, '_g_val', 0.0)),
    }
    with open(path / 'config.json', 'w') as f:
        json.dump(config_data, f, indent=2)

    # 3. Specs
    specs_data = [serialize_spec(sp) for sp in system._specs]
    with open(path / 'specs.json', 'w') as f:
        json.dump(specs_data, f, indent=2)

    # 4. Objects (minimal — motion only; geometry reconstructed from kind)
    objects_data = [serialize_object(obj) for obj in system._objects]
    with open(path / 'objects.json', 'w') as f:
        json.dump(objects_data, f, indent=2)

    # 5. Motion samples (pre-sampled mode)
    motion_sample_arrays = {}
    for i, sp in enumerate(system._specs):
        if sp.motion is not None and sp.motion._sampled_t is not None:
            tag = str(i)
            motion_sample_arrays[f'{tag}_vx']    = sp.motion._sampled_vx
            motion_sample_arrays[f'{tag}_vy']    = sp.motion._sampled_vy
            motion_sample_arrays[f'{tag}_omega'] = sp.motion._sampled_omega
    if motion_sample_arrays:
        np.savez(path / 'motion_samples.npz', **motion_sample_arrays)

    return path


def load_checkpoint(path, system=None):
    """
    Load checkpoint from directory `path` into `system` (or a new System).

    Returns the populated System instance.
    """
    import pathlib
    from src.epd.system import System

    path = pathlib.Path(path)

    # 1. Config
    with open(path / 'config.json') as f:
        cfg = json.load(f)

    # 2. State
    state_npz   = np.load(path / 'state.npz')
    state_keys  = ['x_all', 'x_cm', 'v_cm', 'omega', 'theta', 'u', 'u_dot', 'X_ref']
    import tensorflow as tf
    DTYPE = tf.float64
    state = {k: tf.constant(state_npz[k], dtype=DTYPE) for k in state_keys}
    # Scalar fields
    for k in ['t', 'step']:
        if k in state_npz:
            state[k] = tf.constant(state_npz[k], dtype=DTYPE if k=='t' else tf.int64)

    # 3. Motion samples (optional)
    ms_path = path / 'motion_samples.npz'
    motion_samples = np.load(ms_path) if ms_path.exists() else None

    # 4. Specs
    with open(path / 'specs.json') as f:
        specs_data = json.load(f)
    specs = [deserialize_spec(d) for d in specs_data]

    # 4b. Objects
    objects_path = path / 'objects.json'
    if objects_path.exists():
        with open(objects_path) as f:
            objects_data = json.load(f)
        objects = [deserialize_object(d) for d in objects_data]
    else:
        objects = []

    # 5. Build or reuse System
    if system is None:
        sys_new = System(
            Lx=cfg['Lx'], Ly=cfg['Ly'],
            periodic_x=cfg['periodic_x'],
            periodic_y=cfg['periodic_y'],
            dt_factor=None,
            alpha_damp=cfg['alpha_damp'],
            E_candidates=cfg['E_candidates'],
            skin=cfg['skin'],
        )
    else:
        sys_new = system

    # Set system internals directly (bypass initialize)
    sys_new.Lx = cfg['Lx']
    sys_new.Ly = cfg['Ly']
    sys_new._config = {
        'periodic_x': cfg['periodic_x'],
        'periodic_y': cfg['periodic_y'],
    }
    sys_new._dt         = cfg['dt']
    sys_new._alpha_damp = cfg['alpha_damp']
    sys_new._skin       = cfg['skin']
    sys_new._E_candidates = cfg['E_candidates']
    sys_new._specs   = specs
    sys_new._objects = objects

    # Rebuild particles from specs at positions from state
    x_cm_np   = state_npz['x_cm']       # (P, 2)
    centers   = [(float(x_cm_np[i, 0]), float(x_cm_np[i, 1]))
                 for i in range(len(x_cm_np))]

    particles = []
    idx = 0
    for sp in specs:
        ps = sp.build(centers=centers[idx:idx+sp.count])
        particles.extend(ps)
        idx += sp.count

    sys_new._particles = particles
    sys_new._state     = state

    # Rebuild params from particles
    from src.simulation.tf_sim import make_state, set_periodic_box
    _, params = make_state(particles)
    set_periodic_box(params, cfg['Lx'], cfg['Ly'])

    # Restore per-particle physics params from saved checkpoint (bit-identical resume).
    # The freshly-built particles may have different R0/seed than the original, so we
    # overwrite the critical physics arrays with the saved values.
    _p_keys_restore = ['L0', 'A0', 'K_area', 'El_t', 'EI', 'm_node', 'M_disk', 'I_disk',
                       'rho_d', 'theta0', 'r_c_per_p', 'k_c_per_p', 'gamma_lt']
    for k in _p_keys_restore:
        saved_k = f'params_{k}'
        if saved_k in state_npz:
            params[k] = tf.constant(state_npz[saved_k], dtype=DTYPE)
    sys_new._params = params

    # Store dt/alpha in params for adaptive_swell compatibility
    import src.simulation.tf_sim as tf_sim_mod
    NP_DTYPE = tf_sim_mod.NP_DTYPE
    params['_dt_tf']    = tf.constant(NP_DTYPE(cfg['dt']))
    params['_alpha_tf'] = tf.constant(NP_DTYPE(cfg['alpha_damp']))

    # Set shape_frozen if any spec has frozen_shape
    _frozen = np.zeros(len(particles))
    idx2 = 0
    for sp in specs:
        if sp.frozen_shape:
            _frozen[idx2:idx2+sp.count] = 1.0
        idx2 += sp.count
    params['shape_frozen'] = tf.constant(_frozen, dtype=DTYPE)

    # Rebuild CandidacyManager
    P = len(particles)
    N = particles[0].N
    R0_arr = np.array([p.R0 for p in particles])
    from src.simulation.candidacy_manager import CandidacyManager
    cm_mgr = CandidacyManager(
        P=P, N=N, R0=float(np.mean(R0_arr)), E=cfg['E_candidates'],
        skin=cfg['skin'],
        periodic=cfg['periodic_x'] and cfg['periodic_y'],
        periodic_x=cfg['periodic_x'], periodic_y=cfg['periodic_y'],
        Lx=cfg['Lx'], Ly=cfg['Ly'],
        R0_arr=R0_arr,
    )
    cm_mgr.update(x_cm_np, state_npz['theta'])
    sys_new._cm_mgr = cm_mgr

    # Build prim_data from objects (re-use if objects were re-registered)
    from src.simulation.tf_sim import make_prim_data
    prim_list = []
    for obj in sys_new._objects:
        t_val = float(state_npz['t']) if 't' in state_npz else 0.0
        prim_list.extend(obj.to_make_prim_list(t=t_val))
    sys_new._prim_data = make_prim_data(prim_list)

    # TF scalars
    g_val = float(cfg.get('g', 0.0))
    sys_new._g_val    = g_val
    sys_new._dt_tf    = tf.constant(NP_DTYPE(cfg['dt']))
    sys_new._alpha_tf = tf.constant(NP_DTYPE(cfg['alpha_damp']))
    sys_new._g_tf     = tf.constant(NP_DTYPE(g_val))

    # Wire driven particles from specs (restore set_driven state)
    from src.simulation.tf_sim import set_driven
    _driven_idxs  = []
    _traj_rows    = []
    _frozen_flags = []
    part_offset   = 0
    for sp in specs:
        if sp.motion is not None:
            for k in range(sp.count):
                _driven_idxs.append(part_offset + k)
                _traj_rows.append(sp.motion.to_traj_row())
                _frozen_flags.append(sp.frozen_shape)
        part_offset += sp.count
    if _driven_idxs:
        set_driven(params, _driven_idxs, _traj_rows, frozen=_frozen_flags)

    return sys_new
