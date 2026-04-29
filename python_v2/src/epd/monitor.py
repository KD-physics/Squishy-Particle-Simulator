"""
monitor.py — SimMonitor: real-time health monitoring for EPD simulations.

Checks every time it is called:
  • max_overlap   — max fractional capsule overlap (0=touching, 1=coincident)
  • max_node_spd  — max nodal speed (simulation units / time)
  • min_circ      — min per-particle circularity  4πA/L²  (1=circle, 0=degenerate)
  • max_stretch   — max fractional edge stretch  |L/L0 - 1|
  • max_cm_speed  — max CM speed
  • max_kink      — max exterior turning angle at any perimeter node (radians)
                    0=straight, π/3≈1.05 rad=triangular corner (60° exterior angle)
  • ring_ratio    — max|u_dot| / v_mem  (elastic deformation speed / membrane wave speed)
                    u_dot is the elastic DOF velocity relative to the rigid body frame

Action thresholds (warn / critical):
  overlap  warn=0.45  crit=0.65
  circ     warn=0.55  crit=0.30
  stretch  warn=0.40  crit=0.80
  node_spd scaled by v_ref = mean(r_c) / dt   warn=30×  crit=100×
  kink     warn=π/6≈0.52 rad (30°)  crit=π/3≈1.05 rad (60°, triangular)
  ring     warn=1.0×v_mem  crit=3.0×v_mem  (membrane wave speed reference)
"""

import numpy as np


class SimMonitor:
    """
    Compute and store health metrics for one simulation step.

    Parameters
    ----------
    f_warn, f_crit         : overlap thresholds
    circ_warn, circ_crit   : circularity thresholds (4πA/L²)
    stretch_warn/crit      : edge-stretch thresholds
    spd_warn, spd_crit     : node-speed multiples of v_ref = r_c/dt
    kink_warn, kink_crit   : max exterior angle thresholds (radians)
    ring_warn, ring_crit   : ring_ratio thresholds (multiples of membrane wave speed)
    """

    def __init__(self,
                 f_warn=0.45,      f_crit=0.65,
                 circ_warn=0.55,   circ_crit=0.30,
                 stretch_warn=0.40, stretch_crit=0.80,
                 spd_warn=30.0,    spd_crit=100.0,
                 kink_warn=np.pi/6, kink_crit=np.pi/3,
                 ring_warn=1.0,    ring_crit=3.0):
        self.f_warn        = f_warn
        self.f_crit        = f_crit
        self.circ_warn     = circ_warn
        self.circ_crit     = circ_crit
        self.stretch_warn  = stretch_warn
        self.stretch_crit  = stretch_crit
        self.spd_warn      = spd_warn
        self.spd_crit      = spd_crit
        self.kink_warn     = kink_warn
        self.kink_crit     = kink_crit
        self.ring_warn     = ring_warn
        self.ring_crit     = ring_crit

        # Last computed metrics (dict)
        self.metrics = {}

    # ── geometry helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _circularity(x_all):
        """Per-particle circularity = 4πA/L². Returns (P,) array."""
        P, N, _ = x_all.shape
        circs = np.zeros(P)
        for i in range(P):
            xy = x_all[i]
            xs, ys = xy[:, 0], xy[:, 1]
            A = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(np.roll(xs, -1), ys))
            edges = np.roll(xy, -1, axis=0) - xy
            L = np.sum(np.linalg.norm(edges, axis=1))
            circs[i] = 4.0 * np.pi * A / (L ** 2) if L > 1e-12 else 0.0
        return circs

    @staticmethod
    def _edge_stretch(x_all, L0_arr):
        """Max fractional edge stretch per particle; returns (P,) array."""
        P, N, _ = x_all.shape
        max_s = np.zeros(P)
        for i in range(P):
            xy = x_all[i]
            edges = np.roll(xy, -1, axis=0) - xy
            L = np.linalg.norm(edges, axis=1)
            max_s[i] = float(np.max(np.abs(L / L0_arr[i] - 1.0)))
        return max_s

    @staticmethod
    def _max_kink(x_all):
        """
        Max exterior turning angle at any perimeter node across all particles (radians).

        For a smooth convex polygon with N nodes the natural angle is 2π/N ≈ 0.196 rad
        (N=32).  A "kink" is a sharp corner: angles near π/3 (60°) look triangular.

        Returns scalar (max over all particles and nodes).
        """
        P, N, _ = x_all.shape
        max_kink = 0.0
        for i in range(P):
            xy = x_all[i]
            # Incoming and outgoing edge unit vectors at each node
            e_in  = xy - np.roll(xy,  1, axis=0)   # node[i] - node[i-1]
            e_out = np.roll(xy, -1, axis=0) - xy   # node[i+1] - node[i]
            len_in  = np.linalg.norm(e_in,  axis=1, keepdims=True)
            len_out = np.linalg.norm(e_out, axis=1, keepdims=True)
            # Avoid divide-by-zero for degenerate edges
            n_in  = e_in  / np.maximum(len_in,  1e-12)
            n_out = e_out / np.maximum(len_out, 1e-12)
            # cos of exterior angle (angle between consecutive edge directions)
            cos_t = np.clip(np.sum(n_in * n_out, axis=1), -1.0, 1.0)
            angles = np.arccos(cos_t)   # exterior turning angles in [0, π]
            max_kink = max(max_kink, float(angles.max()))
        return max_kink

    @staticmethod
    def _node_speeds(state):
        """
        Per-node speed = CM translation + rotation + elastic deformation.
        Returns (P, N) float array.
        """
        v_cm  = state['v_cm'].numpy() if hasattr(state['v_cm'], 'numpy') else np.asarray(state['v_cm'])
        omega = state['omega'].numpy() if hasattr(state['omega'], 'numpy') else np.asarray(state['omega'])
        u_dot = state['u_dot'].numpy() if hasattr(state['u_dot'], 'numpy') else np.asarray(state['u_dot'])
        x_all = state['x_all'].numpy() if hasattr(state['x_all'], 'numpy') else np.asarray(state['x_all'])
        x_cm  = state['x_cm'].numpy()  if hasattr(state['x_cm'],  'numpy') else np.asarray(state['x_cm'])
        r     = x_all - x_cm[:, None, :]
        v_rot = np.stack([-omega[:, None] * r[:, :, 1],
                           omega[:, None] * r[:, :, 0]], axis=-1)
        v_node = v_cm[:, None, :] + v_rot + u_dot
        return np.linalg.norm(v_node, axis=-1)

    @staticmethod
    def _ring_ratio(state, params):
        """
        ring_ratio = max|u_dot| / v_mem
        where v_mem = R₀ × ω_mem (membrane wave speed reference, stored in params as _v_mem).
        If _v_mem is not stored, returns 0.0 (metric unavailable).
        """
        if '_v_mem' not in params:
            return 0.0
        v_mem = float(params['_v_mem'].numpy() if hasattr(params['_v_mem'], 'numpy')
                      else params['_v_mem'])
        if v_mem < 1e-12:
            return 0.0
        u_dot = state['u_dot'].numpy() if hasattr(state['u_dot'], 'numpy') else np.asarray(state['u_dot'])
        max_udot = float(np.max(np.linalg.norm(u_dot, axis=-1)))
        return max_udot / v_mem

    # ── main check ────────────────────────────────────────────────────────────

    def check(self, state, params, cm_mgr=None, Lx=1e9, Ly=1e9, dt=0.01):
        """
        Compute health metrics from current state and params.

        Returns (metrics_dict, ok_bool, reason_str).
        """
        x_all  = state['x_all'].numpy() if hasattr(state['x_all'], 'numpy') else np.asarray(state['x_all'])
        x_cm   = state['x_cm'].numpy()  if hasattr(state['x_cm'],  'numpy') else np.asarray(state['x_cm'])
        P, N, _ = x_all.shape

        L0_arr  = params['L0'].numpy()         if hasattr(params['L0'],         'numpy') else np.asarray(params['L0'])
        r_c_arr = params['r_c_per_p'].numpy()  if hasattr(params['r_c_per_p'],  'numpy') else np.asarray(params['r_c_per_p'])

        # Circularity
        circs    = self._circularity(x_all)
        min_circ = float(circs.min())

        # Edge stretch
        stretches   = self._edge_stretch(x_all, L0_arr)
        max_stretch = float(stretches.max())

        # Node speeds
        spds         = self._node_speeds(state)
        max_node_spd = float(spds.max())
        v_cm_np      = state['v_cm'].numpy() if hasattr(state['v_cm'], 'numpy') else np.asarray(state['v_cm'])
        max_cm_speed = float(np.linalg.norm(v_cm_np, axis=1).max())
        v_ref        = float(np.mean(r_c_arr)) / max(dt, 1e-12)
        spd_ratio    = max_node_spd / max(v_ref, 1e-12)

        # Perimeter kink (triangularity)
        max_kink = self._max_kink(x_all)

        # Ringing (elastic u_dot vs membrane wave speed)
        ring_ratio = self._ring_ratio(state, params)

        # Overlap (requires candidacy manager)
        max_overlap = 0.0
        if cm_mgr is not None:
            x_flat    = x_all.reshape(P * N, 2)
            caps      = cm_mgr.CapCandidates
            x_pad     = np.vstack([np.zeros((1, 2)), x_flat])
            E         = caps.shape[1]
            cand_flat = caps.ravel()
            xa        = np.repeat(x_flat, E, axis=0)
            xb        = x_pad[cand_flat]
            active    = cand_flat > 0
            if active.any():
                dx = xa[:, 0] - xb[:, 0];  dx -= Lx * np.round(dx / Lx)
                dy = xa[:, 1] - xb[:, 1];  dy -= Ly * np.round(dy / Ly)
                dist       = np.sqrt(dx**2 + dy**2)
                L0_a_rep   = np.repeat(np.repeat(L0_arr, N), E)
                nb_idx     = np.maximum(cand_flat - 1, 0)
                L0_b_rep   = np.repeat(L0_arr, N)[nb_idx]
                L0_contact = 0.5 * (L0_a_rep + L0_b_rep)
                f_vals     = np.where(active,
                                      np.maximum(0.0, (2*L0_contact - dist) / (2*L0_contact)),
                                      0.0)
                max_overlap = float(f_vals.max())

        m = {
            'max_overlap':  max_overlap,
            'max_node_spd': max_node_spd,
            'max_cm_speed': max_cm_speed,
            'min_circ':     min_circ,
            'max_stretch':  max_stretch,
            'spd_ratio':    spd_ratio,
            'v_ref':        v_ref,
            'max_kink':     max_kink,
            'ring_ratio':   ring_ratio,
            'circs':        circs,
            'stretches':    stretches,
        }
        self.metrics = m

        # Evaluate thresholds (first violation wins)
        ok, reason = True, ''
        if max_overlap >= self.f_crit:
            ok = False; reason = f'overlap={max_overlap:.3f} ≥ f_crit={self.f_crit}'
        elif min_circ <= self.circ_crit:
            ok = False; reason = f'min_circ={min_circ:.3f} ≤ circ_crit={self.circ_crit}'
        elif max_stretch >= self.stretch_crit:
            ok = False; reason = f'max_stretch={max_stretch:.3f} ≥ stretch_crit={self.stretch_crit}'
        elif spd_ratio >= self.spd_crit:
            ok = False; reason = f'spd_ratio={spd_ratio:.1f} ≥ spd_crit={self.spd_crit}'
        elif max_kink >= self.kink_crit:
            ok = False; reason = f'max_kink={max_kink:.3f} rad ≥ kink_crit={self.kink_crit:.3f}'
        elif ring_ratio >= self.ring_crit:
            ok = False; reason = f'ring_ratio={ring_ratio:.2f} ≥ ring_crit={self.ring_crit}'

        return m, ok, reason

    def warn_level(self):
        """Return 'ok', 'warn', or 'crit' based on last computed metrics."""
        m = self.metrics
        if not m:
            return 'ok'
        if (m['max_overlap']  >= self.f_crit         or
                m['min_circ']     <= self.circ_crit      or
                m['max_stretch']  >= self.stretch_crit   or
                m['spd_ratio']    >= self.spd_crit        or
                m['max_kink']     >= self.kink_crit       or
                m['ring_ratio']   >= self.ring_crit):
            return 'crit'
        if (m['max_overlap']  >= self.f_warn          or
                m['min_circ']     <= self.circ_warn       or
                m['max_stretch']  >= self.stretch_warn    or
                m['spd_ratio']    >= self.spd_warn         or
                m['max_kink']     >= self.kink_warn        or
                m['ring_ratio']   >= self.ring_warn):
            return 'warn'
        return 'ok'

    def summary_line(self):
        """One-line human-readable summary of last metrics."""
        m = self.metrics
        if not m:
            return '(no metrics yet)'
        return (f"overlap={m['max_overlap']:.3f}  "
                f"circ={m['min_circ']:.3f}  "
                f"stretch={m['max_stretch']:.3f}  "
                f"kink={m['max_kink']:.2f}rad  "
                f"ring×={m['ring_ratio']:.2f}  "
                f"cm_spd={m['max_cm_speed']:.3f}")
