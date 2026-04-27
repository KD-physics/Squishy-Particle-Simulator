"""
motion.py — MotionSpec for EPD simulation objects and particles.

Two specification modes:

  Mode 1 — parametric (DC + AC):
      MotionSpec(vx=1.0)
      MotionSpec(vx_dc=0.5, vx_ac=1.0, freq_x=2.0)
      MotionSpec(omega_dc=0.3, r_ref=(5.0, 5.0))

  Mode 2 — TF-native callable (preferred for GPU):
      MotionSpec(vx=lambda t: tf.sin(omega*t),
                 vy=0.0, omega=0.0)

  Mode 3 — pre-sampled (Python fallback):
      MotionSpec.from_samples(vx_fn=lambda t: np.cos(t),
                              dt=0.001, duration=10.0)

Interface
---------
  .velocity(t)   → (vx, vy, omega)  Python floats (parametric / callable modes)
  .displacement(t) → (dx, dy, dtheta)  analytic integral (parametric mode only)
  .resolve_tf(t) → (vx_tf, vy_tf, omega_tf, r_ref_tf)  TF tensors (all modes)
  .is_static()   → bool
"""

import numpy as np


# Sentinel to distinguish "not provided" from 0.0
_MISSING = object()


class MotionSpec:
    """
    Rigid-body motion profile.

    Parameters (all keyword)
    -----------------------
    vx, vy      : float | callable(t→scalar) — DC velocity or time-dependent fn
    omega       : float | callable(t→scalar) — spin velocity (rad/s)
    r_ref       : (2,) — reference point for orbital motion (default [0,0])

    DC+AC parametric (used when vx/vy/omega are scalars):
    vx_dc, vy_dc  : DC velocity (alias for scalar vx/vy)
    vx_ac, vy_ac  : AC amplitude
    freq_x, freq_y, freq_omega : AC frequency (rad/s)
    phase_x, phase_y, phase_omega : AC initial phase (rad)
    omega_dc      : DC spin (alias for scalar omega)
    omega_ac      : AC spin amplitude

    When vx/vy/omega are callables, all AC parameters are ignored.
    Callable form: f(t) where t can be a Python float or tf.Tensor.
    """

    def __init__(self,
                 vx=0.0, vy=0.0, omega=0.0,
                 r_ref=(0.0, 0.0),
                 # DC+AC expansion (only used when vx/vy/omega are scalars)
                 vx_dc=None, vy_dc=None, omega_dc=None,
                 vx_ac=0.0, vy_ac=0.0,
                 freq_x=0.0, freq_y=0.0,
                 phase_x=0.0, phase_y=0.0,
                 omega_ac=0.0,
                 freq_omega=0.0,
                 phase_omega=0.0,
                 # Orbital motion: orbit about r_ref at angular velocity omega_orbit_dc
                 omega_orbit_dc=0.0, omega_orbit_ac=0.0,
                 omega_orbit_freq=0.0, omega_orbit_phase=0.0):

        self.r_ref = np.asarray(r_ref, dtype=float)

        # Determine if any component is callable
        self._vx_fn    = vx    if callable(vx)    else None
        self._vy_fn    = vy    if callable(vy)    else None
        self._omega_fn = omega if callable(omega) else None

        # Always initialize scalar params (used when fn is None; ignored otherwise)
        self.vx_dc     = float(vx_dc  if vx_dc  is not None else (0.0 if callable(vx)    else vx))
        self.vx_ac     = float(vx_ac)
        self.freq_x    = float(freq_x)
        self.phase_x   = float(phase_x)
        self.vy_dc     = float(vy_dc  if vy_dc  is not None else (0.0 if callable(vy)    else vy))
        self.vy_ac     = float(vy_ac)
        self.freq_y    = float(freq_y)
        self.phase_y   = float(phase_y)
        self.omega_dc  = float(omega_dc if omega_dc is not None else (0.0 if callable(omega) else omega))
        self.omega_ac  = float(omega_ac)
        self.freq_omega  = float(freq_omega)
        self.phase_omega = float(phase_omega)

        # Orbital motion parameters
        self.omega_orbit_dc    = float(omega_orbit_dc)
        self.omega_orbit_ac    = float(omega_orbit_ac)
        self.omega_orbit_freq  = float(omega_orbit_freq)
        self.omega_orbit_phase = float(omega_orbit_phase)

        # Pre-sampled arrays (filled by from_samples classmethod)
        self._sampled_t     = None   # (M,) float64
        self._sampled_vx    = None   # (M,) float64
        self._sampled_vy    = None   # (M,) float64
        self._sampled_omega = None   # (M,) float64
        self._sampled_dt    = None   # float — for index computation

        # TF tensors (built lazily by resolve_tf)
        self._tf_vx    = None
        self._tf_vy    = None
        self._tf_omega = None

    # ── classmethod constructors ──────────────────────────────────────────────

    @classmethod
    def from_samples(cls, vx_fn=None, vy_fn=None, omega_fn=None,
                     dt=0.001, duration=10.0, r_ref=(0.0, 0.0)):
        """
        Pre-sample velocity functions at dt intervals.

        resolve_tf(t) interpolates linearly in the sampled arrays.

        Parameters
        ----------
        vx_fn, vy_fn, omega_fn : callable f(t) → float  (or None → 0)
        dt       : float — sample interval
        duration : float — total pre-sample duration; behaviour beyond duration
                   clamps to the last sample value.
        """
        ms = cls(r_ref=r_ref)   # bare init
        M  = int(np.ceil(duration / dt)) + 1
        t_arr = np.arange(M, dtype=np.float64) * dt

        ms._sampled_vx    = np.array([vx_fn(t)    if vx_fn    else 0.0 for t in t_arr])
        ms._sampled_vy    = np.array([vy_fn(t)    if vy_fn    else 0.0 for t in t_arr])
        ms._sampled_omega = np.array([omega_fn(t) if omega_fn else 0.0 for t in t_arr])
        ms._sampled_t     = t_arr
        ms._sampled_dt    = float(dt)
        return ms

    # ── velocity / displacement ───────────────────────────────────────────────

    def _eval_component(self, fn, dc, ac, freq, phase, t):
        """Evaluate one velocity component at time t (Python float)."""
        t = float(t)
        if fn is not None:
            return float(fn(t))
        if self._sampled_vx is not None:
            # Use sampled array (first call this method for vx channel)
            pass   # handled separately below
        val = dc
        if freq != 0.0:
            val += ac * np.cos(freq * t + phase)
        elif ac != 0.0:
            val += ac * np.cos(phase)
        return float(val)

    def velocity(self, t):
        """Instantaneous velocity (vx, vy, omega) at time t (Python floats)."""
        t = float(t)

        if self._sampled_t is not None:
            idx   = t / self._sampled_dt
            i0    = int(np.clip(idx, 0, len(self._sampled_t) - 2))
            frac  = idx - i0
            frac  = float(np.clip(frac, 0.0, 1.0))
            vx  = float(self._sampled_vx[i0]    * (1 - frac) + self._sampled_vx[i0+1]    * frac)
            vy  = float(self._sampled_vy[i0]    * (1 - frac) + self._sampled_vy[i0+1]    * frac)
            om  = float(self._sampled_omega[i0] * (1 - frac) + self._sampled_omega[i0+1] * frac)
            return vx, vy, om

        vx  = self._eval_component(self._vx_fn,    self.vx_dc,    self.vx_ac,    self.freq_x,    self.phase_x,    t)
        vy  = self._eval_component(self._vy_fn,    self.vy_dc,    self.vy_ac,    self.freq_y,    self.phase_y,    t)
        om  = self._eval_component(self._omega_fn, self.omega_dc, self.omega_ac, self.freq_omega, self.phase_omega, t)
        return vx, vy, om

    def displacement(self, t):
        """
        Cumulative displacement (dx, dy, dtheta) from t=0 to t.
        Analytic integral (parametric mode only).
        Callable and sampled modes use numerical integration (trapezoid).
        """
        t = float(t)

        if self._sampled_t is not None or self._vx_fn is not None or \
           self._vy_fn is not None or self._omega_fn is not None:
            # Numerical integration via trapezoid rule at dt=0.01 resolution
            dt_integ = min(0.01, t / 100 + 1e-12)
            n_steps  = max(1, int(np.ceil(t / dt_integ)))
            ts       = np.linspace(0.0, t, n_steps + 1)
            vxs, vys, oms = zip(*[self.velocity(ti) for ti in ts])
            dx     = float(np.trapz(vxs, ts))
            dy     = float(np.trapz(vys, ts))
            dtheta = float(np.trapz(oms, ts))
            return dx, dy, dtheta

        # Analytic integral for parametric mode
        dx     = self.vx_dc * t
        dy     = self.vy_dc * t
        dtheta = self.omega_dc * t

        if self.freq_x != 0.0:
            dx += (self.vx_ac / self.freq_x) * (
                np.sin(self.freq_x * t + self.phase_x) - np.sin(self.phase_x))
        elif self.vx_ac != 0.0:
            dx += self.vx_ac * np.cos(self.phase_x) * t

        if self.freq_y != 0.0:
            dy += (self.vy_ac / self.freq_y) * (
                np.sin(self.freq_y * t + self.phase_y) - np.sin(self.phase_y))
        elif self.vy_ac != 0.0:
            dy += self.vy_ac * np.cos(self.phase_y) * t

        if self.freq_omega != 0.0:
            dtheta += (self.omega_ac / self.freq_omega) * (
                np.sin(self.freq_omega * t + self.phase_omega) - np.sin(self.phase_omega))
        elif self.omega_ac != 0.0:
            dtheta += self.omega_ac * np.cos(self.phase_omega) * t

        return float(dx), float(dy), float(dtheta)

    # ── TF interface ──────────────────────────────────────────────────────────

    def resolve_tf(self, t):
        """
        Evaluate motion at absolute time t, returning TF tensors.

        Returns (vx_tf, vy_tf, omega_tf, r_ref_tf) as tf.float64 tensors.
        t can be a Python float or a tf.Tensor (float64).
        """
        import tensorflow as tf

        DTYPE = tf.float64

        r_ref_tf = tf.constant(self.r_ref, dtype=DTYPE)

        if self._sampled_t is not None:
            # Mode 3: pre-sampled linear interpolation in TF
            vx_tf, vy_tf, om_tf = self._resolve_sampled_tf(t, DTYPE)
            return vx_tf, vy_tf, om_tf, r_ref_tf

        # Mode 1/2: callable or parametric
        def _eval_tf(fn, dc, ac, freq, phase):
            if fn is not None:
                val = fn(t)
                return tf.cast(val, DTYPE) if not isinstance(val, tf.Tensor) else val
            t_val = tf.cast(t, DTYPE) if not isinstance(t, tf.Tensor) else t
            result = tf.constant(dc, dtype=DTYPE)
            if freq != 0.0:
                result = result + tf.constant(ac, dtype=DTYPE) * tf.cos(
                    tf.constant(freq, dtype=DTYPE) * t_val + tf.constant(phase, dtype=DTYPE))
            elif ac != 0.0:
                result = result + tf.constant(ac * np.cos(phase), dtype=DTYPE)
            return result

        vx_tf = _eval_tf(self._vx_fn,    self.vx_dc,    self.vx_ac,    self.freq_x,    self.phase_x)
        vy_tf = _eval_tf(self._vy_fn,    self.vy_dc,    self.vy_ac,    self.freq_y,    self.phase_y)
        om_tf = _eval_tf(self._omega_fn, self.omega_dc, self.omega_ac, self.freq_omega, self.phase_omega)

        return vx_tf, vy_tf, om_tf, r_ref_tf

    def _resolve_sampled_tf(self, t, DTYPE):
        """Linear interpolation in pre-sampled arrays using TF ops."""
        import tensorflow as tf

        t_tf  = tf.cast(t, DTYPE) if not isinstance(t, tf.Tensor) else t
        vx_s  = tf.constant(self._sampled_vx,    dtype=DTYPE)
        vy_s  = tf.constant(self._sampled_vy,    dtype=DTYPE)
        om_s  = tf.constant(self._sampled_omega, dtype=DTYPE)
        dt_tf = tf.constant(self._sampled_dt, dtype=DTYPE)
        N     = tf.shape(vx_s)[0]

        idx_f  = t_tf / dt_tf
        i0     = tf.cast(tf.clip_by_value(tf.floor(idx_f), 0, tf.cast(N - 2, DTYPE)), tf.int32)
        frac   = tf.clip_by_value(idx_f - tf.cast(i0, DTYPE), 0.0, 1.0)

        vx_tf = vx_s[i0] * (1.0 - frac) + vx_s[i0 + 1] * frac
        vy_tf = vy_s[i0] * (1.0 - frac) + vy_s[i0 + 1] * frac
        om_tf = om_s[i0] * (1.0 - frac) + om_s[i0 + 1] * frac
        return vx_tf, vy_tf, om_tf

    # ── utilities ─────────────────────────────────────────────────────────────

    def is_static(self):
        """True if motion is identically zero for all t."""
        if self._vx_fn or self._vy_fn or self._omega_fn:
            return False
        if self._sampled_t is not None:
            return (np.all(self._sampled_vx == 0.0) and
                    np.all(self._sampled_vy == 0.0) and
                    np.all(self._sampled_omega == 0.0))
        return (self.vx_dc == 0.0 and self.vy_dc == 0.0 and self.omega_dc == 0.0 and
                self.vx_ac == 0.0 and self.vy_ac == 0.0 and self.omega_ac == 0.0 and
                self.omega_orbit_dc == 0.0 and self.omega_orbit_ac == 0.0)

    def to_traj_row(self, r_ref=None):
        """
        Convert parametric MotionSpec to (18,) traj array for make_traj / set_driven.
        Only works for parametric mode (not callable/sampled).
        """
        from src.simulation.tf_sim import make_traj
        r = self.r_ref if r_ref is None else np.asarray(r_ref, dtype=float)
        return make_traj(
            v_dc=(self.vx_dc, self.vy_dc),
            v_ac=(self.vx_ac, self.vy_ac),
            freq=(self.freq_x, self.freq_y),
            phase=(self.phase_x, self.phase_y),
            omega_spin_dc=self.omega_dc,
            omega_spin_ac=self.omega_ac,
            omega_spin_freq=self.freq_omega,
            omega_spin_phase=self.phase_omega,
            omega_orbit_dc=self.omega_orbit_dc,
            omega_orbit_ac=self.omega_orbit_ac,
            omega_orbit_freq=self.omega_orbit_freq,
            omega_orbit_phase=self.omega_orbit_phase,
            r_ref=tuple(r),
        )

    def is_parametric(self):
        """True if MotionSpec is purely parametric (DC+AC scalars, no callables/samples)."""
        return (self._vx_fn is None and self._vy_fn is None and
                self._omega_fn is None and self._sampled_t is None)

    def rescale(self, f):
        """
        Scale spatial quantities by factor f (called during box compression).
        Positions (r_ref) and translational velocities (vx, vy) scale with box size.
        Angular velocity (omega) is scale-invariant and unchanged.
        """
        self.r_ref   = self.r_ref * float(f)
        self.vx_dc  *= float(f)
        self.vy_dc  *= float(f)
        self.vx_ac  *= float(f)
        self.vy_ac  *= float(f)
