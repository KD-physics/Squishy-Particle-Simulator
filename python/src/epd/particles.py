"""
particles.py — ParticleSpec: high-level particle population specification.

Usage
-----
    from src.epd.particles import ParticleSpec

    # Elastic particles by effective Poisson ratio
    spec = ParticleSpec(count=20, nu=0.5)

    # Elastic particles by low-level q
    spec = ParticleSpec(count=20, q=2.0, tau_b=0.2)

    # Polydisperse
    spec = ParticleSpec(count=50, nu=0.6, poly_dist=0.05)

    # Build CapsuleParticle list (called by System)
    particles = spec.build(seed=42, centers=[(x0,y0), ...])
"""

import os
import json
import numpy as np

# Path to calibration table (N=32, eps_ref=0.08)
_CALIB_PATH = os.path.join(
    os.path.dirname(__file__), '..', '..', 'results',
    'calibration_sweep', 'calibration_data.json')

_calib_cache = None   # (nu_arr, q_arr) loaded once


def _load_calibration():
    global _calib_cache
    if _calib_cache is not None:
        return _calib_cache
    with open(os.path.abspath(_CALIB_PATH)) as f:
        data = json.load(f)
    entries = [(e['q'], e['metrics']['0.08']['nu'])
               for e in data if e['N'] == 32 and '0.08' in e.get('metrics', {})]
    entries.sort(key=lambda x: x[0])    # sort by q
    q_arr  = np.array([e[0] for e in entries])
    nu_arr = np.array([e[1] for e in entries])
    _calib_cache = (nu_arr, q_arr)
    return _calib_cache


def nu_to_q(nu_target):
    """Interpolate q from nu using the N=32, eps_ref=0.08 calibration table."""
    nu_arr, q_arr = _load_calibration()
    return float(np.exp(np.interp(float(nu_target), nu_arr, np.log(q_arr))))


def q_to_nu(q_val):
    """Inverse: interpolate nu from q."""
    nu_arr, q_arr = _load_calibration()
    return float(np.interp(np.log(float(q_val)), np.log(q_arr), nu_arr))


class ParticleSpec:
    """
    Specification for a population of EPD particles.

    Parameters (all keyword except count)
    ----------
    count        : int   — number of particles in this population
    type         : str   — 'elastic' | 'emulsion' | 'rigid'

    High-level physical params (preferred):
      nu         : float — effective Poisson ratio (elastic); derives q
      Bo         : float — Bond number (emulsion only); derives gamma/K_area

    Low-level overrides (used directly if specified):
      q          : float — K_area / El_t; overrides nu if given
      tau_b      : float — bending working point (default 0.2)
      alpha_damp : float — damping coefficient; auto if None
      C          : float — contact hardness; auto (3000*S*(1+q)) if None

    Size distribution:
      N_nodes    : int   — perimeter nodes (default 60)
      R0_mean    : float — mean radius (always normalised to 1.0 exactly)
      poly_dist  : float | dict | None
                   None           → monodisperse (all R0=1.0)
                   0.05           → Gaussian sigma=0.05 (shorthand)
                   {'type': 'gaussian', 'sigma': 0.05}
                   {'type': 'bimodal',  'ratio': 0.5, 'delta': 0.1}
                   {'type': 'explicit', 'values': [...]}  # len must == count

    Motion / shape:
      motion     : MotionSpec — drives CM velocity; None = free particle
      frozen_shape: bool — freeze elastic DOFs (shape_frozen in params)

    Extensible force parameters (all zero/no-op by default):
      extra_forces: dict — e.g. {'drag': 0.1, 'activity': 0.05}
    """

    def __init__(self,
                 count,
                 type='elastic',
                 nu=None,
                 Bo=None,
                 q=None,
                 kappa=None,
                 tau_b=0.2,
                 alpha_damp=None,
                 C=None,
                 gamma=None,
                 Oh=None,
                 N_nodes=60,
                 R0_mean=1.0,
                 poly_dist=None,
                 motion=None,
                 frozen_shape=False,
                 extra_forces=None):
        self.count        = int(count)
        self.type         = str(type)
        self.nu           = float(nu)    if nu    is not None else None
        self.Bo           = float(Bo)    if Bo    is not None else None
        self.q            = float(q)     if q     is not None else None
        self.kappa        = float(kappa) if kappa is not None else None
        self.tau_b        = float(tau_b)
        self.alpha_damp   = float(alpha_damp) if alpha_damp is not None else None
        self.C            = float(C)   if C          is not None else None
        self.gamma        = float(gamma) if gamma    is not None else None
        self.Oh           = float(Oh)  if Oh         is not None else None
        self.N_nodes      = int(N_nodes)
        self.R0_mean      = float(R0_mean)
        self.poly_dist    = poly_dist
        self.motion       = motion
        self.frozen_shape = bool(frozen_shape)
        self.extra_forces = dict(extra_forces) if extra_forces is not None else {}
        self._system      = None   # set by System.add_particles() after init

        # rigid → always frozen shape
        if self.type == 'rigid':
            self.frozen_shape = True

        # Derive material parameters
        self._q_eff, self._TAU, self._El_t, self._K_area, self._C_eff, self._alpha_eff \
            = self._derive_material()
        self._xi_drag = self._derive_xi()

    # ── material derivation ───────────────────────────────────────────────────

    def _derive_material(self, S=1.0):
        """Derive (q, TAU, El_t, K_area, C, alpha) from spec parameters."""
        t = self.type

        if t == 'elastic':
            # q: nu → q (log-interpolate), or q given directly, or default=2.0
            if self.q is not None:
                q = self.q
            elif self.nu is not None:
                q = nu_to_q(self.nu)
            else:
                q = 2.0   # default mid-range
            TAU    = np.sqrt(12.0 * self.tau_b)
            El_t   = 12.0 * S / TAU**2
            K_area = q * El_t
            C      = self.C if self.C is not None else 3000.0 * S * (1.0 + q)
            alpha  = self.alpha_damp if self.alpha_damp is not None else 2.0

        elif t == 'emulsion':
            # gamma (line tension) is the physical energy/force scale
            gamma  = self.gamma if self.gamma is not None else 1.0 * S
            # kappa = gamma / (R0 * K_area) — area compressibility ratio (0 < kappa < 1)
            # User supplies kappa directly; q=1/kappa is only an internal convenience.
            # Default kappa=0.02 (near-incompressible working point from paper).
            if self.kappa is not None:
                kappa = float(np.clip(self.kappa, 1e-3, 0.99))
            elif self.q is not None:
                kappa = 1.0 / max(self.q, 1.01)
            else:
                kappa = 0.02   # paper working point
            q      = 1.0 / kappa
            K_area = gamma / (kappa * self.R0_mean)   # = q * gamma / R0
            TAU    = 0.0    # no bending for emulsion
            El_t   = 0.0    # no edge spring; line tension handled separately
            # C̃ = C·R0/γ = 500 is the paper working point (Table 2, §emulsion).
            # Independent of κ — do NOT scale with q.
            C      = self.C if self.C is not None else 500.0 * gamma / self.R0_mean
            alpha  = self.alpha_damp if self.alpha_damp is not None else 5.0

        elif t == 'rigid':
            q      = 1.0    # irrelevant; shape is frozen
            TAU    = np.sqrt(12.0 * self.tau_b)
            El_t   = 12.0 * S / TAU**2
            K_area = El_t   # moderate area stiffness
            C      = self.C if self.C is not None else 1.0e6 * S
            alpha  = self.alpha_damp if self.alpha_damp is not None else 2.0

        else:
            raise ValueError(f"Unknown particle type: '{t}'. Use 'elastic', 'emulsion', or 'rigid'.")

        return float(q), float(TAU), float(El_t), float(K_area), float(C), float(alpha)

    def _derive_xi(self):
        """Derive drag per unit arc length from Oh and particle reference velocity."""
        if self.Oh is None:
            return 0.0
        Oh = float(self.Oh)
        if self.type == 'emulsion':
            # v_ref = sqrt(gamma / (rho_d * R0)) = 1 at working point
            gamma = float(self.gamma if self.gamma is not None else 1.0)
            v_ref = float(np.sqrt(gamma / self.R0_mean))
            return Oh * v_ref
        elif self.type in ('elastic', 'rigid'):
            # v_ref = sqrt(El_t / (rho_d * R0))
            El_t = self._El_t
            if El_t > 0:
                return Oh * float(np.sqrt(El_t / self.R0_mean))
            return 0.0
        return 0.0

    def _get_gamma_eff(self):
        """Effective line tension (emulsion only; 0 for elastic/rigid)."""
        if self.type == 'emulsion':
            return float(self.gamma if self.gamma is not None else 1.0)
        return 0.0

    @property
    def derived(self):
        """Dict of derived material parameters."""
        d = {
            'TAU':     self._TAU,
            'tau':     self._TAU / np.sqrt(12),
            'El_t':    self._El_t,
            'K_area':  self._K_area,
            'C':       self._C_eff,
            'alpha':   self._alpha_eff,
            'gamma':   self._get_gamma_eff(),
            'Oh':      self.Oh,
            'xi':      self._xi_drag,
        }
        if self.type == 'elastic':
            d['q'] = self._q_eff              # q = K_area / El_t
        elif self.type == 'emulsion':
            d['kappa'] = 1.0 / self._q_eff   # κ = γ/(R0·K_area) — primary emulsion parameter
            d['q']     = self._q_eff          # = 1/κ, kept for backward compatibility
        return d

    # ── drag / Oh interface ───────────────────────────────────────────────────

    def terminal_velocity(self, g):
        """
        Compute terminal velocity for a single particle under gravity g.
        Requires Oh to be set (xi_drag > 0) and g > 0.
        v_t = g * M_disk / (xi * L_perimeter)  ≈  g * rho_d * pi * R0 / (2 * xi)
        """
        if self.Oh is None or self._xi_drag == 0.0:
            raise ValueError("Oh not set — no drag defined for this spec.")
        if g == 0.0:
            raise ValueError("g=0 — terminal velocity undefined without gravity.")
        # For a circle of radius R0: perimeter = 2*pi*R0, mass = rho_d*pi*R0^2
        # v_t = F_gravity / F_drag_coeff = (rho_d*pi*R0^2*g) / (xi*2*pi*R0)
        #      = rho_d * R0 * g / (2 * xi)   [with rho_d=R0=1 at working point]
        return float(g * self.R0_mean / (2.0 * self._xi_drag))

    def set_terminal_velocity(self, v_t, g):
        """
        Set Oh such that terminal velocity equals v_t under gravity g.
        Back-computes xi → Oh and updates the cascade.
        """
        if g == 0.0:
            raise ValueError("g=0 — cannot define terminal velocity without gravity.")
        if v_t <= 0.0:
            raise ValueError("v_t must be positive.")
        # v_t = rho_d * R0 * g / (2 * xi)  →  xi = rho_d * R0 * g / (2 * v_t)
        xi = float(self.R0_mean * g / (2.0 * v_t))
        self.set_xi(xi)

    def set_Oh(self, Oh):
        """Set Ohnesorge number and propagate to xi and downstream params."""
        self.Oh = float(Oh)
        self._xi_drag = self._derive_xi()
        self._push_to_system('xi_drag_per_p', self._xi_drag)

    def set_xi(self, xi):
        """Set drag per unit arc length directly, propagate Oh upward."""
        self._xi_drag = float(xi)
        # Back-compute Oh from xi
        if self.type == 'emulsion':
            gamma = float(self.gamma if self.gamma is not None else 1.0)
            v_ref = float(np.sqrt(gamma / self.R0_mean))
        elif self.type in ('elastic', 'rigid'):
            El_t = self._El_t
            v_ref = float(np.sqrt(El_t / self.R0_mean)) if El_t > 0 else 1.0
        else:
            v_ref = 1.0
        self.Oh = self._xi_drag / v_ref if v_ref > 0 else None
        self._push_to_system('xi_drag_per_p', self._xi_drag)

    def set_kappa(self, kappa):
        """Set kappa (emulsion) and propagate K_area, C, xi, Oh."""
        if self.type != 'emulsion':
            raise ValueError("set_kappa only valid for emulsion particles.")
        self.kappa = float(kappa)
        self.q = None   # kappa takes precedence
        self._q_eff, self._TAU, self._El_t, self._K_area, self._C_eff, self._alpha_eff \
            = self._derive_material()
        self._xi_drag = self._derive_xi()
        self._push_to_system('K_area', self._K_area)
        self._push_to_system('xi_drag_per_p', self._xi_drag)

    def set_nu(self, nu):
        """Set nu (elastic) and propagate q, El_t, K_area, C, xi."""
        if self.type not in ('elastic', 'rigid'):
            raise ValueError("set_nu only valid for elastic/rigid particles.")
        self.nu = float(nu)
        self.q  = None   # nu takes precedence
        self._q_eff, self._TAU, self._El_t, self._K_area, self._C_eff, self._alpha_eff \
            = self._derive_material()
        self._xi_drag = self._derive_xi()
        self._push_to_system('K_area',  self._K_area)
        self._push_to_system('El_t',    self._El_t)
        self._push_to_system('EI',      self._El_t * (self.R0_mean / (2*np.pi/self.N_nodes))**2 / 12.0)
        self._push_to_system('xi_drag_per_p', self._xi_drag)

    def _push_to_system(self, param_name, value):
        """If this spec is attached to a live System, patch the TF params tensor."""
        if self._system is not None:
            self._system._push_param_for_spec(self, param_name, value)

    def set_motion(self, motion_spec):
        """Attach a MotionSpec to drive all particles in this spec."""
        self.motion = motion_spec
        return self

    # ── size distribution ─────────────────────────────────────────────────────

    def sample_R0(self, seed=42):
        """
        Sample per-particle R0 values (length=self.count), normalised so mean=1.0.

        Returns
        -------
        R0_arr : (count,) float64 array, mean=1.0 exactly
        """
        rng = np.random.default_rng(seed)
        pd  = self.poly_dist
        n   = self.count

        if pd is None:
            return np.ones(n)

        # Shorthand: float → Gaussian sigma
        if isinstance(pd, (int, float)):
            pd = {'type': 'gaussian', 'sigma': float(pd)}

        dist_type = pd.get('type', 'gaussian')

        if dist_type == 'gaussian':
            sigma = float(pd.get('sigma', 0.05))
            raw   = rng.normal(1.0, sigma, size=n)
            raw   = np.clip(raw, 1e-2, None)     # no negative radii

        elif dist_type == 'bimodal':
            # ratio = fraction of large particles, delta = half-gap
            ratio = float(pd.get('ratio', 0.5))
            delta = float(pd.get('delta', 0.1))
            # Sizes: 1-delta (small) and 1+delta (large)
            n_large = int(round(ratio * n))
            n_small = n - n_large
            raw = np.concatenate([
                np.full(n_small, 1.0 - delta),
                np.full(n_large, 1.0 + delta),
            ])
            rng.shuffle(raw)

        elif dist_type == 'explicit':
            vals = np.asarray(pd['values'], dtype=float)
            if len(vals) != n:
                raise ValueError(
                    f"poly_dist explicit values length {len(vals)} != count {n}")
            raw = vals.copy()

        else:
            raise ValueError(f"Unknown poly_dist type: '{dist_type}'")

        # Normalise so mean=1.0 exactly
        raw /= raw.mean()
        return raw.astype(np.float64)

    # ── build CapsuleParticle list ────────────────────────────────────────────

    def build(self, seed=42, centers=None):
        """
        Instantiate CapsuleParticle (or EmulsionParticle) objects.

        Parameters
        ----------
        seed    : int — RNG seed for size distribution sampling
        centers : list of (x, y) tuples, length == count, or None
                  (if None, particles placed at origin; caller must set positions)

        Returns
        -------
        list of CapsuleParticle (or subtype)
        """
        from src.simulation.capsule_shell import CapsuleParticle

        R0_arr  = self.sample_R0(seed=seed)
        d       = self.derived

        if centers is None:
            centers = [(0.0, 0.0)] * self.count

        particles = []
        for i in range(self.count):
            R0_i = float(R0_arr[i])
            cx, cy = float(centers[i][0]), float(centers[i][1])

            if self.type == 'emulsion':
                from src.simulation.emulsion_particle import EmulsionParticle
                # γ ∝ R0^{+1}; K_area and C are R0-independent (preserves κ and C̃).
                gamma_eff = d['gamma'] * R0_i / self.R0_mean
                p = EmulsionParticle(
                    N=self.N_nodes,
                    R0=R0_i,
                    gamma=gamma_eff,
                    K_area=d['K_area'],          # constant — κ = γ/(R0·K_area) preserved
                    C=d['C'],                    # constant — C̃ = C·R0/γ preserved
                    rho_d=1.0,
                    center=(cx, cy),
                )
                # Stamp dimensionless targets for future recompute calls
                p._kappa_target = 1.0 / self._q_eff
                p._C_tilde      = d['C'] * self.R0_mean / d['gamma']  # C̃ = 500
                p._Oh_target    = self.Oh
            else:
                p = CapsuleParticle(
                    N=self.N_nodes,
                    R0=R0_i,
                    tau=d['tau'],
                    S=1.0,
                    C=d['C'],
                    K_area=d['K_area'],
                    center=(cx, cy),
                )
                # Stamp dimensionless targets for future recompute calls
                p._nu_target  = self.nu
                p._q_target   = self._q_eff
                p._Oh_target  = self.Oh

            # Per-particle alpha_damp (∝ R0^{-1}) and xi_drag (∝ R0^{+1})
            p.alpha_damp = self._alpha_eff * self.R0_mean / R0_i
            p.xi_drag    = self._xi_drag   * R0_i / self.R0_mean

            particles.append(p)

        return particles

    # ── display ───────────────────────────────────────────────────────────────

    def summary(self):
        """Return a human-readable summary string."""
        d = self.derived
        lines = [
            f"ParticleSpec(count={self.count}, type='{self.type}')",
            f"  nu={self.nu}  →  q={d['q']:.4f}",
            f"  TAU={d['TAU']:.4f}  El_t={d['El_t']:.4f}  K_area={d['K_area']:.4f}",
            f"  C={d['C']:.2f}  alpha={d['alpha']:.4f}",
            f"  N_nodes={self.N_nodes}  poly_dist={self.poly_dist}",
            f"  frozen_shape={self.frozen_shape}  extra_forces={self.extra_forces}",
        ]
        return '\n'.join(lines)
