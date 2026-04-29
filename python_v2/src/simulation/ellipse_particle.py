"""
ellipse_particle.py — EllipseParticle: CapsuleParticle with elliptical reference shape.

Equal-arc-length discretisation of the ellipse perimeter ensures every edge has
the same reference length L0 = perimeter/N, so the membrane tension is uniform at
rest.  The first node is pinned to the positive semi-major axis: (a, 0).

Usage
-----
    from src.simulation.ellipse_particle import EllipseParticle, ellipse_arclength_nodes

    p = EllipseParticle(a=1.5, b=0.8, N=32, tau=0.05, S=1.0, C=3000.0, rho_d=1.0)
    p.set_center([-3.0, 0.0])
    p.set_velocity([0.1, 0.0], omega=0.5)

    sim = CapsuleSim([p], primitives=[])

Everything downstream (step_rb, contact forces, movies) works unchanged because
EllipseParticle inherits CapsuleParticle's full interface.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.simulation.capsule_shell import CapsuleParticle


# ── Equal-arc-length parameterisation ────────────────────────────────────────

def ellipse_arclength_nodes(a, b, N, n_quad=8192):
    """
    Return N positions on ellipse x=a cos(t), y=b sin(t) with equal arc-length
    spacing.  First node is at t=0: (a, 0) — the positive semi-major axis tip.

    Algorithm
    ---------
    1. Tabulate ds/dt = sqrt((a sin t)^2 + (b cos t)^2) on a fine grid.
    2. Integrate cumulatively (trapezoidal) to get s(t).
    3. Invert s(t) to find t at each equally-spaced s value (linear interpolation).

    Parameters
    ----------
    a, b   : float — semi-major and semi-minor axes  (a >= b > 0)
    N      : int   — number of nodes
    n_quad : int   — quadrature points (8192 gives <1e-6 relative error in node pos)

    Returns
    -------
    xy       : (N, 2) node positions in body frame (centred at origin)
    perimeter: float — true arc length ≈ trapezoidal integral
    L0       : float — arc-length per segment = perimeter / N
    """
    from scipy.interpolate import interp1d

    t  = np.linspace(0.0, 2.0 * np.pi, n_quad, endpoint=False)
    dt = 2.0 * np.pi / n_quad

    # ds/dt: arc-length element
    dsdt = np.sqrt((a * np.sin(t)) ** 2 + (b * np.cos(t)) ** 2)

    # Cumulative arc length via trapezoidal rule
    s_cum = np.zeros(n_quad + 1)
    s_cum[1:] = np.cumsum(dsdt) * dt
    perimeter = float(s_cum[-1])
    L0 = perimeter / N

    # t_extended[k] corresponds to s_cum[k]  (s goes from 0 to perimeter)
    t_extended = np.append(t, 2.0 * np.pi)
    t_of_s = interp1d(s_cum, t_extended, kind='linear', assume_sorted=True)

    # Node arc-length targets: 0, L0, 2*L0, ..., (N-1)*L0
    s_nodes = np.arange(N, dtype=float) * L0
    thetas  = t_of_s(s_nodes)

    xy = np.column_stack([a * np.cos(thetas), b * np.sin(thetas)])
    return xy, perimeter, L0


# ── EllipseParticle ──────────────────────────────────────────────────────────

class EllipseParticle(CapsuleParticle):
    """
    CapsuleParticle whose reference shape is an ellipse (a >= b).

    The parent class is initialised with an *effective radius* R_eff = sqrt(a*b)
    (so that the elastic parameters scale correctly), then the reference geometry
    is patched to the equal-arc-length ellipse.

    Overridden attributes
    ---------------------
    X_ref   : (N,2) — body-frame reference positions on ellipse
    L0      : arc-length per segment = perimeter(a,b) / N
    r_c     : contact bead radius = L0  (same convention as CapsuleParticle)
    A0      : shoelace area of reference polygon
    theta0  : reference turning angles (vary around the ellipse)
    M_disk  : rho_d * pi * a * b
    I_disk  : M_disk * (a^2 + b^2) / 4  (solid uniform ellipse disk)
    R0      : max(a, b)  — used only for the broad-phase bounding check

    Parameters
    ----------
    a, b   : float — semi-major and semi-minor axes
    N      : int   — node count
    tau    : float — shell thickness ratio (controls elastic stiffness)
    S      : float — force scale
    C      : float — contact spring constant (penalty stiffness)
    rho_d  : float — area density (mass per unit area)
    center : (2,)  — initial centre position
    """

    def __init__(self, a, b, N, tau=0.05, S=1.0, C=3000.0, rho_d=1.0,
                 center=(0.0, 0.0)):
        if a < b:
            raise ValueError("Convention: a >= b (a = semi-major axis)")

        # Initialise parent at origin with effective radius R_eff = sqrt(a*b)
        R_eff = float(np.sqrt(a * b))
        super().__init__(N=N, R0=R_eff, tau=tau, S=S, C=C, rho_d=rho_d,
                         center=(0.0, 0.0))

        self.a = float(a)
        self.b = float(b)

        # ── Patch reference geometry to ellipse (body frame, origin) ──────────
        xy_body, perimeter, L0 = ellipse_arclength_nodes(a, b, N)

        self.X_ref  = xy_body.copy()
        self.x      = xy_body.copy()         # world positions at origin
        self.L0     = float(L0)
        self.r_c    = float(L0)
        self.A0     = self.area()            # shoelace area of ellipse polygon
        self.theta0 = self._compute_turning_angles()   # reference curvature angles

        # ── Patch mass / inertia ──────────────────────────────────────────────
        self.M_disk = float(rho_d * np.pi * a * b)
        self.I_disk = float(self.M_disk * (a ** 2 + b ** 2) / 4.0)

        # Broadphase in _capsule_capsule_forces uses R0 as bounding radius:
        # use semi-major axis to be conservative (never miss a collision)
        self.R0 = float(max(a, b))

        # ── Reset rigid-body state at origin ──────────────────────────────────
        self.x_cm  = np.zeros(2)
        self.v_cm  = np.zeros(2)
        self.theta = 0.0
        self.omega = 0.0
        self.u     = np.zeros((N, 2))
        self.u_dot = np.zeros((N, 2))

        # Move to requested centre
        cx, cy = float(center[0]), float(center[1])
        if cx != 0.0 or cy != 0.0:
            self.translate(np.array([cx, cy]))
