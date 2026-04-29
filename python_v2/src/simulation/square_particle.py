"""
square_particle.py — SquareParticle: CapsuleParticle with square reference shape.

N nodes equally spaced around the perimeter of a square with half-side a.
Exactly 4 nodes land on the corners: indices 0, N//4, N//2, 3*N//4.
Node 0 is at the bottom-right corner (a, -a); nodes proceed counter-clockwise.

Mass/inertia for a uniform solid square:
    M = 4 * rho_d * a^2
    I = M * (2a)^2 / 6 = 2 * M * a^2 / 3    (solid square, side 2a)
    R0 = a * sqrt(2)   (circumradius, for broadphase)
    L0 = 8a / N        (arc-length per segment, perimeter = 8a)
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.simulation.capsule_shell import CapsuleParticle


def square_nodes(a, N):
    """
    Return N equally-spaced nodes around square perimeter (half-side a),
    counter-clockwise from bottom-right corner (a, -a).

    N must be divisible by 4. Corners are at indices 0, N//4, N//2, 3*N//4.

    Returns
    -------
    xy : (N, 2) — body-frame reference positions
    L0 : float  — arc-length per segment = 8a / N
    """
    if N % 4 != 0:
        raise ValueError(f"N must be divisible by 4, got {N}")

    n = N // 4          # nodes per side
    L0 = 8.0 * a / N   # = 2a / n

    nodes = []
    # Side 1: right edge, bottom to top  (a, -a) → (a, a)
    for k in range(n):
        nodes.append([a,  -a + k * L0])
    # Side 2: top edge, right to left  (a, a) → (-a, a)
    for k in range(n):
        nodes.append([a - k * L0,  a])
    # Side 3: left edge, top to bottom  (-a, a) → (-a, -a)
    for k in range(n):
        nodes.append([-a,  a - k * L0])
    # Side 4: bottom edge, left to right  (-a, -a) → (a, -a)
    for k in range(n):
        nodes.append([-a + k * L0,  -a])

    return np.array(nodes, dtype=float), float(L0)


class SquareParticle(CapsuleParticle):
    """
    CapsuleParticle whose reference shape is a square with half-side a.

    The parent is initialised with R_eff = a (so elastic parameters are sensible),
    then reference geometry and mass/inertia are patched.

    Parameters
    ----------
    a      : float — half-side length
    N      : int   — node count (must be divisible by 4)
    tau    : float — shell thickness ratio
    S      : float — force scale
    C      : float — contact spring constant
    rho_d  : float — area density (mass per unit area)
    center : (2,)  — initial centre position
    theta  : float — initial orientation (radians)
    """

    def __init__(self, a, N=32, tau=0.05, S=1.0, C=3000.0, rho_d=1.0,
                 center=(0.0, 0.0), theta=0.0):
        if N % 4 != 0:
            raise ValueError(f"N must be divisible by 4, got {N}")

        self.a = float(a)

        # Initialise parent at origin; R_eff = a
        super().__init__(N=N, R0=float(a), tau=tau, S=S, C=C, rho_d=rho_d,
                         center=(0.0, 0.0))

        # ── Patch reference geometry ──────────────────────────────────────────
        xy_body, L0 = square_nodes(a, N)

        self.X_ref  = xy_body.copy()
        self.x      = xy_body.copy()
        self.L0     = L0
        self.r_c    = L0
        self.A0     = self.area()
        self.theta0 = self._compute_turning_angles()

        # ── Patch mass / inertia ──────────────────────────────────────────────
        side = 2.0 * a
        self.M_disk = float(rho_d * side ** 2)          # = 4 rho a^2
        self.I_disk = float(self.M_disk * side ** 2 / 6.0)  # = 2 M a^2 / 3

        # Circumradius for broadphase bounding check
        self.R0 = float(a * np.sqrt(2.0))

        # ── Reset rigid-body state ────────────────────────────────────────────
        self.x_cm  = np.zeros(2)
        self.v_cm  = np.zeros(2)
        self.theta = 0.0
        self.omega = 0.0
        self.u     = np.zeros((N, 2))
        self.u_dot = np.zeros((N, 2))

        # Apply initial orientation then translate to requested centre
        if theta != 0.0:
            self.rotate_reference(theta)

        cx, cy = float(center[0]), float(center[1])
        if cx != 0.0 or cy != 0.0:
            self.translate(np.array([cx, cy]))

    def rotate_reference(self, angle):
        """Set initial orientation: rotates world positions self.x; X_ref stays canonical."""
        c, s = np.cos(angle), np.sin(angle)
        R = np.array([[c, -s], [s, c]])
        # X_ref stays in body frame; x is reconstructed from x_cm + R(theta) @ X_ref
        self.x     = self.X_ref @ R.T   # world positions (x_cm=0 here)
        self.theta = float(angle)
