"""
capsule_shell.py — Fluid-filled elastic shell particle model (Capsule Shell DEM)

Each particle is a 2D elastic shell enclosing an incompressible fluid.
The shell is discretized as N perimeter nodes; all forces are purely local
(no K matrix, no global solve). Time integration is explicit Verlet.

Physical model
--------------
Forces per node (all described in CAPSULE_SHELL.md):
  1. Fluid pressure     — global area → instantaneous normal force on all nodes
  2. Hydrostatic gravity — optional ρ_f g h_i term (regularized to zero net x-force)
  3. Edge elasticity    — tangential spring resisting perimeter stretch/compression
  4. Bending            — energy-consistent 3-node hinge resisting curvature change
  5. Capsule–capsule    — 2-point Gauss quadrature over edge pairs between particles
  6. Capsule–primitive  — same quadrature, using gap_and_normal() interface from
                          contact_primitives.py (rigid walls, arcs, polygons)

Dimensionless reference scales: R0 = 1, K_fluid = 1.
Three physical knobs: τ (shell thickness), S (squishiness), C (contact hardness).

Usage
-----
    from capsule_shell import CapsuleParticle, CapsuleSim

    p1 = CapsuleParticle(N=64, R0=1.0, tau=0.05, S=1.0, C=500.0, rho_d=1.0)
    p2 = CapsuleParticle(N=64, R0=1.0, tau=0.05, S=1.0, C=500.0, rho_d=1.0)
    p2.translate([0.0, 2.05])   # place second disk just above first

    sim = CapsuleSim([p1, p2], primitives=[], g=0.0)
    sim.step(dt=1e-3, alpha_damp=5.0)
"""

import numpy as np


# ── geometry helpers ──────────────────────────────────────────────────────────

def _rot90ccw(v):
    """Rotate 2D vector 90° counter-clockwise."""
    return np.array([-v[1], v[0]])


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-15 else v


def _signed_angle(t_prev, t_next):
    """
    Signed turning angle from t_prev to t_next (CCW positive).
    Uses atan2(cross, dot) for robustness near 0 and π.
    """
    cross = t_prev[0] * t_next[1] - t_prev[1] * t_next[0]
    dot   = t_prev[0] * t_next[0] + t_prev[1] * t_next[1]
    return np.arctan2(cross, dot)


def _closest_point_on_segment(p, a, b):
    """Return (closest_point, t) where t ∈ [0,1] is parameter along a→b."""
    ab  = b - a
    ab2 = np.dot(ab, ab)
    if ab2 < 1e-30:
        return a.copy(), 0.0
    t = float(np.clip(np.dot(p - a, ab) / ab2, 0.0, 1.0))
    return a + t * ab, t


def _shoelace_area(x):
    """Signed shoelace area of a polygon with vertices x (N,2)."""
    xr = np.roll(x, -1, axis=0)
    return 0.5 * np.sum(x[:, 0] * xr[:, 1] - xr[:, 0] * x[:, 1])


# ── Gauss quadrature points on [0, 1] ────────────────────────────────────────

_GAUSS2_S = np.array([(1.0 - 1.0 / np.sqrt(3)) / 2.0,
                       (1.0 + 1.0 / np.sqrt(3)) / 2.0])
_GAUSS2_W = np.array([0.5, 0.5])


# ─────────────────────────────────────────────────────────────────────────────
# CapsuleParticle
# ─────────────────────────────────────────────────────────────────────────────

class CapsuleParticle:
    """
    A single fluid-filled elastic shell particle.

    Parameters
    ----------
    N : int
        Number of perimeter nodes.
    R0 : float
        Reference radius (length unit; default 1).
    tau : float
        Shell thickness ratio t/R0.
    S : float
        Squishiness = EI / (K_fluid R0^4). Controls bending vs fluid compliance.
    C : float
        Contact hardness = k_c R0 / K_fluid.
    rho_d : float
        Fluid density (sets node mass and wave timescales).
    K_fluid : float
        Fluid bulk modulus (default 1 — the reference pressure scale).
    center : array-like (2,)
        Initial center position.
    """

    def __init__(self, N=64, R0=1.0, tau=0.05, S=1.0, C=500.0,
                 rho_d=1.0, K_fluid=1.0, K_area=None, center=(0.0, 0.0)):
        self.N       = N
        self._R0     = float(R0)    # backing store; use R0 property for access
        self.tau     = tau
        self.S       = S
        self.C       = C
        self.rho_d   = rho_d
        self.K_fluid = K_fluid
        # K_area: stiffness of the area (fluid pressure) restoring force.
        # Decoupled from K_fluid so the ratio K_area / El_t can be varied
        # independently of the membrane and bending stiffness.
        # Default = K_fluid (backward compatible).
        self.K_area  = K_fluid if K_area is None else float(K_area)

        # Derived shell parameters
        # For a thin shell of thickness t = tau*R0 and plane-strain modulus E_l:
        #   El_t = E_l * t = E_l * tau * R0               (membrane stiffness, ∝ R0)
        #   EI   = E_l * t^3 / 12 = El_t * (tau*R0)^2/12 (bending stiffness, ∝ R0^3)
        #   beta = El_t / (K_fluid * R0) = E_l * tau / K_fluid = 12 S / tau^2
        #   S    = E_l * tau^3 / (12 * K_fluid)           (dimensionless squishiness knob)
        #
        # Correct formula: EI = S * K_fluid * R0**3   (∝ R0^3 for self-similar scaling).
        # Prior code had EI = S * K_fluid * R0**4, which made bending grow as R0 and
        # broke R0 self-similarity.  At R0=1 both give EI=S — calibration is unchanged.
        self.beta    = 12.0 * S / (tau ** 2)            # membrane stiffness ratio
        self.EI      = S * K_fluid * R0 ** 3            # bending stiffness (∝ R0^3 ✓)
        self.El_t    = self.beta * K_fluid * R0          # E_l * t (membrane stiffness ∝ R0 ✓)
        self.k_c     = C * K_fluid / R0                  # contact spring constant

        # Initialize positions on a circle
        theta_ref = 2.0 * np.pi * np.arange(N) / N
        cx, cy = float(center[0]), float(center[1])
        self.x = np.column_stack([
            cx + R0 * np.cos(theta_ref),
            cy + R0 * np.sin(theta_ref)
        ])                                               # (N, 2)
        self.v = np.zeros((N, 2))                        # velocities
        self.f = np.zeros((N, 2))                        # forces (reset each step)

        # Reference geometry — computed from actual initial positions so the
        # undeformed polygon is exactly at equilibrium (chord ≠ arc length).
        _, _, L_init = self.edge_tangents_normals()
        self.L0      = float(L_init.mean())             # mean chord length (all equal for regular polygon)
        self.A0      = self.area()                       # shoelace area of initial polygon

        # Node mass: m_i = rho_d * t * L0  (t = tau * R0)
        self.m_node  = rho_d * (tau * R0) * self.L0

        # Capsule contact radius: r_c = L0 (each edge is a sausage of half-width L0)
        self.r_c     = self.L0

        # Reference turning angles at each node (computed from initial positions)
        self.theta0 = self._compute_turning_angles()

        # Reference center (for diagnostics)
        self.center_ref = np.array([cx, cy])

        # ── Rigid-body decomposition state ────────────────────────────────────
        # These attributes support step_rb() in CapsuleSim.  They are also
        # kept consistent by translate() / set_center() so old code that calls
        # those helpers before switching to step_rb() stays valid.
        #
        # Body-frame reference positions (theta=0, u=0): shape (N, 2)
        self.X_ref   = self.x - np.array([cx, cy])
        # Rigid-body degrees of freedom
        self.x_cm    = np.array([cx, cy], dtype=float)
        self.v_cm    = np.zeros(2)
        self.theta   = 0.0          # orientation angle (rad)
        self.omega   = 0.0          # angular velocity (rad / time)
        # Elastic degrees of freedom (body frame)
        self.u       = np.zeros((N, 2))   # elastic displacement of each node
        self.u_dot   = np.zeros((N, 2))   # elastic velocity of each node
        # Correct solid-disk mass and moment of inertia
        # (node-only gives M_shell = rho_l * 2π R0 and I_shell = M_shell R0^2;
        #  a solid disk of the same "fluid" density has M and I below)
        self.M_disk  = rho_d * np.pi * R0 ** 2
        self.I_disk  = 0.5 * self.M_disk * R0 ** 2

        # Per-particle damping / drag (set by ParticleSpec.build; 0 = not set)
        self.alpha_damp = 0.0   # internal elastic damping coefficient (∝ R0^{-1})
        self.xi_drag    = 0.0   # rigid-body fluid drag per arc-length (∝ R0^{+1})

    # ── R0 property — setter triggers force-constant recompute ───────────────

    @property
    def R0(self):
        return self._R0

    @R0.setter
    def R0(self, value):
        self._R0 = float(value)
        self._recompute_params()

    def _recompute_params(self):
        """
        Recompute all R0-dependent force constants from current _R0 and geometry.
        Called by the R0 property setter; also called by System.adjust_params_for_size().
        Assumes node positions and L0 already reflect the new R0 (geometry is correct).
        """
        R0 = self._R0
        self.EI     = self.S * self.K_fluid * R0 ** 3
        self.El_t   = self.beta * self.K_fluid * R0
        self.k_c    = self.C * self.K_fluid / R0
        self.m_node = self.rho_d * (self.tau * R0) * self.L0
        self.M_disk = self.rho_d * np.pi * R0 ** 2
        self.I_disk = 0.5 * self.M_disk * R0 ** 2

    # ── geometry queries ─────────────────────────────────────────────────────

    def _compute_turning_angles(self):
        """Signed turning angle at each node from current positions (vectorized)."""
        t_hat, _, _ = self.edge_tangents_normals()
        t_prev = np.roll(t_hat, 1, axis=0)
        cross  = t_prev[:, 0] * t_hat[:, 1] - t_prev[:, 1] * t_hat[:, 0]
        dot    = t_prev[:, 0] * t_hat[:, 0] + t_prev[:, 1] * t_hat[:, 1]
        return np.arctan2(cross, dot)

    def edge_tangents_normals(self):
        """
        Returns:
            t_hat (N, 2) — unit tangent of each edge i→(i+1)
            n_hat (N, 2) — outward normal of each edge (rot90ccw of t_hat)
            L     (N,)   — current edge lengths
        """
        N = self.N
        xnext = np.roll(self.x, -1, axis=0)
        e = xnext - self.x                               # edge vectors (N, 2)
        L = np.linalg.norm(e, axis=1)                   # (N,)
        t_hat = e / L[:, None]                           # unit tangents
        n_hat = np.column_stack([-t_hat[:, 1], t_hat[:, 0]])  # rot90ccw
        return t_hat, n_hat, L

    def node_normals(self, n_hat):
        """
        Outward node normals = average of adjacent edge normals (normalized).
        n_hat_node_i = normalize(n_hat_{i-1} + n_hat_i)
        """
        n_prev = np.roll(n_hat, 1, axis=0)
        avg = n_hat + n_prev
        norms = np.linalg.norm(avg, axis=1, keepdims=True)
        norms = np.where(norms < 1e-15, 1.0, norms)
        return avg / norms

    def area(self):
        """Current enclosed area (shoelace)."""
        return abs(_shoelace_area(self.x))

    def center_of_mass(self):
        return self.x.mean(axis=0)

    # ── translation / placement ───────────────────────────────────────────────

    def translate(self, delta):
        """Move all nodes by delta (2,).  Also updates x_cm."""
        delta = np.asarray(delta, dtype=float)
        self.x    += delta
        self.x_cm += delta

    def set_center(self, pos):
        """Move particle so its centroid is at pos.  Also updates x_cm."""
        self.translate(np.asarray(pos) - self.center_of_mass())

    def set_velocity(self, v_cm, omega=0.0):
        """
        Set rigid-body velocity (for rigid-body integrator).
        Does NOT touch nodal velocities self.v (used by old step()).
        """
        self.v_cm  = np.asarray(v_cm, dtype=float).copy()
        self.omega = float(omega)

    # ── force accumulation ────────────────────────────────────────────────────

    def accumulate_internal_forces(self, g=0.0):
        """
        Compute and add to self.f:
          - fluid pressure
          - hydrostatic gravity (if g > 0)
          - edge elasticity
          - bending
        Does NOT touch contact forces.
        """
        t_hat, n_hat, L = self.edge_tangents_normals()
        n_node = self.node_normals(n_hat)
        N = self.N

        # ── 1. Fluid pressure ──────────────────────────────────────────────
        # Physical convention: compressed fluid (A < A0) pushes outward (P > 0).
        # P = K_area * (A0 - A) / A0  →  positive when compressed, outward force.
        # Equivalently: F_i = -dU/dx_i where U = (K_area/2) * (A-A0)^2 / A0
        # and dA/dx_i ≈ L0 * n_hat_node_i_OUTWARD → F_i = -K*(A-A0)/A0 * L0 * n̂_i_out
        #
        # NOTE ON SIGN: n_node computed from rot90CCW(t_hat) gives the INWARD
        # normal for a CCW polygon (see edge_tangents_normals).  The restoring
        # force must be OUTWARD, so we subtract rather than add Fp.
        A  = self.area()
        P  = self.K_area * (self.A0 - A) / self.A0       # positive when compressed
        Fp = P * self.L0 * n_node                         # (N, 2) — n_node is INWARD
        self.f -= Fp                                      # subtract → outward for P>0

        # ── 2. Internal hydrostatic gravity (optional) ────────────────────
        # Internal fluid (density rho_d) exerts outward pressure on the shell.
        # P_k = rho_d * g * (y_top - y_k): pressure from fluid column above node k.
        # Force is outward (+P * n̂_out * ds_k).  Since n_node is INWARD,
        # self.f -= Fg gives the outward contribution.
        # Net result: Σ F_y = -rho_d * g * A = -m*g (gravity), Σ F_x = 0.
        # Reference y_top cancels in net force (∮ n̂ ds = 0 on closed curve).
        if g != 0.0:
            # Edge-centred integration: exact Green's theorem on the actual polygon.
            # For each edge k (node k → node k+1):
            #   y_mid  = (y_k + y_{k+1}) / 2
            #   P_mid  = rho_d * g * (y_top - y_mid)   [pressure at edge midpoint]
            #   dF     = P_mid * L_k * n̂_out_k          [outward force on edge]
            # Distribute dF equally to the two endpoint nodes.
            # Net: Σ F_y = -rho_d * g * A_polygon  (exact Green's theorem identity).
            y_top  = self.x[:, 1].max()
            x_next = np.roll(self.x, -1, axis=0)        # x_{k+1}
            y_mid  = 0.5 * (self.x[:, 1] + x_next[:, 1])  # (N,) edge midpoint heights
            P_mid  = self.rho_d * g * (y_top - y_mid)   # (N,) pressure at each edge
            # n̂_out for edge k = -n_hat[k] (n_hat is INWARD; outward = negated)
            dF_edge = P_mid[:, None] * L[:, None] * (-n_hat)  # (N,2) outward per edge
            # Enforce Σ Fx = 0 exactly (removes discretisation drift)
            dF_edge[:, 0] -= dF_edge[:, 0].mean()
            # Distribute half to each endpoint node
            self.f += 0.5 * dF_edge
            self.f += 0.5 * np.roll(dF_edge, 1, axis=0)

        # ── 3. Edge elasticity ─────────────────────────────────────────────
        # For edge i (node i → node i+1):
        #   eps_i = (L_i - L0) / L0
        #   F_edge = El_t * eps_i * t_hat_i
        #   node i   += F_edge  (pulled toward i+1 if stretched)
        #   node i+1 -= F_edge
        eps = (L - self.L0) / self.L0                    # (N,)
        F_edge = self.El_t * eps[:, None] * t_hat        # (N, 2)
        self.f += F_edge
        self.f -= np.roll(F_edge, 1, axis=0)            # subtract at node i+1 (roll so index matches)
        # Note: np.roll(F_edge, 1, axis=0)[i] = F_edge[i-1], so node i gets -F_edge[i-1]
        # which is -F_edge applied at node i from edge i-1. This is correct.

        # ── 4. Bending (3-node hinge, energy-consistent) ──────────────────
        # For each node i (hinge between edges i-1 and i):
        #   theta_i = signed angle from t_hat_{i-1} to t_hat_i
        #   M_i = (EI / L0) * (theta_i - theta0_i)
        #   F_{i-1} += M_i / L0 * n_hat_{i-1}    (outward normal of left edge)
        #   F_i     -= M_i / L0 * (n_hat_{i-1} + n_hat_i)
        #   F_{i+1} += M_i / L0 * n_hat_i         (outward normal of right edge)
        t_prev = np.roll(t_hat, 1, axis=0)               # t_hat_{i-1}
        n_prev = np.roll(n_hat, 1, axis=0)               # n_hat_{i-1}

        cross = t_prev[:, 0] * t_hat[:, 1] - t_prev[:, 1] * t_hat[:, 0]
        dot   = t_prev[:, 0] * t_hat[:, 0] + t_prev[:, 1] * t_hat[:, 1]
        theta = np.arctan2(cross, dot)                   # (N,) current turning angles

        M = (self.EI / self.L0) * (theta - self.theta0)  # (N,) moments

        Mb = M / self.L0                                  # (N,)
        # Energy-consistent bending forces: F = -dU/dx.
        # Numerical verification shows the correct sign for restoring forces is:
        #   F_{i-1} -= Mb_i * n_prev_i
        #   F_i     += Mb_i * (n_prev_i + n_hat_i)
        #   F_{i+1} -= Mb_i * n_hat_i
        # (Opposite to the CAPSULE_SHELL.md spec, which has a sign error.)
        dF_im1 = -Mb[:, None] * n_prev
        dF_i   = +Mb[:, None] * (n_prev + n_hat)
        dF_ip1 = -Mb[:, None] * n_hat

        self.f += np.roll(dF_im1, -1, axis=0)  # contribution to node i from hinge i+1
        self.f += dF_i
        self.f += np.roll(dF_ip1,  1, axis=0)  # contribution to node i from hinge i-1

    def zero_forces(self):
        self.f[:] = 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Contact force computation (between particles and with primitives)
# ─────────────────────────────────────────────────────────────────────────────

def _vcp_segments_batch(points, seg_a, seg_b):
    """
    Vectorized: closest point on each segment (seg_a[j]→seg_b[j]) to each point in points.

    points: (M, 2)
    seg_a:  (N, 2)
    seg_b:  (N, 2)

    Returns:
        cp: (M, N, 2) — closest points
        t:  (M, N)    — parameter in [0, 1] along each segment
    """
    ab   = seg_b - seg_a                              # (N, 2)
    ab2  = np.sum(ab * ab, axis=1)                   # (N,)
    pa   = points[:, None, :] - seg_a[None, :, :]   # (M, N, 2)
    dots = np.einsum('mnk,nk->mn', pa, ab)           # (M, N)
    safe = ab2 > 1e-30
    t    = np.where(safe[None, :],
                     np.clip(dots / np.where(safe, ab2, 1.0)[None, :], 0.0, 1.0),
                     np.zeros_like(dots))             # (M, N)
    cp   = seg_a[None, :, :] + t[:, :, None] * ab[None, :, :]  # (M, N, 2)
    return cp, t


def _capsule_capsule_forces(pA, pB):
    """
    Compute capsule–capsule contact forces (vectorized, sparse active-contact).
    All edges of pA vs all edges of pB, 2-point Gauss quadrature.
    Forces accumulated into pA.f and pB.f.
    """
    NA = pA.N;  NB = pB.N
    xA = pA.x;  xA_next = np.roll(xA, -1, axis=0)   # (NA, 2)
    xB = pB.x;  xB_next = np.roll(xB, -1, axis=0)   # (NB, 2)

    contact_r = pA.r_c + pB.r_c
    k_c       = pA.k_c
    L0A       = pA.L0

    # Broad-phase: check bounding box of the two particles
    cA = xA.mean(axis=0);  cB = xB.mean(axis=0)
    if np.linalg.norm(cA - cB) > 2.5 * (pA.R0 + pB.R0):
        return

    # Gauss quadrature points on all edges of A: (NA*2, 2)
    s   = _GAUSS2_S
    w   = _GAUSS2_W
    xq_flat = ((1 - s)[np.newaxis, :, np.newaxis] * xA[:, np.newaxis, :]
               + s[np.newaxis, :, np.newaxis] * xA_next[:, np.newaxis, :]).reshape(NA * 2, 2)

    # Closest point on each edge of B to each quadrature point: (M, NB, 2), (M, NB)
    cp, t_B = _vcp_segments_batch(xq_flat, xB, xB_next)

    diff = xq_flat[:, np.newaxis, :] - cp              # (M, NB, 2)
    d2   = np.sum(diff * diff, axis=2)                  # (M, NB)
    gap2 = d2 - contact_r * contact_r                   # squared gap criterion (faster)
    # Active contacts (gap < 0): where d2 < contact_r^2
    active = gap2 < 0.0

    if not np.any(active):
        return

    # Compute only for active pairs
    m_act, j_act = np.nonzero(active)                   # sparse indices
    d_act   = np.sqrt(d2[m_act, j_act])                 # (K,)
    gap_act = d_act - contact_r                          # (K,) — negative
    diff_act = diff[m_act, j_act]                       # (K, 2)

    safe_d = np.where(d_act > 1e-15, d_act, 1.0)
    n_hat_act = diff_act / safe_d[:, np.newaxis]        # (K, 2)

    # Force magnitudes: weight depends on which Gauss point (m % 2)
    w_m_act = w[m_act % 2]                              # (K,)
    F_mag   = k_c * (-gap_act) * L0A * w_m_act          # (K,)
    F_vec   = F_mag[:, np.newaxis] * n_hat_act           # (K, 2)

    # A-side shape functions
    i_act   = m_act // 2                                 # edge index on A (K,)
    s_m_act = s[m_act % 2]                              # (K,)
    np.add.at(pA.f, i_act,            (1.0 - s_m_act)[:, np.newaxis] * F_vec)
    np.add.at(pA.f, (i_act + 1) % NA, s_m_act[:, np.newaxis]         * F_vec)

    # B-side shape functions (equal and opposite)
    t_B_act = t_B[m_act, j_act]                         # (K,)
    np.add.at(pB.f, j_act,            -(1.0 - t_B_act)[:, np.newaxis] * F_vec)
    np.add.at(pB.f, (j_act + 1) % NB, -t_B_act[:, np.newaxis]         * F_vec)


def _capsule_primitive_forces(particle, primitives):
    """
    Compute capsule–primitive contact forces.
    Uses gap_and_normal_batch() if available on the primitive, else falls back
    to per-point gap_and_normal() calls.
    """
    if not primitives:
        return

    N      = particle.N
    x      = particle.x
    x_next = np.roll(x, -1, axis=0)
    r_c    = particle.r_c
    k_c    = particle.k_c
    L0     = particle.L0
    s      = _GAUSS2_S
    w      = _GAUSS2_W

    # Gauss quadrature points: (N*2, 2)
    xq_flat = ((1 - s)[np.newaxis, :, np.newaxis] * x[:, np.newaxis, :]
               + s[np.newaxis, :, np.newaxis] * x_next[:, np.newaxis, :]).reshape(N * 2, 2)

    i_idx = np.repeat(np.arange(N), 2)                 # edge index for each quad pt
    s_m   = np.tile(s, N)                              # Gauss param for each quad pt
    w_m   = np.tile(w, N)                              # weight for each quad pt

    for prim in primitives:
        # Use batch interface if available, else fall back to per-point loop
        if hasattr(prim, 'gap_and_normal_batch'):
            gaps, normals = prim.gap_and_normal_batch(xq_flat)
        else:
            gaps    = np.empty(N * 2)
            normals = np.empty((N * 2, 2))
            for m, xq_m in enumerate(xq_flat):
                g, n = prim.gap_and_normal(xq_m)
                gaps[m]    = g
                normals[m] = n

        gap_capsule = gaps - r_c - getattr(prim, 'r_c', 0.0)
        act = gap_capsule < 0.0
        if not np.any(act):
            continue

        F_mag = k_c * (-gap_capsule[act]) * L0 * w_m[act]    # (K,)
        F_vec = F_mag[:, np.newaxis] * normals[act]            # (K, 2)
        i_a   = i_idx[act]
        s_a   = s_m[act]

        np.add.at(particle.f, i_a,            (1.0 - s_a)[:, np.newaxis] * F_vec)
        np.add.at(particle.f, (i_a + 1) % N,  s_a[:, np.newaxis]         * F_vec)


# ─────────────────────────────────────────────────────────────────────────────
# CapsuleSim — integrator and simulation manager
# ─────────────────────────────────────────────────────────────────────────────

class CapsuleSim:
    """
    Manages a collection of CapsuleParticle objects and rigid primitives.
    Advances state with explicit Verlet integration.

    Parameters
    ----------
    particles : list of CapsuleParticle
    primitives : list of contact primitives (LineSegment, Arc, Polygon, ...)
    g : float
        Gravitational acceleration (0 = no gravity).
    """

    def __init__(self, particles, primitives=None, g=0.0):
        self.particles  = particles
        self.primitives = primitives or []
        self.g          = g
        self.t          = 0.0

    def compute_forces(self):
        """Compute all forces on all particles (resets f first)."""
        for p in self.particles:
            p.zero_forces()

        # Internal forces (per particle)
        for p in self.particles:
            p.accumulate_internal_forces(g=self.g)

        # Particle–particle contact
        for i in range(len(self.particles)):
            for j in range(i + 1, len(self.particles)):
                _capsule_capsule_forces(self.particles[i], self.particles[j])

        # Particle–primitive contact
        if self.primitives:
            for p in self.particles:
                _capsule_primitive_forces(p, self.primitives)

    def step(self, dt, alpha_damp=0.0):
        """
        One explicit Verlet step with optional linear velocity damping.

        F_total_i = F_internal_i + F_contact_i - alpha_damp * v_i * m_i
        a_i = F_total_i / m_i
        v_i += a_i * dt
        x_i += v_i * dt

        alpha_damp is a dimensionless damping coefficient: F_damp = -alpha_damp * m * v.
        Large alpha_damp (>> 1/dt) → over-damped / quasi-static limit.
        """
        self.compute_forces()

        for p in self.particles:
            m   = p.m_node
            if alpha_damp > 0.0:
                p.f -= alpha_damp * m * p.v
            a   = p.f / m
            p.v += a * dt
            p.x += p.v * dt

        self.t += dt

    def step_rb(self, dt, alpha_damp=0.0, f_ext=None, contact_log=None):
        """
        Rigid-body decomposition integrator.

        Splits each particle's motion into:
          • Rigid part (x_cm, theta) — integrated with undamped Verlet,
            correct solid-disk mass M = ρ_f π R₀² and I = ½ M R₀².
          • Elastic part (u_i in body frame) — integrated with damped Verlet;
            only deformation-mode velocity u_dot_i is damped, never v_cm or ω.

        This fixes two errors in the legacy step():
          E1 — node-only integration uses shell inertia I = m R₀² (factor-2 error)
          E2 — nodal damping spuriously damps rigid-body translation

        Conservation properties:
          • Total linear momentum changes only due to net external (contact + wall) force.
          • Between contacts: p_cm is constant to machine precision.
          • Energy dissipated = α_d · m_node · Σ|u_dot_i|² · dt per step (elastic modes only).

        Algorithm
        ---------
        1. Compute all forces (internal elastic + contact) → p.f
        2. Rigid-body resultants: F = Σ f_i, T = Σ r_i × f_i  (elastic terms cancel)
        3. Undamped Verlet for (v_cm, x_cm) and (omega, theta)
        4. Deformation force = deviatoric total force: f_i^dev = f_i − F/N
           (zero-sum by construction → preserves Σu_i = 0)
        5. Damped Verlet for (u_dot_i, u_i) driven by f_i^dev in body frame
        6. Reconstruct perimeter: x_i = x_cm + R(theta) @ (X_ref_i + u_i)

        Why f_i^dev = f_i − F/N drives deformation:
          • Elastic-only (f_contact = 0): f_i^dev = f_i^el  (correct; elastic forces sum to 0)
          • During contact: deviatoric contact force (bottom nodes pushed up, rest react)
            drives the shape deformation that stores and then releases elastic energy.
          • F/N is the per-node rigid-body share; subtracting it keeps Σu_i = 0 invariant.
        """
        # ── 1. Elastic (internal) forces — saved separately ───────────────────
        # This is required to isolate contact forces for the rigid-body equation.
        # The discrete fluid-pressure force does NOT sum to exactly zero for a
        # non-circular shape (it uses normalized average normals, not exact area
        # gradients). Keeping f_el separate lets us compute F = Σf_contact = 0
        # exactly during free flight, ensuring machine-precision momentum conservation.
        for p in self.particles:
            p.zero_forces()
            p.accumulate_internal_forces(g=self.g)
        for p in self.particles:
            # Save physical elastic forces.  For EmulsionParticle the tangential
            # regularisation pseudo-force (_f_reg) is stored separately and must
            # be excluded from the rigid-body resultants F and T.
            f_el = p.f.copy()
            if hasattr(p, '_f_reg'):
                f_el -= p._f_reg
            p.f_el = f_el

        # ── 2. Contact forces (added on top of elastic) ───────────────────────
        for i in range(len(self.particles)):
            for j in range(i + 1, len(self.particles)):
                _capsule_capsule_forces(self.particles[i], self.particles[j])
        if self.primitives:
            for p in self.particles:
                _capsule_primitive_forces(p, self.primitives)

        # ── 2b. External forces (validation/testing only) ─────────────────────
        # Injected AFTER f_el is saved → appear in F = F_total − F_el exactly.
        # f_ext: dict {particle_index: (N,2) force array for one step}
        if f_ext is not None:
            for pidx, farray in f_ext.items():
                self.particles[pidx].f += np.asarray(farray, dtype=float)

        # ── 3–6. Per-particle rigid + elastic integration ─────────────────────
        for p_idx, p in enumerate(self.particles):
            # ── 3. Rigid-body resultants from CONTACT forces only ──────────────
            # f_contact = f_total − f_el − f_reg; analytically Σf_el = 0 and Σr×f_el = 0
            # (elastic energy is rigid-motion invariant), but numerically there is
            # a small imbalance in the discrete pressure force for non-circular
            # shapes.  Computing F = Σf_total − Σf_el = Σf_contact cancels this
            # error exactly: F = 0 during free flight regardless of shape distortion.
            #
            # NOTE: hydrostatic gravity is stored in f_el (via accumulate_internal_forces)
            # and is therefore subtracted here.  We add it back as an explicit rigid-body
            # body force so that the particle falls under gravity correctly.
            #
            # NOTE: for EmulsionParticle the tangential regularisation pseudo-force
            # (_f_reg) has a non-zero sum for deformed shapes.  It was already
            # excluded from f_el (f_el = p.f_before − _f_reg), so it appears in
            # F = F_total − F_el as a spurious rigid-body force.  We subtract it
            # explicitly here to keep F = Σf_contact exactly.
            F_el     = p.f_el.sum(axis=0)                   # Σf_el (physical, no reg)
            F_total  = p.f.sum(axis=0)                      # Σ(f_el + f_reg + f_contact + f_ext)
            F        = F_total - F_el                        # Σf_reg + Σf_contact + Σf_ext
            if hasattr(p, '_f_reg'):
                F -= p._f_reg.sum(axis=0)                   # subtract reg pseudo-force → F = Σf_contact
            # Re-add gravitational body force (was in f_el via hydrostatic; cancel = -M*g)
            if self.g != 0.0:
                F[1] -= p.M_disk * self.g   # F_y += -M*g  (downward)

            r        = p.x - p.x_cm                         # (N,2) moment arms
            T_total  = float(np.sum(r[:, 0] * p.f[:, 1]
                                  - r[:, 1] * p.f[:, 0]))
            T_el     = float(np.sum(r[:, 0] * p.f_el[:, 1]
                                  - r[:, 1] * p.f_el[:, 0]))
            T        = T_total - T_el                        # torque from reg+contact+ext
            if hasattr(p, '_f_reg'):
                T -= float(np.sum(r[:, 0] * p._f_reg[:, 1]
                                - r[:, 1] * p._f_reg[:, 0]))  # subtract reg torque

            # ── Optional contact-force recording (for validation) ──────────────
            if contact_log is not None:
                contact_log.append(dict(
                    pidx=p_idx,
                    t=self.t,
                    F=F.copy(),          # net contact (+ext) force on rigid body
                    T=T,                 # net torque on rigid body
                    f_contact=(p.f - p.f_el).copy(),  # per-node contact force (N,2)
                    r=r.copy(),          # moment arms (N,2)
                ))

            # ── 4. Undamped Verlet for rigid body ──────────────────────────────
            p.v_cm  += (F / p.M_disk) * dt
            p.x_cm  += p.v_cm * dt
            p.omega += (T / p.I_disk) * dt
            p.theta += p.omega * dt

            # ── 5. Deformation force — deviatoric of physical forces ──────────────
            # f_i^dev = (f_i − f_reg_i) − (F_total − F_reg)/N  (zero-mean)
            # Excludes the regularisation pseudo-force from the elastic driver so
            # only physical forces (elastic restoring + contact) deform the shape.
            # Subtracting the corrected mean keeps Σf_dev = 0 → Σu_i stays zero.
            if hasattr(p, '_f_reg'):
                F_reg_sum  = p._f_reg.sum(axis=0)
                f_phys     = p.f - p._f_reg                 # (N,2) physical forces only
                F_phys     = F_total - F_reg_sum             # Σ physical forces
            else:
                f_phys     = p.f
                F_phys     = F_total
            f_dev = f_phys - F_phys[None, :] / p.N          # (N,2), zero-sum

            # Rotate to body frame: R(-θ) = [[c, s], [-s, c]]
            c = np.cos(p.theta);  s = np.sin(p.theta)
            f_dev_body = np.column_stack([
                 c * f_dev[:, 0] + s * f_dev[:, 1],
                -s * f_dev[:, 0] + c * f_dev[:, 1],
            ])                                               # (N,2)

            # Damping acts only on u_dot — elastic velocity, not rigid-body velocity.
            accel_el = f_dev_body / p.m_node - alpha_damp * p.u_dot
            p.u_dot += accel_el * dt
            p.u     += p.u_dot * dt

            # ── 6. Reconstruct perimeter ──────────────────────────────────────
            # x_i = x_cm + R(+θ)(X_ref_i + u_i),   R(θ) = [[c, -s], [s, c]]
            body = p.X_ref + p.u                             # (N,2)
            p.x = p.x_cm + np.column_stack([
                c * body[:, 0] - s * body[:, 1],
                s * body[:, 0] + c * body[:, 1],
            ])

        self.t += dt

    def step_rb_to(self, t_end, dt, alpha_damp=0.0,
                   callback=None, callback_interval=1):
        """Integrate with step_rb() from current time to t_end."""
        n = 0
        while self.t < t_end - 0.5 * dt:
            self.step_rb(dt, alpha_damp=alpha_damp)
            n += 1
            if callback is not None and n % callback_interval == 0:
                callback(self, n)
        return n

    def step_to(self, t_end, dt, alpha_damp=0.0, callback=None, callback_interval=1):
        """
        Integrate from current time to t_end.
        callback(sim, step_idx) called every callback_interval steps if provided.
        Returns number of steps taken.
        """
        n = 0
        while self.t < t_end - 0.5 * dt:
            self.step(dt, alpha_damp=alpha_damp)
            n += 1
            if callback is not None and n % callback_interval == 0:
                callback(self, n)
        return n

    # ── diagnostics ───────────────────────────────────────────────────────────

    def kinetic_energy(self):
        """Total kinetic energy of all nodes."""
        ke = 0.0
        for p in self.particles:
            ke += 0.5 * p.m_node * np.sum(p.v ** 2)
        return ke

    def potential_energy(self):
        """
        Total elastic potential energy:
        E_fluid + E_edge + E_bend (summed over all particles).
        Contact energy not included (penalty formulation — not a true potential).
        """
        pe = 0.0
        for p in self.particles:
            # Fluid (area): U = (K_fluid/2) * (A-A0)^2 / A0
            A   = p.area()
            dA  = (A - p.A0) / p.A0
            pe += 0.5 * p.K_fluid * p.A0 * dA ** 2  # symmetric in dA, sign doesn't matter

            # Edge elasticity
            t_hat, n_hat, L = p.edge_tangents_normals()
            eps = (L - p.L0) / p.L0
            pe += 0.5 * p.El_t * p.L0 * np.sum(eps ** 2)

            # Bending
            t_prev = np.roll(t_hat, 1, axis=0)
            cross = t_prev[:, 0] * t_hat[:, 1] - t_prev[:, 1] * t_hat[:, 0]
            dot   = t_prev[:, 0] * t_hat[:, 0] + t_prev[:, 1] * t_hat[:, 1]
            theta = np.arctan2(cross, dot)
            dtheta = theta - p.theta0
            pe += 0.5 * (p.EI / p.L0) * np.sum(dtheta ** 2)

        return pe

    def total_energy(self):
        return self.kinetic_energy() + self.potential_energy()

    def net_force_on_particle(self, p):
        """Sum of all nodal forces on particle (should be zero for internal forces)."""
        self.compute_forces()
        return p.f.sum(axis=0)

    def wall_forces(self):
        """
        Compute net force exerted by each primitive on the particle system.
        Returns list of (2,) arrays, one per primitive.
        """
        # Re-compute forces but track per-primitive contribution
        # We do this by temporarily replacing primitives and computing delta-f
        for p in self.particles:
            p.zero_forces()
            p.accumulate_internal_forces(g=self.g)
        for i in range(len(self.particles)):
            for j in range(i + 1, len(self.particles)):
                _capsule_capsule_forces(self.particles[i], self.particles[j])

        result = []
        for prim in self.primitives:
            for p in self.particles:
                p.f[:] = 0.0
            for p in self.particles:
                _capsule_primitive_forces(p, [prim])
            total = sum(p.f.sum(axis=0) for p in self.particles)
            result.append(total.copy())
        return result

    # ── rigid-body diagnostics (for step_rb) ──────────────────────────────────

    def kinetic_energy_rb(self):
        """
        Rigid-body kinetic energy: ½ M v_cm² + ½ I ω²  (per particle, summed).
        Use with step_rb(); meaningless after legacy step().
        """
        ke = 0.0
        for p in self.particles:
            ke += 0.5 * p.M_disk * np.dot(p.v_cm, p.v_cm)
            ke += 0.5 * p.I_disk * p.omega ** 2
        return ke

    def elastic_kinetic_energy(self):
        """
        Kinetic energy of elastic (deformation) modes: ½ m_node Σ|u_dot_i|².
        Use with step_rb().
        """
        ke = 0.0
        for p in self.particles:
            ke += 0.5 * p.m_node * np.sum(p.u_dot ** 2)
        return ke

    def total_energy_rb(self):
        """Total energy under rigid-body integrator: KE_rigid + KE_elastic + PE."""
        return self.kinetic_energy_rb() + self.elastic_kinetic_energy() + self.potential_energy()

    def momentum(self):
        """Total linear momentum Σ M_disk v_cm_i (rigid-body integrator)."""
        return sum(p.M_disk * p.v_cm for p in self.particles)

    def angular_momentum(self):
        """Total angular momentum Σ (I_disk ω_i + M_disk |x_cm_i| cross v_cm_i)."""
        L = 0.0
        for p in self.particles:
            L += p.I_disk * p.omega
            L += float(p.x_cm[0] * p.v_cm[1] - p.x_cm[1] * p.v_cm[0]) * p.M_disk
        return L

    def dissipation_rate_rb(self, alpha_damp):
        """
        Instantaneous power dissipated by elastic damping: α_d m_node Σ|u_dot_i|².
        Integrating this over time gives the total energy lost to dissipation.
        """
        rate = 0.0
        for p in self.particles:
            rate += alpha_damp * p.m_node * np.sum(p.u_dot ** 2)
        return rate

    def elastic_displacement_rms(self):
        """RMS elastic displacement Σ|u_i|/N per particle (shape distortion measure)."""
        return [float(np.sqrt(np.mean(np.sum(p.u ** 2, axis=1)))) for p in self.particles]

    # ── stability estimate ─────────────────────────────────────────────────────

    def estimate_dt_max(self):
        """
        Estimate stable timestep from the three mode frequencies.
        Returns (dt_max, dominant_mode_name).

        Modes:
          ω_breath  = sqrt(2 K_fluid / (rho_d t R0))
          c_edge    = sqrt(El_t / (rho_d L0))  →  omega_edge = N * c_edge / R0
          ω_contact = sqrt(k_c / m_node)
        """
        dt_mins = []
        names   = []
        for p in self.particles:
            t_thick = p.tau * p.R0
            rho_d   = p.rho_d

            # Breathing mode (uses K_area, which may differ from K_fluid)
            om_breath = np.sqrt(2.0 * p.K_area / (rho_d * t_thick * p.R0))
            dt_mins.append(2.0 / om_breath)
            names.append(f"breath(p{id(p)%1000})")

            # Edge wave
            c_edge   = np.sqrt(p.El_t / (rho_d * p.L0))
            om_edge  = p.N * c_edge / p.R0
            dt_mins.append(2.0 / om_edge)
            names.append(f"edge(p{id(p)%1000})")

            # Contact
            om_contact = np.sqrt(p.k_c / p.m_node)
            dt_mins.append(2.0 / om_contact)
            names.append(f"contact(p{id(p)%1000})")

        idx    = int(np.argmin(dt_mins))
        dt_max = dt_mins[idx]
        return dt_max, names[idx]
