"""
EmulsionParticle — a 2D fluid droplet for the emulsion DEM model.

Inherits CapsuleParticle but replaces membrane + bending elasticity with:
  1. Line tension  γ L  (total perimeter length energy)
  2. Tangential regularization pseudo-force (keeps nodes evenly spaced;
     applied inside the shell integrator ONLY — never added to rigid-body sums)

Area incompressibility and internal hydrostatic gravity are inherited unchanged.

Dimensionless parameters (natural scales: γ force, R0 length):
  κ      = γ / (R0 * K_area)          pre-compression ratio (must be < 1)
  Bo     = rho_d * g * R0² / γ        Bond number
  C̃      = C * R0 / γ                 dimensionless contact stiffness
  k̃_reg  = k_reg / γ                  numerical regularization (not physical)

Derived:  A0 = π R0² / (1 − κ)   (equilibrium circular droplet at radius R0)
          m  = rho_d * π R0²      (total fluid mass, distributed evenly to nodes)
"""

import numpy as np
from src.simulation.capsule_shell import CapsuleParticle


class EmulsionParticle(CapsuleParticle):
    """
    Parameters
    ----------
    N : int       Number of perimeter nodes (default 120).
    R0 : float    Equilibrium radius.
    gamma : float Line tension (force scale; sets γ).
    K_area : float Area incompressibility stiffness.
    C : float     Contact hardness (same role as CapsuleParticle).
    rho_d : float Internal fluid density.
    k_reg : float Tangential regularization stiffness (numerical, not physical).
                  Recommended: k_reg = 10 * gamma.  Set 0 to disable.
    center : array-like (2,)
    """

    def __init__(self, N=120, R0=1.0, gamma=1.0, K_area=5.0,
                 C=500.0, rho_d=1.0, k_reg=None, center=(0.0, 0.0)):

        # κ = γ / (R0 * K_area) — must be < 1 for a valid equilibrium
        kappa = gamma / (R0 * K_area)
        if kappa >= 1.0:
            raise ValueError(
                f"κ = γ/(R0·K_area) = {kappa:.3f} ≥ 1: no equilibrium exists. "
                f"Increase K_area or decrease gamma.")

        # Initialise parent with tau→0, S→0 so El_t=EI=0
        # Use K_fluid=K_area so the area term is set correctly.
        super().__init__(N=N, R0=R0, tau=1e-6, S=0.0, C=C,
                         rho_d=rho_d, K_fluid=K_area, K_area=K_area,
                         center=center)

        # Zero out elastic shell terms explicitly (tau→0 already makes them ~0)
        self.El_t   = 0.0
        self.EI     = 0.0
        self.gamma  = float(gamma)
        self.kappa  = kappa

        # Reference values for R0-scaling recompute (γ ∝ R0, k_c ∝ R0^{-1})
        self._gamma0 = float(gamma) / float(R0)  # γ per unit R0 (= γ at R0=1)

        # Override A0: use the actual polygon area (set by parent as self.area()),
        # scaled by 1/(1-κ) so the polygon equilibrium condition is satisfied exactly:
        #   γ/R0 = K_area * (A0 - A_polygon) / A0  →  A0 = A_polygon / (1 - κ)
        # This avoids a drift when A_polygon < π R0² (polygon vs circle correction).
        A_polygon = self.A0   # parent set this to self.area() at initialisation
        self.A0   = A_polygon / (1.0 - kappa)

        # Override node mass: fluid mass distributed evenly (no shell thickness)
        m_total     = rho_d * np.pi * R0**2
        self.m_node = m_total / N

        # Tangential regularization stiffness
        self.k_reg  = float(10.0 * gamma if k_reg is None else k_reg)

        # Store kappa for diagnostics
        self._kappa = kappa

    # ------------------------------------------------------------------
    def _recompute_params(self):
        """
        Emulsion override: γ ∝ R0^{+1}, k_c ∝ R0^{-1}, m_node ∝ R0^{+2}.
        K_area and C are constants (dimensionless κ and C̃ preserved).
        Assumes node positions and L0 already reflect the new R0.
        """
        R0 = self._R0
        self.gamma  = self._gamma0 * R0
        self.k_c    = self.C * self.K_fluid / R0
        m_total     = self.rho_d * np.pi * R0 ** 2
        self.m_node = m_total / self.N
        self.M_disk = m_total
        self.I_disk = 0.5 * m_total * R0 ** 2
        # A0 tracks current geometry: re-derive from actual polygon area
        self.A0     = self.area() / (1.0 - self.kappa)
        # k_reg tracks gamma
        self.k_reg  = 10.0 * self.gamma
        # Elastic terms stay zero for emulsion
        self.El_t = 0.0
        self.EI   = 0.0

    # ------------------------------------------------------------------
    def accumulate_internal_forces(self, g=0.0):
        """
        Shell forces for the emulsion model:
          1. Area pressure  (inherited logic, using overridden A0)
          2. Internal hydrostatic gravity (inherited, edge-centred)
          3. Line tension  γ L
          4. Tangential regularisation  (NOT added to rigid-body force sum)
        """
        t_hat, n_hat, L = self.edge_tangents_normals()
        n_node = self.node_normals(n_hat)

        # ── 1. Area pressure (same as CapsuleParticle) ─────────────────
        A  = self.area()
        P  = self.K_area * (self.A0 - A) / self.A0
        Fp = P * self.L0 * n_node          # n_node is INWARD
        self.f -= Fp                        # subtract → outward

        # ── 2. Internal hydrostatic gravity (inherited) ─────────────────
        if g != 0.0:
            y_top  = self.x[:, 1].max()
            x_next = np.roll(self.x, -1, axis=0)
            y_mid  = 0.5 * (self.x[:, 1] + x_next[:, 1])
            P_mid  = self.rho_d * g * (y_top - y_mid)
            dF_edge = P_mid[:, None] * L[:, None] * (-n_hat)
            dF_edge[:, 0] -= dF_edge[:, 0].mean()
            self.f += 0.5 * dF_edge
            self.f += 0.5 * np.roll(dF_edge, 1, axis=0)

        # ── 3. Line tension  γ L ───────────────────────────────────────
        # Force on node k: F_k = -γ (t̂_{k,k+1} − t̂_{k−1,k})
        # = γ (t̂_{k−1,k} − t̂_{k,k+1})
        # This is the discrete curvature force: pulls each node toward
        # its neighbours, producing inward pressure γ/R on a circle.
        if self.gamma != 0.0:
            # F_k = -∂(γL)/∂x_k = γ (t̂_{k,k+1} − t̂_{k−1,k})
            # Points INWARD for a convex polygon (curvature pressure γ/R inward).
            t_next = t_hat                       # t̂_{k, k+1}
            t_prev = np.roll(t_hat, 1, axis=0)  # t̂_{k-1, k}
            self.f += self.gamma * (t_next - t_prev)

        # ── 4. Tangential regularisation (pseudo-force, shell-only) ────
        # F_reg,k = k_reg * (ℓ_k+ − ℓ_k−) / 2 * t̂_k
        # Drives nodes toward equal arc-length spacing.
        # Excluded from rigid-body force/torque sums by design.
        if self.k_reg != 0.0:
            L_plus  = L                             # ℓ_k+ = length of edge (k → k+1)
            L_minus = np.roll(L, 1)                 # ℓ_k- = length of edge (k-1 → k)
            t_node  = 0.5 * (t_hat + np.roll(t_hat, 1, axis=0))  # tangent at node k
            t_mag   = np.linalg.norm(t_node, axis=1, keepdims=True)
            t_node  = t_node / np.maximum(t_mag, 1e-30)
            F_reg   = (self.k_reg * 0.5 * (L_plus - L_minus))[:, None] * t_node
            self.f += F_reg
            # Store separately so step_rb can subtract it from the RB sum
            self._f_reg = F_reg.copy()
        else:
            self._f_reg = np.zeros_like(self.f)

    def rb_force_torque(self, x_cm):
        """
        Return (F_total, T_total) excluding the regularisation pseudo-force.
        Called by the rigid-body integrator instead of summing self.f directly.
        """
        f_phys = self.f - self._f_reg   # physical forces only
        F = f_phys.sum(axis=0)
        r = self.x - x_cm
        T = float((r[:, 0] * f_phys[:, 1] - r[:, 1] * f_phys[:, 0]).sum())
        return F, T
