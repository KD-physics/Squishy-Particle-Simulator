"""
rcp_utils.py — Bridge between rcpgenerator and the EPD simulation engine.

Workflow
--------
1. Use rcpgenerator to obtain particle positions at the 2D jamming point φ_J.
2. Scale to EPD units (R0=1, r_c = 2π/N) so the physical packing fraction matches.
3. Optionally rescale the box (and particle CMs affinely) to a target φ != φ_J.

Key insight: rcpgenerator outputs positions in arbitrary units with diameter d_rcp.
The EPD effective radius is R_eff = R0 + r_c.  We rescale so 2*R_eff = d_rcp_scaled,
i.e., the physical touching distance in EPD matches the touching distance in rcpgenerator.
"""

import numpy as np
from src.simulation.capsule_shell import CapsuleParticle


# ── public API ────────────────────────────────────────────────────────────────

def rcp_seed(N_particles, N_nodes=60, polydispersity=0.0, seed=42,
             walls=None, q=2.0, tau_b=0.2, S=1.0, C_scale=1.0, verbose=True):
    """
    Generate a jammed 2D packing via rcpgenerator and return EPD particles.

    Parameters
    ----------
    N_particles    : int   — number of particles
    N_nodes        : int   — perimeter nodes per EPD particle (default 60)
    polydispersity : float — radius polydispersity (std/mean of diameter distribution)
                             0 = monodisperse; rcpgenerator 'bi' type uses ±delta
    seed           : int   — random seed
    walls          : list  — [wx, wy] where 0=periodic, 1=hard wall (default [0,0])
    q              : float — squishiness parameter (K_area / El_t_base)
    tau_b          : float — bending working point (default 0.2)
    S              : float — force scale (default 1.0)
    C_scale        : float — multiplier on the canonical C = 3000·S·(1+q).
                             Use C_scale = 1/n_contacts to compensate for multi-contact
                             loading in dense packings (default 1.0 = two-particle calibration)
    verbose        : bool

    Returns
    -------
    particles : list of CapsuleParticle at EPD units (R0≈1, positions scaled to φ_J)
    Lx, Ly    : float — periodic box dimensions in EPD units
    phi_J     : float — measured jamming packing fraction from rcpgenerator
    """
    from rcpgenerator import Packing

    if walls is None:
        walls = [0, 0]

    # ── build rcpgenerator packing ────────────────────────────────────────────
    if polydispersity > 0:
        # Gaussian: sigma/mean = polydispersity (e.g. 0.05 = 5% PDI)
        dist = {'type': 'gaussian', 'd': 1.0, 'sigma': polydispersity}
    else:
        dist = {'type': 'mono', 'd': 1.0}

    # Initial box: large enough to start at low φ (~0.3)
    # rcpgenerator will compress during pack() to reach φ_J
    box_init = np.sqrt(N_particles * np.pi * 0.25 / 0.3)   # side of square at φ=0.3
    packing = Packing(N=N_particles, Ndim=2,
                      box=[box_init, box_init],
                      walls=walls,
                      dist=dist,
                      seed=seed)
    packing.pack(verbose=verbose)

    phi_J   = packing.phi_final
    d_rcp   = np.mean(packing.diameters)   # mean diameter in rcp units
    pos_rcp = np.array(packing.positions)  # (N_p, 2)
    dia_rcp = np.array(packing.diameters)  # (N_p,)
    box_rcp = packing.box                  # [Lx_rcp, Ly_rcp]

    if verbose:
        print(f"  rcpgenerator: N={N_particles}, φ_J={phi_J:.4f}, "
              f"box={box_rcp[0]:.3f}×{box_rcp[1]:.3f}, d_rcp={d_rcp:.4f}")

    # ── scale to EPD units ────────────────────────────────────────────────────
    # Just-touching condition: d_center = 2*R0 + 2*r_c  where r_c = L0 = 2*R0*sin(π/N).
    # Normalise so <R0> = 1.  Then d_center_touch = 2*(1 + 2*sin(π/N)) ≡ 2*R_eff.
    # This gives zero contact-force penetration at the rcpgenerator φ_J configuration.
    R_eff = 1.0 + 2.0 * np.sin(np.pi / N_nodes)   # effective radius for R0=1
    scale = 2.0 * R_eff / d_rcp                    # position / box scale

    pos_epd = pos_rcp * scale                       # (N_p, 2)
    Lx_epd  = box_rcp[0] * scale
    Ly_epd  = box_rcp[1] * scale

    # R0 normalised: relative sizes preserved, mean = 1 exactly.
    # Particles are SMALLER than the rcp discs so the r_c halos just touch at φ_J.
    R0_arr = dia_rcp / np.mean(dia_rcp)             # (N_p,)  <R0> = 1

    if verbose:
        phi_check = N_particles * np.pi * R_eff**2 / (Lx_epd * Ly_epd)
        print(f"  EPD scale: R_eff={R_eff:.4f}, scale={scale:.4f}, "
              f"box={Lx_epd:.3f}×{Ly_epd:.3f}, φ_check={phi_check:.4f}")

    # ── build CapsuleParticle list ────────────────────────────────────────────
    TAU     = np.sqrt(12.0 * tau_b)
    K_area  = q * (12.0 * S / TAU**2)

    particles = []
    for i in range(N_particles):
        R0_i = float(R0_arr[i])        # ≈ 1.0 for monodisperse
        C_i  = 3000.0 * S * (1.0 + q) * C_scale
        p = CapsuleParticle(
            N=N_nodes, R0=R0_i, tau=TAU / np.sqrt(12),
            S=S, C=C_i, K_area=K_area,
            center=(float(pos_epd[i, 0]), float(pos_epd[i, 1])),
        )
        particles.append(p)

    return particles, Lx_epd, Ly_epd, phi_J


def scale_packing(particles, Lx, Ly, phi_target, phi_current=None, center=None):
    """
    Affinely rescale all particle CM positions and the box to reach phi_target.

    Does NOT change R0 or any internal particle parameters — only x_cm moves.
    phi_target > phi_current → compression (box shrinks, particles move inward).
    phi_target < phi_current → dilation  (box grows,  particles spread out).

    Parameters
    ----------
    particles    : list of CapsuleParticle
    Lx, Ly       : float — current box dimensions
    phi_target   : float — desired packing fraction
    phi_current  : float — current packing fraction (computed if None)
    center       : (2,) — scaling origin (default: box centre)

    Returns
    -------
    particles    : same list, CMs updated in place
    Lx_new, Ly_new : float — new box dimensions
    """
    if phi_current is None:
        R_eff = np.mean([p.R0 + 2.0 * np.pi * p.R0 / p.N for p in particles])
        phi_current = len(particles) * np.pi * R_eff**2 / (Lx * Ly)

    scale   = np.sqrt(phi_current / phi_target)   # < 1 → compress
    Lx_new  = Lx * scale
    Ly_new  = Ly * scale

    if center is None:
        cx, cy = Lx / 2.0, Ly / 2.0
    else:
        cx, cy = center

    # Offset to shift scaled box to [0, Lx_new] × [0, Ly_new]
    x_off = cx - Lx_new / 2.0
    y_off = cy - Ly_new / 2.0

    for p in particles:
        p.x_cm[0] = cx + scale * (p.x_cm[0] - cx) - x_off
        p.x_cm[1] = cy + scale * (p.x_cm[1] - cy) - y_off
        # Rebuild x from rigid-body frame (CM + rotation + reference shape + elastic displacement)
        c = np.cos(p.theta);  s = np.sin(p.theta)
        body = p.X_ref + p.u                          # (N, 2)
        p.x[:, 0] = p.x_cm[0] + c * body[:, 0] - s * body[:, 1]
        p.x[:, 1] = p.x_cm[1] + s * body[:, 0] + c * body[:, 1]

    return particles, Lx_new, Ly_new


def identify_layers(particles, Ly, n_layers=1, layer_frac=0.20):
    """
    Identify top and bottom particle layers by y_cm position.

    Parameters
    ----------
    particles  : list of CapsuleParticle
    Ly         : float — box height (y extent, 0 to Ly)
    n_layers   : int   — not used; kept for API compatibility
    layer_frac : float — fraction of box height to define each layer (default 0.20)

    Returns
    -------
    top_idx    : list of int — indices of top-layer particles
    bot_idx    : list of int — indices of bottom-layer particles
    y_thresh_top, y_thresh_bot : float — threshold y positions used
    """
    y_cms = np.array([p.x_cm[1] for p in particles])
    y_mid = np.mean(y_cms)

    # Compute box extent from particle positions (may not start at 0)
    y_min = y_cms.min() - (1.0 + 2.0 * np.pi / particles[0].N)
    y_max = y_cms.max() + (1.0 + 2.0 * np.pi / particles[0].N)
    height = y_max - y_min

    y_thresh_top = y_max - layer_frac * height
    y_thresh_bot = y_min + layer_frac * height

    top_idx = [i for i, p in enumerate(particles) if p.x_cm[1] >= y_thresh_top]
    bot_idx = [i for i, p in enumerate(particles) if p.x_cm[1] <= y_thresh_bot]

    return top_idx, bot_idx, y_thresh_top, y_thresh_bot
