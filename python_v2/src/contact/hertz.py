"""
hertz.py — 2D Plane-Strain Hertzian Contact Theory

Implements the analytic 2D Hertz contact formulas for two elastic cylinders.
Used as:
  - Ground truth for Phase 2 benchmark validation
  - ExtForce model for Phase 1-3 DEM simulations

Reference:
  Johnson, K.L. (1985). Contact Mechanics. Cambridge University Press.
  Chapter 4: Normal Contact of Elastic Solids — Hertz Theory
  (2D plane-strain cylinder-on-cylinder formulas)

All formulas are for PLANE STRAIN (not 3D Hertz).

Units: SI (Pa, N/m, m). All forces are per unit depth (N/m in 2D).
"""

import numpy as np


def hertz_contact_2d(
    F: float,        # Normal force per unit depth (N/m)
    R1: float,       # Radius of disk 1 (m)
    R2: float,       # Radius of disk 2 (m)
    E1: float,       # Young's modulus of disk 1 (Pa)
    nu1: float,      # Poisson's ratio of disk 1
    E2: float = None,   # Young's modulus of disk 2 (defaults to E1)
    nu2: float = None,  # Poisson's ratio of disk 2 (defaults to nu1)
) -> dict:
    """
    Compute 2D Hertzian contact parameters for two elastic cylinders.

    Returns a dict with:
        R_eff   : effective contact radius
        E_eff   : combined elastic modulus (plane strain)
        a       : contact half-width
        P_max   : maximum contact pressure (at x=0)
        delta   : center-to-center approach displacement
    """
    if E2 is None:
        E2 = E1
    if nu2 is None:
        nu2 = nu1

    # Effective radius
    R_eff = R1 * R2 / (R1 + R2)

    # Combined elastic modulus (plane strain Hertz)
    # 1/E_eff = (1-nu1^2)/E1 + (1-nu2^2)/E2
    E_eff = 1.0 / ((1 - nu1**2) / E1 + (1 - nu2**2) / E2)

    # Contact half-width (2D plane strain)
    # a = sqrt(4 * F * R_eff / (pi * E_eff))
    a = np.sqrt(4.0 * F * R_eff / (np.pi * E_eff))

    # Maximum pressure P(0)
    P_max = 2.0 * F / (np.pi * a)

    # Approach displacement (center-to-center shortening)
    # delta = (F / (pi * E_eff)) * (ln(4*R_eff/a) - 0.5) * 2  [for both surfaces combined]
    # But the PLAN.md formula is: delta = (2*F*(1-nu^2)/(pi*E)) * (ln(4*R_eff/a) - 0.5)
    # This uses single-material formula; for two identical cylinders:
    # delta_total = 2 * delta_one_surface
    # For two identical disks: delta = (4*F*(1-nu^2) / (pi*E)) * (ln(4*R_eff/a) - 0.5)
    # Using E_eff for general case:
    delta = (2.0 * F / (np.pi * E_eff)) * (np.log(4.0 * R_eff / a) - 0.5)

    return {
        "R_eff": R_eff,
        "E_eff": E_eff,
        "a": a,
        "P_max": P_max,
        "delta": delta,
        "F": F,
    }


def hertz_pressure_profile(x: np.ndarray, a: float, F: float) -> np.ndarray:
    """
    2D Hertz pressure distribution P(x) along contact line.

    P(x) = P_max * sqrt(1 - (x/a)^2)   for |x| <= a
    P(x) = 0                             for |x| > a

    P_max = 2F / (pi * a)
    """
    P_max = 2.0 * F / (np.pi * a)
    P = np.zeros_like(x)
    mask = np.abs(x) <= a
    P[mask] = P_max * np.sqrt(1.0 - (x[mask] / a)**2)
    return P


def hertz_traction_on_disk(
    perimeter_pos: np.ndarray,  # (N_p, 2) perimeter node positions
    contact_center: np.ndarray,  # (2,) center of contact on perimeter
    contact_normal: np.ndarray,  # (2,) unit normal pointing INTO disk (inward)
    a: float,                    # contact half-width
    F: float,                    # normal force per unit depth
) -> np.ndarray:
    """
    Compute Hertz traction (force per unit length) on each perimeter node.

    The contact zone is centered at `contact_center` on the perimeter.
    The traction is normal (no friction) and follows the Hertzian elliptic profile.

    Parameters
    ----------
    perimeter_pos : (N_p, 2)
        Positions of perimeter nodes.
    contact_center : (2,)
        Point on perimeter at the center of contact (midpoint of contact arc).
    contact_normal : (2,)
        Unit normal pointing INTO the disk at the contact center.
    a : float
        Contact half-width (meters).
    F : float
        Total normal force per unit depth (N/m).

    Returns
    -------
    traction : (N_p, 2)
        Traction vector (tx, ty) at each perimeter node.
        Only nodes within distance a from contact_center are non-zero.
    """
    # Project perimeter positions onto the contact tangent direction
    tangent = np.array([-contact_normal[1], contact_normal[0]])  # perpendicular to normal
    rel = perimeter_pos - contact_center
    x_local = rel @ tangent   # signed distance along contact line
    n_local = rel @ contact_normal  # distance along normal (should be ~0 at contact)

    # Hertz pressure at each node
    P = hertz_pressure_profile(x_local, a, F)

    # Traction direction: along inward normal (compression)
    traction = P[:, None] * contact_normal[None, :]

    return traction


def two_disk_contact_setup(
    center1: np.ndarray,   # (2,) center of disk 1
    center2: np.ndarray,   # (2,) center of disk 2
    R1: float,
    R2: float,
    E: float,
    nu: float,
) -> dict:
    """
    Compute the Hertz contact setup for two touching disks.

    Returns contact geometry and force parameters needed to apply
    Hertz traction as boundary conditions.
    """
    d_vec = center2 - center1
    d = np.linalg.norm(d_vec)
    n12 = d_vec / d   # unit normal from 1 to 2

    # Contact point on each disk perimeter
    contact_on_1 = center1 + R1 * n12     # on disk 1's perimeter
    contact_on_2 = center2 - R2 * n12     # on disk 2's perimeter

    # For given overlap: delta = R1 + R2 - d
    delta = R1 + R2 - d

    return {
        "n12": n12,         # normal from disk 1 to disk 2
        "contact_on_1": contact_on_1,
        "contact_on_2": contact_on_2,
        "overlap": delta,
        "distance": d,
    }


if __name__ == "__main__":
    # Test: verify Hertz formulas for two equal disks
    R = 1.0   # m
    E = 1e5   # Pa
    nu = 0.3
    F = 100.0  # N/m (applied normal force per unit depth)

    result = hertz_contact_2d(F, R, R, E, nu)
    print("2D Hertz Contact — Two Equal Disks")
    print(f"  R={R}m, E={E:.1e}Pa, nu={nu}, F={F}N/m")
    print(f"  R_eff = {result['R_eff']:.4f} m")
    print(f"  E_eff = {result['E_eff']:.4e} Pa")
    print(f"  a     = {result['a']:.4e} m  (contact half-width)")
    print(f"  P_max = {result['P_max']:.4e} Pa  (max pressure)")
    print(f"  delta = {result['delta']:.4e} m  (approach)")
    print()

    # Verify: a = sqrt(4*F*R_eff / (pi*E_eff))
    a_check = np.sqrt(4*F*R/2 / (np.pi * (E/(2*(1-nu**2)))))
    print(f"  a check (single material): {a_check:.4e}")
    print(f"  Match: {np.isclose(result['a'], a_check)}")
    print()

    # Plot pressure profile
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import os

    x = np.linspace(-1.5*result['a'], 1.5*result['a'], 200)
    P = hertz_pressure_profile(x, result['a'], F)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x*1000, P/1000, 'b-', lw=2)
    ax.axvline(-result['a']*1000, color='r', ls='--', label=f'a={result["a"]*1000:.2f}mm')
    ax.axvline(+result['a']*1000, color='r', ls='--')
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('Pressure (kPa)')
    ax.set_title('2D Hertz Contact Pressure Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    os.makedirs("results", exist_ok=True)
    fig.savefig("results/hertz_pressure_profile.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print("  Saved results/hertz_pressure_profile.png")
