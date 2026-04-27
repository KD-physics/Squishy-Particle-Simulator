"""
two_disk.py — Two-disk Hertzian contact simulation using FEM + analytic Hertz forces.

For Phase 1-2: simulates a single disk with Hertz contact traction applied as
Neumann BC. The contact force is computed from analytic Hertz theory.
For Phase 2: validates that the FEM deformation matches Hertz predictions.

Workflow:
    1. Prescribe contact force F between two disks (Hertz input)
    2. Compute Hertz contact geometry (a, P_max) from F and material params
    3. Apply Hertz traction profile as Neumann BC to each disk
    4. Solve FEM for each disk
    5. Extract perimeter displacements
    6. Verify against Hertz approach displacement formula
"""

import numpy as np
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

from src.simulation.disk_mesh import make_disk_mesh
from src.simulation.fem_elastic import solve_disk_neumann, analytic_uniform_pressure_disp
from src.contact.hertz import hertz_contact_2d, hertz_pressure_profile, hertz_traction_on_disk
from src.contact.extract import extract_forces, perimeter_forces


def make_hertz_traction(
    contact_normal: np.ndarray,  # (2,) unit normal pointing INTO disk 1
    contact_center: np.ndarray,  # (2,) contact point on disk 1 boundary
    a: float,
    F: float,
):
    """
    Build a traction_fn for Hertz contact traction on a disk.

    The traction is applied only near the contact point, following the
    elliptic Hertz pressure profile along the tangential direction.

    The function signature matches what solve_disk_neumann expects:
        traction_fn(x_q, y_q, nx_q, ny_q) -> (tx, ty)
    where inputs are quadrature point arrays.
    """
    tangent = np.array([-contact_normal[1], contact_normal[0]])
    P_max = 2.0 * F / (np.pi * a)
    # Cutoff radius around contact center: traction is zero beyond 2a from the center
    cutoff = 2.0 * a

    def traction_fn(x_q, y_q, nx_q, ny_q):
        # Relative position from contact center
        rel_x = x_q - contact_center[0]
        rel_y = y_q - contact_center[1]

        # Distance from contact center (must be small for contact zone)
        dist = np.sqrt(rel_x**2 + rel_y**2)

        # Tangential coordinate along contact line
        s = rel_x * tangent[0] + rel_y * tangent[1]

        # Hertz pressure: non-zero only within contact half-width AND near contact center
        in_zone = (np.abs(s) <= a) & (dist <= cutoff)
        P = np.where(in_zone,
                     P_max * np.sqrt(np.maximum(0.0, 1.0 - (s / a)**2)),
                     0.0)

        # Traction direction: inward normal (contact_normal points INTO disk)
        tx = P * contact_normal[0]
        ty = P * contact_normal[1]
        return tx, ty

    return traction_fn


def run_two_disk_hertz(
    R: float = 1.0,
    E: float = 1e5,
    nu: float = 0.3,
    F: float = 100.0,
    N_perimeter: int = 32,
    save_dir: str = "results",
) -> dict:
    """
    Simulate two equal elastic disks in Hertzian contact.

    Disk 1 is centered at (-R, 0); disk 2 at (+R, 0).
    Contact at origin. Normal direction: x-axis.

    Parameters
    ----------
    R : disk radius (m)
    E : Young's modulus (Pa)
    nu : Poisson's ratio
    F : normal contact force per unit depth (N/m)
    N_perimeter : number of perimeter nodes
    save_dir : directory for output plots

    Returns
    -------
    results : dict with simulation and analytic comparison data
    """
    # ── Hertz theory ──────────────────────────────────────────────────────
    hertz = hertz_contact_2d(F, R, R, E, nu)
    a = hertz["a"]
    delta_theory = hertz["delta"]

    # ── Mesh for each disk (identical by symmetry) ─────────────────────────
    mesh = make_disk_mesh(R=R, N_perimeter=N_perimeter)

    # ── Disk 1: at (-R, 0), contact on right (+x direction) ───────────────
    # Contact on disk 1: contact normal pointing in +x (from disk 1 toward disk 2)
    n1 = np.array([+1.0, 0.0])   # outward normal toward contact
    c1 = np.array([+R, 0.0])     # contact point on disk 1 perimeter

    traction_1 = make_hertz_traction(
        contact_normal=-n1,   # inward: force on disk 1 points in -x
        contact_center=c1,
        a=a, F=F,
    )

    u1, f_int1, meta1 = solve_disk_neumann(mesh, traction_1, E, nu)

    # ── Disk 2: by symmetry same as disk 1, mirrored ──────────────────────
    # (Contact normal on disk 2 points in -x direction, contact at -R)
    n2 = np.array([-1.0, 0.0])
    c2 = np.array([-R, 0.0])

    traction_2 = make_hertz_traction(
        contact_normal=-n2,   # inward: force on disk 2 points in +x
        contact_center=c2,
        a=a, F=F,
    )

    u2, f_int2, meta2 = solve_disk_neumann(mesh, traction_2, E, nu)

    # ── Approach displacement from simulation ──────────────────────────────
    # The approach delta = (displacement of disk 1 contact point) - (disk 2 contact point)
    # Disk 1's contact node: the perimeter node closest to (+R, 0)
    peri_pos1 = mesh["perimeter_pos"]
    peri_ids1 = mesh["perimeter_ids"]
    idx_c1 = np.argmin(np.linalg.norm(peri_pos1 - c1, axis=1))
    u_contact_1_x = u1[peri_ids1[idx_c1], 0]  # x-displacement at contact node of disk 1

    # By symmetry disk 2 has same magnitude, opposite sign at its contact node
    idx_c2 = np.argmin(np.linalg.norm(peri_pos1 - c2, axis=1))
    u_contact_2_x = u2[peri_ids1[idx_c2], 0]  # u_x at (-R, 0) for disk 2

    # Approach = compression of gap = u1_x(contact) - u2_x(contact)
    # Disk 1 at -R: contact node approaches in +x direction (positive ux)
    # Wait, disk 1 is at origin (we use reference mesh), contact on the RIGHT
    # In reference frame: contact node of disk 1 at (+R,0) moves in -x (compression)
    # contact node of disk 2 at (-R,0) moves in +x (compression)
    # Total approach = |u1_contact_x| + |u2_contact_x|
    delta_sim = abs(u_contact_1_x) + abs(u_contact_2_x)

    # ── Contact half-width from simulation ────────────────────────────────
    # Measure contact half-width from the non-zero traction zone at perimeter
    # Use the perimeter node displacements as proxy for contact zone extent
    forces1 = extract_forces(mesh, u1, traction_1, E, nu)
    pf1 = perimeter_forces(forces1)

    Fext_mag1 = np.linalg.norm(pf1["Fext"], axis=1)
    # Contact zone: nodes where |Fext| > threshold (1% of max)
    thresh = 0.01 * np.max(Fext_mag1)
    contact_nodes_1 = pf1["X_ref"][Fext_mag1 > thresh]
    if len(contact_nodes_1) >= 2:
        # Project onto tangent direction (y-axis for x-aligned contact)
        y_contact = contact_nodes_1[:, 1]
        a_sim = (np.max(y_contact) - np.min(y_contact)) / 2.0
    else:
        a_sim = 0.0

    # ── Summary ───────────────────────────────────────────────────────────
    if a > 0:
        a_err = abs(a_sim - a) / a
        delta_err = abs(delta_sim - abs(delta_theory)) / abs(delta_theory)
    else:
        a_err = delta_err = np.nan

    results = {
        "R": R, "E": E, "nu": nu, "F": F,
        "N_perimeter": N_perimeter,
        # Hertz theory
        "a_theory": a,
        "delta_theory": abs(delta_theory),
        "P_max_theory": hertz["P_max"],
        # Simulation
        "a_sim": a_sim,
        "delta_sim": delta_sim,
        # Errors
        "a_err": a_err,
        "delta_err": delta_err,
        # Solver metadata
        "meta1": meta1,
        "forces1": forces1,
        "u1": u1,
        "mesh": mesh,
    }

    return results


def visualize_two_disk(results: dict, save_dir: str = "results") -> None:
    """Plot the deformed disk and contact force distribution."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    os.makedirs(save_dir, exist_ok=True)

    mesh = results["mesh"]
    u1 = results["u1"]
    v = mesh["vertices"]
    t = mesh["triangles"]

    # Scale displacement for visualization
    scale = 5.0
    v_def = v + scale * u1

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: deformed disk
    ax = axes[0]
    tri_orig = mtri.Triangulation(v[:, 0], v[:, 1], t)
    tri_def = mtri.Triangulation(v_def[:, 0], v_def[:, 1], t)
    ax.triplot(tri_orig, 'b-', lw=0.3, alpha=0.3, label='reference')
    ax.triplot(tri_def, 'r-', lw=0.5, label=f'deformed (×{scale})')
    ax.set_aspect('equal')
    ax.set_title(f'Deformed Disk 1 (R={results["R"]}, E={results["E"]:.1e}, F={results["F"]})')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.2)

    # Right: contact force distribution at perimeter
    ax = axes[1]
    forces1 = results["forces1"]
    pf1 = perimeter_forces(forces1)
    peri_pos = pf1["X_ref"]
    fext_mag = np.linalg.norm(pf1["Fext"], axis=1)
    angles = np.arctan2(peri_pos[:, 1], peri_pos[:, 0])
    sort_idx = np.argsort(angles)
    ax.bar(np.degrees(angles[sort_idx]),
           fext_mag[sort_idx],
           width=np.degrees(2*np.pi/len(angles)), align='center')
    ax.set_xlabel('Perimeter angle (degrees)')
    ax.set_ylabel('|Fext| (N/m)')
    ax.set_title(f'Contact Force at Perimeter\na_theory={results["a_theory"]:.4f}m')
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f'Two-Disk Hertz Contact: a_err={results["a_err"]:.3%}, '
        f'δ_err={results["delta_err"]:.3%}',
        fontsize=12)
    plt.tight_layout()
    fname = os.path.join(save_dir, f'two_disk_N{results["N_perimeter"]}.png')
    fig.savefig(fname, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {fname}")


if __name__ == "__main__":
    print("Two-Disk Hertz Contact Simulation")
    print("=" * 50)

    R, E, nu, F = 1.0, 1e5, 0.3, 100.0

    for N in [16, 32]:
        print(f"\nN_perimeter = {N}:")
        res = run_two_disk_hertz(R=R, E=E, nu=nu, F=F, N_perimeter=N)
        print(f"  a_theory = {res['a_theory']:.4e} m")
        print(f"  a_sim    = {res['a_sim']:.4e} m  err={res['a_err']:.3%}")
        print(f"  delta_theory = {res['delta_theory']:.4e} m")
        print(f"  delta_sim    = {res['delta_sim']:.4e} m  err={res['delta_err']:.3%}")
        visualize_two_disk(res)
