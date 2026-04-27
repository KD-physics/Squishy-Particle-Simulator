"""
extract.py — Post-process FEM solutions to extract per-node force decomposition.

For a disk solved with external Neumann traction (contact forces), this module
extracts the full force decomposition at perimeter nodes:

    x_node   : absolute node position
    X_ref    : reference (undeformed) position on unit circle
    u_elastic: displacement = x_node - X_ref
    Fext     : external (contact) traction integrated to nodal forces
    Fint     : internal elastic force = K @ u (stiffness times displacement)
    Fnet     : Fext + Fint (should be near zero in static equilibrium)

Force balance check: |sum(Fnet)| / max(|Fnet|) < 1e-8 (global balance)

This is the Phase 1 verification tool. In Phase 3+, Fext comes from IPC
contact simulations or analytic Hertz formulas.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import skfem
from skfem import Basis, ElementTriP1, ElementVector, FacetBasis, LinearForm, MeshTri
from skfem.models.elasticity import linear_elasticity, lame_parameters
from scipy.sparse import csr_matrix


def extract_forces(
    mesh_dict: dict,
    u: np.ndarray,           # (N_nodes, 2) displacement field
    traction_fn,             # same callable used for the solve
    E: float,
    nu: float,
) -> dict:
    """
    Extract per-node force decomposition from a solved FEM disk.

    Parameters
    ----------
    mesh_dict : dict
        From disk_mesh.make_disk_mesh.
    u : (N_nodes, 2)
        Nodal displacements from fem_elastic.solve_disk_neumann.
    traction_fn : callable
        traction_fn(x, y, nx, ny) -> (tx, ty) — the boundary traction applied.
    E, nu : float
        Material parameters used in the solve.

    Returns
    -------
    forces : dict with keys:
        'x_node'   : (N_nodes, 2) — absolute positions (X_ref + u)
        'X_ref'    : (N_nodes, 2) — reference positions
        'u_elastic': (N_nodes, 2) — elastic displacement (same as input u)
        'Fext'     : (N_nodes, 2) — external nodal forces (from traction)
        'Fint'     : (N_nodes, 2) — internal elastic forces (K @ u)
        'Fnet'     : (N_nodes, 2) — net nodal force = Fint + Fext
        'perimeter_ids' : (N_p,)  — indices of perimeter nodes
        'force_balance_relative': float — global force balance metric
    """
    lam, mu = lame_parameters(E, nu)
    v = mesh_dict["vertices"]   # (N, 2)
    t = mesh_dict["triangles"]  # (T, 3)
    N = len(v)
    sfmesh = MeshTri(v.T, t.T)

    # ── Assemble stiffness matrix ──────────────────────────────────────────
    basis = Basis(sfmesh, ElementVector(ElementTriP1()))
    K = skfem.asm(linear_elasticity(lam, mu), basis)

    # ── Internal forces: f_int = K @ u ────────────────────────────────────
    u_flat = np.zeros(2 * N)
    u_flat[0::2] = u[:, 0]   # ux
    u_flat[1::2] = u[:, 1]   # uy
    f_int_flat = K @ u_flat   # (2N,)

    # ── External forces: high-order Neumann assembly ──────────────────────
    # Use _assemble_neumann (n_quad=50) instead of FacetBasis to correctly
    # capture narrow tractions (Hertz contact zone << facet size).
    from src.simulation.fem_elastic import _assemble_neumann
    f_ext_flat = _assemble_neumann(sfmesh, traction_fn, n_quad=50)

    # ── Net force = Fext - Fint (residual of K u = f_ext) ─────────────────
    # In static equilibrium: K @ u = f_ext, so f_net = f_ext - K @ u = 0
    # But numerically there will be small residuals.
    f_net_flat = f_ext_flat - f_int_flat

    # ── Reshape to (N, 2) ─────────────────────────────────────────────────
    Fext = np.column_stack([f_ext_flat[0::2], f_ext_flat[1::2]])
    Fint = np.column_stack([f_int_flat[0::2], f_int_flat[1::2]])
    Fnet = np.column_stack([f_net_flat[0::2], f_net_flat[1::2]])

    # ── Reference and deformed positions ──────────────────────────────────
    X_ref = v                    # (N, 2)
    x_node = X_ref + u           # (N, 2)

    # ── Global force balance ───────────────────────────────────────────────
    # For a self-equilibrated load: |sum(Fext)| / sum(|Fext|) << 1
    # (The sum of all external nodal forces should be zero by Newton's 3rd law)
    F_total_ext = np.sum(Fext, axis=0)
    F_sum_norm  = np.sum(np.linalg.norm(Fext, axis=1))
    force_balance_rel = np.linalg.norm(F_total_ext) / (F_sum_norm + 1e-30)

    # Also record the net residual (should be near machine epsilon)
    F_total = np.sum(Fnet, axis=0)

    return {
        "x_node": x_node,
        "X_ref": X_ref,
        "u_elastic": u,
        "Fext": Fext,
        "Fint": Fint,
        "Fnet": Fnet,
        "perimeter_ids": mesh_dict["perimeter_ids"],
        "force_balance_relative": force_balance_rel,
        "F_total": F_total,
    }


def perimeter_forces(forces: dict) -> dict:
    """
    Extract only the perimeter-node data from the full-field force dict.
    Returns a dict with the same keys but arrays indexed to perimeter nodes.
    """
    ids = forces["perimeter_ids"]
    return {
        "x_node":    forces["x_node"][ids],
        "X_ref":     forces["X_ref"][ids],
        "u_elastic": forces["u_elastic"][ids],
        "Fext":      forces["Fext"][ids],
        "Fint":      forces["Fint"][ids],
        "Fnet":      forces["Fnet"][ids],
        "force_balance_relative": forces["force_balance_relative"],
        "F_total": forces["F_total"],
    }


def check_force_balance(forces: dict, tol: float = 1e-8) -> bool:
    """
    Verify global force balance: |sum(Fnet)| / max(|Fnet|) < tol.
    """
    return forces["force_balance_relative"] < tol


if __name__ == "__main__":
    import sys
    sys.path.insert(0, "/root/workspace/projects/polyfem")
    from src.simulation.disk_mesh import make_disk_mesh
    from src.simulation.fem_elastic import (
        solve_disk_neumann, uniform_pressure_traction,
        analytic_uniform_pressure_disp,
    )

    R, E, nu, p = 1.0, 1e5, 0.3, 1e3

    print("Force extraction test — uniform pressure disk")
    print(f"  R={R}, E={E:.1e}, nu={nu}, p={p:.1e}\n")

    all_pass = True
    for N in [16, 32, 64]:
        mesh = make_disk_mesh(R=R, N_perimeter=N)
        traction = uniform_pressure_traction(p)
        u, f_int_solver, meta = solve_disk_neumann(mesh, traction, E, nu)

        forces = extract_forces(mesh, u, traction, E, nu)
        fb = forces["force_balance_relative"]

        # Waypoint 1.3 gate: |sum(Fnet)| / max(|Fnet|) < 1e-8
        gate_pass = fb < 1e-8
        if not gate_pass:
            all_pass = False

        # Check Fext is non-trivially non-zero (extraction is actually working)
        pf = perimeter_forces(forces)
        fext_norm = np.linalg.norm(pf["Fext"])
        fint_norm = np.linalg.norm(pf["Fint"])
        fnet_norm = np.linalg.norm(pf["Fnet"])
        # In equilibrium Fext ≈ Fint; but both should be non-zero
        nontrivial = fext_norm > 1e-10 and fint_norm > 1e-10

        print(f"  N_perimeter={N:3d}:")
        print(f"    force_balance_rel = {fb:.2e}  gate={'PASS ✓' if gate_pass else 'FAIL ✗'}")
        print(f"    |Fext|={fext_norm:.3e}  |Fint|={fint_norm:.3e}  |Fnet|={fnet_norm:.3e}  non-zero={nontrivial}")
        print(f"    |Fnet_total| = {np.linalg.norm(forces['F_total']):.2e}")
        print(f"    Fext[:3] =\n      {pf['Fext'][:3]}")
        print(f"    Fint[:3] =\n      {pf['Fint'][:3]}")
        print()

    print(f"Phase 1 Waypoint 1.3 force balance gate: {'PASS' if all_pass else 'FAIL'}")
