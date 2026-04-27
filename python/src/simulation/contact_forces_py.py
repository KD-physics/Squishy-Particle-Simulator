"""
contact_forces_py.py — Python inter-capsule force kernel using CapCandidates.

ONE-DIRECTIONAL design (matches frozen branch to machine precision):
  - CapCandidates is asymmetric: row eA on particle A lists candidates on B only
    when pA < pB.  Row eB on B does NOT list eA (no double-counting).
  - Forces computed from SOURCE Gauss points only (matches frozen _capsule_capsule_forces).
  - Newton 3rd applied to candidate nodes (scatter in Python, segment_sum in TF Phase 3.3).

Why this matches frozen:
  The frozen _capsule_capsule_forces uses Gauss points on A to find closest points on B
  and applies +F to A, -F to B.  This kernel does the same per row, giving bit-identical
  forces when candidates match the frozen brute-force pairs.
"""

import numpy as np

_GAUSS2_S = np.array([(1.0 - 1.0 / np.sqrt(3)) / 2.0,
                       (1.0 + 1.0 / np.sqrt(3)) / 2.0])
_GAUSS2_W = np.array([0.5, 0.5])


def inter_capsule_forces_py(x_all, CapCandidates, r_c, k_c, L0):
    """
    Compute inter-capsule contact forces from CapCandidates index matrix.

    Parameters
    ----------
    x_all         : (P, N, 2) float64
    CapCandidates : (K, E)    int32   — 0 = ghost; one-directional (pA < pB rows only)
    r_c, k_c, L0  : float — contact radius, spring constant, edge length (uniform)

    Returns
    -------
    f_contact : (P, N, 2) float64
    """
    P, N, _ = x_all.shape
    K, E    = CapCandidates.shape

    # Build flat ghost-prepended edge endpoint arrays (index 0 = ghost)
    ghost   = np.array([[1e9, 1e9]], dtype=np.float64)
    x_flat      = x_all.reshape(P * N, 2)
    x_flat_next = np.roll(x_all, -1, axis=1).reshape(P * N, 2)
    x_start = np.vstack([ghost, x_flat])        # (K+1, 2); index 0=ghost, 1..K = edges
    x_end   = np.vstack([ghost, x_flat_next])

    f_contact = np.zeros((P * N, 2), dtype=np.float64)

    contact_r = 2.0 * r_c
    s = _GAUSS2_S
    w = _GAUSS2_W

    for k in range(K):
        row = CapCandidates[k]
        active = row[row != 0]
        if len(active) == 0:
            continue

        p_src = k // N
        e_src = k  % N
        a0 = x_all[p_src, e_src]
        a1 = x_all[p_src, (e_src + 1) % N]

        # 2-pt Gauss on source edge (same as frozen branch)
        for g_idx in range(2):
            sg  = s[g_idx]
            wg  = w[g_idx]
            xq  = (1.0 - sg) * a0 + sg * a1          # Gauss point on A

            # Gather candidate edge endpoints
            b0_arr = x_start[active]                  # (n_act, 2)
            b1_arr = x_end[active]                    # (n_act, 2)
            cand_0idx = active - 1                    # 0-indexed into x_flat

            # Closest point on each candidate edge to xq
            ab   = b1_arr - b0_arr                    # (n_act, 2)
            ab2  = np.einsum('ij,ij->i', ab, ab)      # (n_act,)
            t_B  = np.clip(
                np.where(ab2 > 1e-30,
                         np.dot(b0_arr - xq, -ab.T).diagonal() / np.maximum(ab2, 1e-30),
                         0.5 * np.ones(len(active))),
                0.0, 1.0)
            # Correct: t_B[i] = dot(xq - b0[i], ab[i]) / ab2[i]
            t_B = np.array([
                float(np.clip(np.dot(xq - b0_arr[i], ab[i]) / max(ab2[i], 1e-30), 0, 1))
                for i in range(len(active))
            ])
            cp   = b0_arr + t_B[:, None] * ab         # (n_act, 2)
            diff = xq[None, :] - cp                    # push source away: xq - cp
            d    = np.linalg.norm(diff, axis=1)        # (n_act,)
            gap  = d - contact_r                       # (n_act,)

            act2 = gap < 0.0
            if not np.any(act2):
                continue

            d_act  = np.where(d[act2] > 1e-15, d[act2], 1.0)
            n_hat  = diff[act2] / d_act[:, None]       # (k2, 2)
            F_mag  = k_c * (-gap[act2]) * L0 * wg     # (k2,)
            F_vec  = F_mag[:, None] * n_hat             # (k2, 2) — pushes A away

            # ── Source (A) side ──────────────────────────────────────────────
            # No scatter: accumulate to source edge endpoints
            node_a0 = p_src * N + e_src
            node_a1 = p_src * N + (e_src + 1) % N
            np.add.at(f_contact, node_a0, (1.0 - sg) * F_vec.sum(axis=0))
            np.add.at(f_contact, node_a1, sg          * F_vec.sum(axis=0))

            # ── Candidate (B) side — Newton 3rd ─────────────────────────────
            # Scatter to candidate edge endpoints (unavoidable for 1-sided matching)
            active_act   = cand_0idx[act2]             # (k2,) 0-indexed
            t_B_act      = t_B[act2]
            p_b  = active_act // N
            e_b  = active_act  % N
            for idx in range(len(active_act)):
                node_b0 = p_b[idx] * N + e_b[idx]
                node_b1 = p_b[idx] * N + (e_b[idx] + 1) % N
                np.add.at(f_contact, node_b0, -(1.0 - t_B_act[idx]) * F_vec[idx])
                np.add.at(f_contact, node_b1, -t_B_act[idx]          * F_vec[idx])

    return f_contact.reshape(P, N, 2)
