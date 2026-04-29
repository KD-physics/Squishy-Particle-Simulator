"""
candidacy_manager.py — Python candidacy manager for TF sim pipeline.

Maintains CapCandidates (K, E) int32 array where:
  K = P * N  — total edges (edge k = particle p * N + edge e, k in [1..K])
  E          — max candidates per edge (skin-buffered, set at init)
  0          — ghost index (far-away, always zero force)

Three-level filtering:
  Level 1: center-center distance (cell list, O(P))
  Level 2: contact-normal registration — find facing edge pair via angle arithmetic
  Level 3: index-sliding fill ±dj around the registered contact edge

C++ will replace this class in Phase 3.4. Interface is identical.
"""

import numpy as np
import json
import os


class CandidacyManager:
    """
    Python implementation of the capsule candidacy manager.

    Parameters
    ----------
    P        : int   — number of particles
    N        : int   — nodes per particle (uniform)
    R0       : float — particle radius (uniform)
    E        : int   — max candidates per edge (must be >= actual max contacts)
    skin     : float — skin distance beyond contact threshold for candidate inclusion
    dj       : int   — half-width of index window around registered contact edge
    periodic : bool  — if True, use periodic BC (minimum-image) in BOTH x and y
    Lx, Ly   : float — box dimensions (required when periodic=True)
    log_dir  : str or None — if set, log every update call for Phase 3.4 corpus
    """

    def __init__(self, P, N, R0, E=32, skin=0.3, dj=None,
                 periodic=False, periodic_x=None, periodic_y=None,
                 Lx=None, Ly=None, log_dir=None,
                 R0_arr=None, skin_arr=None):
        # dj default: covers ~45° half-arc per contact; non-overlapping for
        # typical multi-neighbor configurations (equilateral, hex packing).
        # Formula: dj = max(3, N//8) → window = 2*dj+1 edges ≈ N/4 arc.
        if dj is None:
            dj = max(3, N // 8)
        # Per-axis periodicity: periodic_x/y override the global `periodic` flag.
        # Use periodic_x=True, periodic_y=False for walls=[0,1] (x-periodic, hard y).
        if periodic_x is None:
            periodic_x = periodic
        if periodic_y is None:
            periodic_y = periodic
        self.P          = P
        self.N          = N
        self.E          = E
        self.skin       = skin
        self.dj         = dj
        self.periodic   = periodic
        self.periodic_x = periodic_x
        self.periodic_y = periodic_y
        self.Lx         = float(Lx) if Lx is not None else 1e9
        self.Ly         = float(Ly) if Ly is not None else 1e9
        self.log_dir    = log_dir
        if (periodic_x or periodic_y) and (Lx is None or Ly is None):
            raise ValueError("periodic_x/y=True requires Lx and Ly")

        # Patch A: per-particle R0 (polydisperse / mixed-type support)
        if R0_arr is None:
            self._R0_arr   = np.full(P, float(R0))
        else:
            self._R0_arr   = np.asarray(R0_arr, dtype=np.float64)
        self.R0 = float(np.mean(self._R0_arr))  # display / needs_update heuristic

        # Patch D: per-particle skin (mixed-type support, defaults to scalar skin)
        if skin_arr is None:
            self._skin_arr = np.full(P, float(skin))
        else:
            self._skin_arr = np.asarray(skin_arr, dtype=np.float64)

        # CapCandidates: (K, E) int32, initialized to ghost (0)
        K                   = P * N
        self.K              = K
        self.CapCandidates  = np.zeros((K, E), dtype=np.int32)

        # Update-trigger tracking
        self._last_x_cm  = None   # (P, 2) float64
        self._last_theta = None   # (P,)  float64
        self._n_updates  = 0

        # Debug / corpus logging
        self._corpus     = []     # list of (x_cm, theta, CapCandidates) snapshots
        self._debug      = False  # enable brute-force missing-contact check

        # Stats
        self.utilization_history = []

        if log_dir is not None:
            os.makedirs(log_dir, exist_ok=True)

    # ── public interface ──────────────────────────────────────────────────────

    def needs_update(self, x_cm, theta):
        """
        Return True if CapCandidates should be refreshed.
        Trigger: max(|Δx_cm_mi| + R0 * |Δtheta|) > skin/2 for any particle.
        Uses minimum-image displacement when periodic=True.
        On first call always returns True.
        """
        if self._last_x_cm is None:
            return True
        delta = x_cm - self._last_x_cm          # (P, 2)
        if self.periodic:
            delta[:, 0] -= self.Lx * np.round(delta[:, 0] / self.Lx)
            delta[:, 1] -= self.Ly * np.round(delta[:, 1] / self.Ly)
        dx      = np.max(np.linalg.norm(delta, axis=1))
        dtheta  = np.max(np.abs(theta - self._last_theta))
        return (dx + self.R0 * dtheta) > self.skin * 0.5

    def update(self, x_cm, theta):
        """
        Recompute CapCandidates from current particle centers and orientations.
        Writes CapCandidates in place.

        Parameters
        ----------
        x_cm  : (P, 2) float64 — particle centers
        theta : (P,)   float64 — particle orientations
        """
        x_cm  = np.asarray(x_cm,  dtype=np.float64)
        theta = np.asarray(theta, dtype=np.float64)

        self.CapCandidates[:] = 0   # reset to ghost

        neighbor_pairs = self._level1_pairs(x_cm)
        for (pA, pB) in neighbor_pairs:
            self._fill_pair(pA, pB, x_cm, theta)

        self._last_x_cm  = x_cm.copy()
        self._last_theta = theta.copy()
        self._n_updates += 1

        # Log corpus entry
        if self.log_dir is not None:
            self._corpus.append({
                'x_cm':           x_cm.tolist(),
                'theta':          theta.tolist(),
                'CapCandidates':  self.CapCandidates.tolist(),
            })

        # Utilization stats
        active = np.any(self.CapCandidates != 0, axis=1)
        used   = np.sum(self.CapCandidates != 0, axis=1)
        self.utilization_history.append(float(np.mean(used[active])) if active.any() else 0.0)

        # Warn if any row is near capacity
        max_used = np.max(np.sum(self.CapCandidates != 0, axis=1))
        if max_used >= self.E - 2:
            print(f"WARNING: CandidacyManager E={self.E} near capacity "
                  f"(max_used={max_used}). Increase E.")

    def step(self, x_cm, theta):
        """Alias for update() matching the C++ pybind11 interface."""
        self.update(x_cm, theta)

    # ── level 1: center-center cell list ─────────────────────────────────────

    def _min_image(self, v):
        """Apply minimum-image convention to displacement vector v (2,).
        Wraps only the periodic axes (periodic_x / periodic_y)."""
        v = v.copy()
        if self.periodic_x:
            v[0] -= self.Lx * np.round(v[0] / self.Lx)
        if self.periodic_y:
            v[1] -= self.Ly * np.round(v[1] / self.Ly)
        return v

    def _image_offsets(self):
        """Return list of (nx, ny) periodic image offsets to check.
        For non-periodic axes only nx/ny=0 is included.
        Offsets beyond ±1 are never needed when Lx, Ly > threshold."""
        nx_vals = [-1, 0, 1] if self.periodic_x else [0]
        ny_vals = [-1, 0, 1] if self.periodic_y else [0]
        return [(nx, ny) for nx in nx_vals for ny in ny_vals]

    def _pair_threshold(self, pA, pB):
        """Contact+skin threshold for pair (pA, pB) using per-particle R0 and skin."""
        return self._R0_arr[pA] + self._R0_arr[pB] + 0.5 * (self._skin_arr[pA] + self._skin_arr[pB])

    def _level1_pairs(self, x_cm):
        """
        Return list of (pA, pB) pairs with pA < pB that have at least one
        periodic image within contact+skin range (center-to-center).
        Checks all relevant periodic image offsets so that pairs near Lx/2
        (or Ly/2) are caught even when both direct and periodic paths are active.
        """
        pairs = []
        offsets = self._image_offsets()
        for pA in range(self.P):
            for pB in range(pA + 1, self.P):
                threshold = self._pair_threshold(pA, pB)
                raw = x_cm[pB] - x_cm[pA]
                for nx, ny in offsets:
                    d_vec = raw + np.array([nx * self.Lx, ny * self.Ly])
                    if np.linalg.norm(d_vec) < threshold:
                        pairs.append((pA, pB))
                        break
        return pairs

    # ── level 2 + 3: contact-normal registration and fill ────────────────────

    def _contact_facing_edge(self, p_idx, n_contact, theta):
        """
        Index of the edge on particle p_idx whose outward normal is closest
        to n_contact direction. O(1) via angle arithmetic.

        Edge k of a regular N-gon has reference normal at angle:
          phi_k = 2*pi*k/N + pi/2   (90° ahead of the k-th node)
        In the world frame: phi_k + theta[p_idx].
        The facing edge is argmin |angle_diff(phi_k_world, angle(n_contact))|.
        """
        angle_contact = np.arctan2(n_contact[1], n_contact[0])
        # Outward normal of edge k points at 2π*(k+0.5)/N + theta (world frame).
        k_arr      = np.arange(self.N)
        phi_world  = 2.0 * np.pi * (k_arr + 0.5) / self.N + theta[p_idx]
        # Signed angular distance (wrapped to [-pi, pi])
        diff       = (angle_contact - phi_world + np.pi) % (2 * np.pi) - np.pi
        return int(np.argmin(np.abs(diff)))

    def _fill_pair(self, pA, pB, x_cm, theta):
        """
        Fill CapCandidates rows for all periodic images of (pA, pB) within range.

        Iterates over every image offset (nx, ny) so that when both the direct
        and a periodic-image contact are active simultaneously (e.g., when CM
        separation is close to Lx/2), both contact zones get candidates.

        One-directional: only pA rows are written.  The force kernel applies
        Newton 3rd law via scatter, so pB rows are not needed.
        """
        N         = self.N
        dj        = self.dj
        threshold = self._pair_threshold(pA, pB)
        raw       = x_cm[pB] - x_cm[pA]

        for nx, ny in self._image_offsets():
            d_vec = raw + np.array([nx * self.Lx, ny * self.Ly])
            dist  = np.linalg.norm(d_vec)
            if dist < 1e-12 or dist >= threshold:
                continue
            n_AB = d_vec / dist

            # Facing edge on A (toward this image of B) and on B (toward A)
            jA = self._contact_facing_edge(pA,  n_AB, theta)
            jB = self._contact_facing_edge(pB, -n_AB, theta)

            # Window of candidate edges on B seen by edges near jA on A
            for de in range(-dj, dj + 1):
                eA  = (jA + de) % N
                row = pA * N + eA
                slots_used = int(np.sum(self.CapCandidates[row] != 0))
                for dc in range(-dj, dj + 1):
                    eB       = (jB + dc) % N
                    cand_idx = pB * N + eB + 1   # +1 because 0 is ghost
                    if slots_used >= self.E:
                        print(f"WARNING: row {row} full (E={self.E}). Skipping candidate.")
                        break
                    if cand_idx not in self.CapCandidates[row]:
                        self.CapCandidates[row, slots_used] = cand_idx
                        slots_used += 1

    # ── debug tools ───────────────────────────────────────────────────────────

    def check_missing_contacts(self, x_all, r_contact):
        """
        Brute-force scan: for every edge pair (eA, eB) between neighboring
        particles, check if the closest approach < r_contact.  Report any
        such pair absent from CapCandidates.

        x_all     : (P, N, 2) float64 — current node positions
        r_contact : float — contact detection radius (= r_cA + r_cB)

        Returns list of (eA_global, eB_global, gap) for missing contacts.
        """
        N   = self.N
        missing = []
        for pA in range(self.P):
            for pB in range(pA + 1, self.P):
                xA      = x_all[pA]
                xB      = x_all[pB]
                xA_next = np.roll(xA, -1, axis=0)
                xB_next = np.roll(xB, -1, axis=0)
                for eA in range(N):
                    for eB in range(N):
                        # Closest approach between two segments
                        gap = _segment_segment_min_dist(
                            xA[eA], xA_next[eA], xB[eB], xB_next[eB]) - r_contact
                        if gap < 0.0:
                            gA = pA * N + eA + 1
                            gB = pB * N + eB + 1
                            # Check if gB is in row gA-1 of CapCandidates
                            row_A = gA - 1
                            if gB not in self.CapCandidates[row_A]:
                                missing.append((gA, gB, gap))
        return missing

    def candidate_utilization(self):
        """Fraction of E slots that are non-ghost, per row."""
        used = np.sum(self.CapCandidates != 0, axis=1)
        return used / self.E

    def contact_count_histogram(self):
        """Number of active candidates per edge (non-zero slots)."""
        return np.sum(self.CapCandidates != 0, axis=1)

    def save_corpus(self):
        """Save logged input/output corpus to log_dir for Phase 3.4 C++ validation."""
        if self.log_dir is None or not self._corpus:
            return
        path = os.path.join(self.log_dir, 'candidacy_corpus.json')
        with open(path, 'w') as f:
            json.dump(self._corpus, f)
        print(f"  Corpus saved: {len(self._corpus)} entries → {path}")

    @staticmethod
    def load_corpus(path):
        with open(path) as f:
            return json.load(f)


# ── segment–segment minimum distance helper ──────────────────────────────────

def _segment_segment_min_dist(a0, a1, b0, b1):
    """
    Minimum distance between segments a0-a1 and b0-b1 (2D).
    Used only in debug check_missing_contacts().
    """
    def _pt_seg(p, s0, s1):
        d = s1 - s0
        t = np.dot(p - s0, d) / max(np.dot(d, d), 1e-30)
        t = np.clip(t, 0, 1)
        return np.linalg.norm(p - (s0 + t * d))

    # Check all four endpoint-to-segment distances and segment midpoint
    dists = [
        _pt_seg(a0, b0, b1),
        _pt_seg(a1, b0, b1),
        _pt_seg(b0, a0, a1),
        _pt_seg(b1, a0, a1),
    ]
    # Also check closest point on each segment to midpoint of the other
    am = 0.5 * (a0 + a1)
    bm = 0.5 * (b0 + b1)
    dists.append(_pt_seg(am, b0, b1))
    dists.append(_pt_seg(bm, a0, a1))
    return min(dists)
