"""
PairRegistrationCM — Python wrapper around the C++ pair-shell registered
candidacy manager (`_pair_registration_cm.PairRegistrationCM`).

Exposes the same surface as the production CandidacyManager so the System
+ tf_sim path can use either via `cm.update(x_cm, theta, x_all=...)`.

Async mode: when `async_mode=True`, update() returns immediately while a
worker thread runs the C++ PRCM (which releases the GIL during compute).
The worker writes to a back buffer; the next update() call commits that
buffer into CapCandidates before spawning the next worker. Net effect:
the kernel reads candidacy from one cycle ago, but PRCM cost is fully
hidden behind kernel work — usable to lower cand_check_interval (fresher
candidacy) without paying speed cost.

PRCM advantages over production sliding-fill CM:
  - Per-pair shell registration with O(1) index transport (no per-node search).
  - Tier-based surgical updates (L1/L2/Q) with bounded one-hop affected set.
  - Direct C++ implementation; ~6 ms/step at P=1290 vs ~5 ms for production.
  - Designed for narrow E (e.g. 9 with delta=4) without sliding-fill overflow
    warnings — keeps force-kernel cost low.

Construction requires per-particle R0 and r_c arrays (no cell-side
calibration), so build PRCM after particles are realized:

    from src.simulation.pair_registration_cm import PairRegistrationCM
    R0_arr  = np.array([p.R0  for p in sys._particles])
    r_c_arr = np.array([p.r_c for p in sys._particles])
    cm = PairRegistrationCM(P=len(sys._particles), N=sys._particles[0].N,
                              R0_arr=R0_arr, r_c_per_p=r_c_arr,
                              Lx=sys.Lx, Ly=sys.Ly, delta=4)

Or use System(candidacy_kind='prcm') — System constructs it for you at
initialize-time after particle realization.
"""
import threading
import numpy as np

from src.simulation import _pair_registration_cm as _ext


class PairRegistrationCM:
    """Drop-in replacement for production CandidacyManager. Same .update()
    + .CapCandidates surface so it works inside tf_sim.run_simulation_tf."""

    def __init__(self, P, N, R0_arr, r_c_per_p, Lx, Ly,
                 M1=20, M2=10, delta=4, skin=0.5,
                 L1_radius_scale=1.5, L2_radius_scale=1.0,
                 periodic_x=False, periodic_y=False,
                 L1_trigger_frac=0.5, L2_trigger_frac=0.1,
                 Q_trigger_frac=0.25,
                 async_mode=False,
                 cand_log_dir=None):
        self._cpp = _ext.PairRegistrationCM(
            P=P, N=N,
            R0_arr=np.asarray(R0_arr, dtype=np.float64),
            r_c_per_p=np.asarray(r_c_per_p, dtype=np.float64),
            Lx=float(Lx), Ly=float(Ly),
            M1=M1, M2=M2, delta=delta,
            L1_radius_scale=L1_radius_scale,
            L2_radius_scale=L2_radius_scale,
            skin=skin,
            periodic_x=periodic_x, periodic_y=periodic_y,
            L1_trigger_frac=L1_trigger_frac,
            L2_trigger_frac=L2_trigger_frac,
            Q_trigger_frac=Q_trigger_frac)
        # Mirror attributes the rest of the system reads
        self.P = P
        self.N = N
        self.E = self._cpp.E
        self.K = self._cpp.K
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self.skin = float(skin)
        self.R0 = float(np.mean(R0_arr))
        # CapCandidates buffer the force kernel reads — mirrors C++ side
        self.CapCandidates = np.asarray(self._cpp.CapCandidates).copy()
        self._initialized = False
        self.log_dir = cand_log_dir   # for compat (no logging implemented)
        self._n_updates = 0
        self._last_x_cm  = None
        self._last_theta = None
        # Async-mode plumbing (off by default — flip via set_async_mode(True))
        self._async_mode = bool(async_mode)
        self._worker_thread = None
        self._back_buffer = np.empty_like(self.CapCandidates)

    # ── Same public surface as CandidacyManager ────────────────────────────

    def needs_update(self, x_cm, theta):
        """Always returns True. PRCM has its own internal per-particle motion
        triggers — letting it run on every poll is cheap when nothing moves."""
        return True

    def set_async_mode(self, on):
        """Toggle async refresh. Safe to call any time. Spawned worker threads
        complete naturally; switching to sync first joins the in-flight worker."""
        if (not on) and (self._worker_thread is not None):
            self._worker_thread.join()
            self.CapCandidates[:] = self._back_buffer
            self._worker_thread = None
        self._async_mode = bool(on)

    def update(self, x_cm, theta, x_all=None):
        """Full Q refresh + L1/L2 surgical (per-particle motion triggered).

        Sync path: blocks until C++ PRCM completes, copies to CapCandidates.
        Async path: joins prior worker (committing its result), then spawns
                    new worker. Kernel sees one-cycle-old candidacy.

        Parameters
        ----------
        x_cm  : (P, 2)    float64
        theta : (P,)      float64
        x_all : (P, N, 2) float64  REQUIRED — PRCM needs per-node positions
        """
        if x_all is None:
            raise ValueError(
                "PairRegistrationCM.update() requires x_all (per-node positions). "
                "Make sure tf_sim's candidacy_step plumbs it through the callback.")
        x_cm  = np.ascontiguousarray(x_cm,  dtype=np.float64)
        theta = np.ascontiguousarray(theta, dtype=np.float64)
        x_all = np.ascontiguousarray(x_all, dtype=np.float64)

        if not self._async_mode:
            self._cpp.update_L1_L2_then_full_Q(x_all, x_cm, theta)
            self.CapCandidates[:] = self._cpp.CapCandidates
            self._n_updates += 1
            self._last_x_cm  = x_cm.copy()
            self._last_theta = theta.copy()
            return

        # Async: join prior worker → commit its result → spawn new worker
        if self._worker_thread is not None:
            self._worker_thread.join()
            self.CapCandidates[:] = self._back_buffer
            self._n_updates += 1
        # Snapshot args (worker reads these arrays after we return)
        x_all_local = x_all
        x_cm_local  = x_cm
        th_local    = theta
        def _worker():
            self._cpp.update_L1_L2_then_full_Q(x_all_local, x_cm_local, th_local)
            np.copyto(self._back_buffer, np.asarray(self._cpp.CapCandidates))
        self._worker_thread = threading.Thread(target=_worker)
        self._worker_thread.start()
        self._last_x_cm  = x_cm.copy()
        self._last_theta = theta.copy()

    # ── Diagnostics ─────────────────────────────────────────────────────────

    @property
    def diag_tier_L1(self): return self._cpp.diag_tier_L1
    @property
    def diag_tier_L2(self): return self._cpp.diag_tier_L2
    @property
    def diag_tier_Q(self):  return self._cpp.diag_tier_Q
    @property
    def diag_affected(self): return self._cpp.diag_affected
