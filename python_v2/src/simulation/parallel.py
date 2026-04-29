"""
parallel.py — run N independent Systems concurrently on the GPU via
XLA-fused step calls.

The N count is fixed at trace time: each unique N pays one ~10s XLA compile
cost (cached for the rest of the session).  Beyond N, only `sample_every`,
`cand_check_interval`, and the per-system tensor *shapes* affect the trace
key — so a 4-experiment Bo sweep at P=300 reuses the trace across chunks.

Hard requirements across all systems in one parallel_run call:
  - same dtype (fp64 throughout)
  - same n_steps and sample_every (set by the call itself)

Soft requirement (efficiency only — XLA SIMD-fusion is best when shapes match):
  - same P, N_nodes
Mismatched shapes still work; XLA emits separate sub-kernels per branch and
they run concurrently across SMs but without the SIMD vectorisation gain.

Free to vary per system: Bo / γ / κ / Oh / ν / dt / walls / box / periodic
config / initial state / particle type / E_candidates.
"""
from __future__ import annotations

import time
import numpy as np
import tensorflow as tf

from src.simulation.tf_sim import step_full_tf, DTYPE, NP_DTYPE


# ──────────────────────────────────────────────────────────────────────────────
# CPU-side per-system candidacy driver
# ──────────────────────────────────────────────────────────────────────────────

class _CandidacyDriver:
    """Tracks last (x_cm, theta) and rebuilds the C++ candidacy when the skin
    threshold is exceeded — the same logic as run_simulation_tf, but pulled
    out of the TF graph so the GPU side stays jit-compilable."""

    def __init__(self, system):
        self.mgr     = system._cm_mgr
        R0_arr       = np.array([p.R0 for p in system._particles])
        self.R0_max  = float(np.max(R0_arr))
        self.skin_th = float(system._skin * np.mean(R0_arr) * 0.5)
        self.last_xc = system._state['x_cm'].numpy().copy()
        self.last_th = system._state['theta'].numpy().copy()
        self.n_updates = 0
        self.n_checks  = 0
        # Force one initial update so CapCandidates is current
        self.mgr.update(np.ascontiguousarray(self.last_xc, dtype=np.float64),
                        np.ascontiguousarray(self.last_th, dtype=np.float64))

    def maybe_update(self, x_cm_np, theta_np):
        self.n_checks += 1
        dx  = float(np.max(np.linalg.norm(x_cm_np - self.last_xc, axis=1)))
        dth = float(np.max(np.abs(theta_np - self.last_th)))
        if (dx + self.R0_max * dth) > self.skin_th:
            self.mgr.update(np.ascontiguousarray(x_cm_np, dtype=np.float64),
                            np.ascontiguousarray(theta_np, dtype=np.float64))
            self.last_xc = x_cm_np.copy()
            self.last_th = theta_np.copy()
            self.n_updates += 1
            return True
        return False


# ──────────────────────────────────────────────────────────────────────────────
# GPU-side fused N-system step runner (jit-compiled, no Python callbacks)
# ──────────────────────────────────────────────────────────────────────────────

# Cache: maps N → tf.function.  Each unique N gets compiled once per session.
_RUNNER_CACHE: dict[int, callable] = {}


def _make_fused_runner(N: int):
    """Return a tf.function that runs `n_steps` GPU-side steps for N systems
    in lockstep, fully XLA-compiled (no py_function inside).  Cached per N.

    The function takes flat per-system tensors and returns updated state +
    per-system max closing-rate observed during the chunk.
    """
    if N in _RUNNER_CACHE:
        return _RUNNER_CACHE[N]

    @tf.function(jit_compile=True)
    def runner(states_flat, cands_flat, params_flat, prim_flat,
               n_steps, step_offset_flat):
        """
        Args (all flat lists/tuples of length N):
          states_flat  : list of state dicts (without 't', 'step')
          cands_flat   : list of (K, E) int32 CapCandidates tensors
          params_flat  : list of params dicts
          prim_flat    : list of prim_data dicts (or None per slot)
          n_steps      : scalar int32 — number of steps to run this chunk
          step_offset_flat : list of scalar int32 — global step index per sys

        Returns:
          new_states_flat  : list of updated state dicts
          max_ratios_flat  : list of scalar — per-system max closing-rate
        """
        dtype = DTYPE
        max_init = [tf.constant(0.0, dtype=dtype) for _ in range(N)]
        # local step counter in [0, n_steps)
        local_init = tf.constant(0, dtype=tf.int32)

        def body(local, states, maxes):
            new_states = []
            new_maxes  = []
            for i in range(N):
                t_i = (tf.cast(step_offset_flat[i], dtype) +
                       tf.cast(local, dtype)) * params_flat[i]['_dt_tf']
                ns_i, m_i = step_full_tf(
                    states[i], cands_flat[i],
                    params_flat[i]['_dt_tf'],
                    params_flat[i]['_alpha_tf'],
                    params_flat[i]['_g_tf'],
                    params_flat[i],
                    t=t_i,
                    prim_data=prim_flat[i],
                )
                new_states.append(ns_i)
                new_maxes.append(tf.maximum(maxes[i], m_i['max_closing_ratio']))
            return local + 1, new_states, new_maxes

        def cond(local, states, maxes):
            return local < n_steps

        _, final_states, final_maxes = tf.while_loop(
            cond, body,
            loop_vars=(local_init, states_flat, max_init),
            parallel_iterations=1,
        )
        return final_states, final_maxes

    _RUNNER_CACHE[N] = runner
    return runner


# ──────────────────────────────────────────────────────────────────────────────
# Public driver
# ──────────────────────────────────────────────────────────────────────────────

def parallel_run(systems, n_steps, sample_every=None,
                 cand_check_interval=10, verbose=True,
                 record_initial=True):
    """Run N independent Systems concurrently on the GPU.

    Each System advances exactly `n_steps` steps; sim time per system depends
    on its own dt.  Every `sample_every` steps each system's `frames` is
    appended a snapshot.  The per-system stability watchdog is updated.

    Returns
    -------
    diag : dict  — aggregate diagnostics including per-system ρ_max and
                   total wall time.
    """
    if not systems:
        raise ValueError("parallel_run: empty systems list")

    N = len(systems)
    if sample_every is None:
        sample_every = n_steps
    if sample_every <= 0 or n_steps % sample_every != 0:
        raise ValueError(
            f"n_steps={n_steps} must be a positive multiple of "
            f"sample_every={sample_every}")
    if sample_every % cand_check_interval != 0:
        raise ValueError(
            f"sample_every={sample_every} must be a multiple of "
            f"cand_check_interval={cand_check_interval}")

    # Validate cross-system constraints
    dtype = systems[0]._state['x_all'].dtype
    for s in systems[1:]:
        if s._state['x_all'].dtype != dtype:
            raise ValueError("parallel_run: all systems must share dtype")

    # Build candidacy drivers (also performs initial mgr.update)
    drivers = [_CandidacyDriver(s) for s in systems]
    runner  = _make_fused_runner(N)

    # Optional: record the initial state (step 0) before any integration
    if record_initial:
        for s in systems:
            s.frames.append(s.snapshot())

    cand_interval_tf = tf.constant(cand_check_interval, dtype=tf.int32)
    n_big_chunks   = n_steps // sample_every
    inner_per_big  = sample_every // cand_check_interval

    diag = {
        'n_steps':            n_steps,
        'sample_every':       sample_every,
        'cand_check_interval': cand_check_interval,
        'n_systems':          N,
        'rho_max_per_sys':    [0.0] * N,
        'cand_updates_per_sys': [0] * N,
        'wall_time_s':        0.0,
    }

    t_start = time.time()
    for big in range(n_big_chunks):
        for inner in range(inner_per_big):
            # 1. CPU-side: per-system candidacy refresh if skin breached
            cands = []
            for i, (s, drv) in enumerate(zip(systems, drivers)):
                xc = s._state['x_cm'].numpy()
                th = s._state['theta'].numpy()
                drv.maybe_update(xc, th)
                cands.append(tf.constant(drv.mgr.CapCandidates))

            # 2. GPU-side: cand_check_interval steps fused for all N systems
            phys_states = [
                {k: v for k, v in s._state.items() if k not in ('t', 'step')}
                for s in systems
            ]
            step_offsets = [tf.constant(int(s._state['step'].numpy()), dtype=tf.int32)
                            for s in systems]

            new_states, new_maxes = runner(
                phys_states, cands,
                [s._params for s in systems],
                [s._prim_data for s in systems],
                cand_interval_tf,
                step_offsets,
            )

            # 3. Write back per-system state with updated t and step
            for s, ns, nm in zip(systems, new_states, new_maxes):
                step_val = int(s._state['step'].numpy()) + cand_check_interval
                t_val    = step_val * float(s._dt)
                s._state = {
                    **ns,
                    't':    tf.constant(NP_DTYPE(t_val), dtype=DTYPE),
                    'step': tf.constant(step_val,        dtype=tf.int64),
                }
                rho = float(nm)
                if rho > s._max_disp_ratio_run:
                    s._max_disp_ratio_run = rho

        # End of big chunk: snapshot
        for s in systems:
            s.frames.append(s.snapshot())

        if verbose:
            rhos = [s._max_disp_ratio_run for s in systems]
            print(f"  parallel chunk {big+1}/{n_big_chunks}  "
                  f"step={int(systems[0]._state['step'].numpy())}  "
                  f"rho_max=[{', '.join(f'{r:.4f}' for r in rhos)}]",
                  flush=True)

    diag['wall_time_s'] = time.time() - t_start
    diag['rho_max_per_sys'] = [s._max_disp_ratio_run for s in systems]
    diag['cand_updates_per_sys'] = [d.n_updates for d in drivers]
    return diag
