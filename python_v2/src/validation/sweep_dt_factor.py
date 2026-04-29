"""
Phase 2 — sweep dt_factor on the hopper test bed and watch for the new
displacement-stability watchdog warnings. Goal: find the largest dt_factor
that runs cleanly (no advisory/strong/critical fires, no NaN, no blow-up).

For each dt_factor, builds the hopper test bed, restores the 50k checkpoint,
overrides dt_factor (rebuilding dt and the integrator constants), runs
2000 steps, and reports:
  - max_disp_ratio observed
  - any watchdog warnings (captured via stdout)
  - NaN-free check on final state
  - mean y_cm shift (sanity: particles still falling/discharging)
  - wall time
"""
import os
import sys
import time
import argparse
import contextlib
import io
import pathlib

import numpy as np
import tensorflow as tf

sys.path.insert(0, os.getcwd())

import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.validation.profile_hopper import build_test_bed

OUT_DIR = pathlib.Path("results/profiling")


def run_one(P, scale, dt_factor, T_target=None, n_steps=2000,
            tag='emulsion_Bo005', Bo=0.05, kappa=0.02, Oh=0.15,
            ptype='emulsion', nu=0.5, E_candidates=None):
    """Run from 50k checkpoint at given dt_factor.

    If T_target is given, runs ceil(T_target / dt) steps so all dt_factors
    cover the SAME simulated time (correct way to compare physics).
    Otherwise runs n_steps fixed (different simulated times across dts).
    """
    state_path = OUT_DIR / f"hopper_N{P}_{tag}_state_50k.npz"
    if not state_path.exists():
        print(f"  missing {state_path}")
        return None

    # Build with dt_factor override (and matching physics knobs!)
    sys_h, info = build_test_bed(P=P, N_nodes=60, scale=scale, verbose=False,
                                 Bo=Bo, kappa=kappa, Oh=Oh, ptype=ptype, nu=nu,
                                 E_candidates=E_candidates)
    # sys_h._dt was set during initialize as sys_h._dt_factor * dt_max.
    # Recover dt_max and rescale to the requested dt_factor.
    prev_factor      = float(sys_h._dt_factor)
    dt_max_inferred  = sys_h._dt / prev_factor
    sys_h._dt_factor = float(dt_factor)
    new_dt           = float(dt_factor) * dt_max_inferred
    sys_h._dt        = new_dt
    from src.simulation.tf_sim import NP_DTYPE, DTYPE
    sys_h._dt_tf = tf.constant(NP_DTYPE(new_dt), dtype=DTYPE)

    sys_h.restore_state(state_path)
    # restore_state rewrites dt to the checkpoint's value — re-apply.
    # IMPORTANT: sys.run uses run_fast → run_simulation_tf which reads dt
    # from params['_dt_tf'], NOT self._dt_tf. Must update BOTH.
    sys_h._dt = new_dt
    sys_h._dt_tf = tf.constant(NP_DTYPE(new_dt), dtype=DTYPE)
    sys_h._params['_dt_tf'] = sys_h._dt_tf

    # Pick step count: matched simulated time if T_target given
    if T_target is not None:
        n_steps = max(1, int(np.ceil(T_target / new_dt)))

    # Initial state for sanity comparison
    y_cm_initial = float(sys_h._state['x_cm'].numpy()[:, 1].mean())

    # Capture stdout to count warnings
    buf = io.StringIO()
    t0 = time.time()
    nan_safe = True
    try:
        with contextlib.redirect_stdout(buf):
            sys_h.run(n_steps, sample_every=n_steps, verbose=False)
    except Exception as exc:
        elapsed = time.time() - t0
        return dict(P=P, dt_factor=dt_factor, dt=float(new_dt),
                    n_steps=n_steps, elapsed_s=elapsed,
                    error=str(exc), nan_safe=False)
    elapsed = time.time() - t0

    # Inspect final state
    x_cm  = sys_h._state['x_cm'].numpy()
    x_all = sys_h._state['x_all'].numpy()
    nan_safe = (not np.isnan(x_cm).any()) and (not np.isnan(x_all).any())
    y_cm_final = float(x_cm[:, 1].mean())
    sim_time   = n_steps * float(new_dt)

    log = buf.getvalue()
    # Watchdog signatures (only the displacement watchdog, not candidacy warnings)
    n_advisory = log.count('NOTE (watchdog)')
    n_strong   = log.count('WARNING (watchdog)')
    n_critical = log.count('CRITICAL')   # watchdog-only label
    n_cand_warn = log.count('CandidacyManager') + log.count('Skipping candidate')

    return dict(
        P=P, dt_factor=dt_factor, dt=float(new_dt),
        n_steps=n_steps, elapsed_s=elapsed, sim_time=sim_time,
        nan_safe=nan_safe,
        y_cm_initial=y_cm_initial, y_cm_final=y_cm_final,
        x_cm_final=x_cm,                             # (P, 2) for per-particle diff
        max_disp_ratio_run=sys_h._max_disp_ratio_run,
        n_advisory=n_advisory, n_strong=n_strong, n_critical=n_critical,
        n_cand_warn=n_cand_warn,
        log_excerpt=log[:2000],
    )


def print_row(r, ref=None):
    if r is None or 'error' in r:
        print(f"  ERROR  {r.get('error', '')}")
        return
    flag = ''
    if not r['nan_safe']:
        flag = '  ⚠ NaN'
    elif r['n_critical'] > 0:
        flag = '  CRITICAL'
    elif r['n_strong'] > 0:
        flag = f"  STRONG×{r['n_strong']}"
    elif r['n_advisory'] > 0:
        flag = f"  advisory×{r['n_advisory']}"

    # Per-particle drift relative to reference (matched simulated time)
    drift_str = ""
    if ref is not None and r['nan_safe']:
        diff = r['x_cm_final'] - ref['x_cm_final']                  # (P, 2)
        dist = np.linalg.norm(diff, axis=1)                         # (P,)
        drift_str = (f"  drift mean={dist.mean():.4f} "
                     f"std={dist.std():.4f} max={dist.max():.4f}")

    print(f"  dt_factor={r['dt_factor']:.2f}  dt={r['dt']:.5f}  "
          f"sim_t={r['sim_time']:.3f}  steps={r['n_steps']}  "
          f"rho_contact={r['max_disp_ratio_run']:.4f}  "
          f"({r['elapsed_s']:.1f}s){flag}{drift_str}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p', type=int, default=300)
    ap.add_argument('--T-target', type=float, default=2.0,
                    help='Target simulated time per dt_factor (matched-T comparison)')
    ap.add_argument('--only', type=float, default=None,
                    help='Run only this dt_factor (used by driver subprocess)')
    ap.add_argument('--tag',     type=str, default='emulsion_Bo005')
    ap.add_argument('--Bo',      type=float, default=0.05)
    ap.add_argument('--kappa',   type=float, default=0.02)
    ap.add_argument('--Oh',      type=float, default=0.15)
    ap.add_argument('--ptype',   type=str, default='emulsion',
                    choices=['emulsion', 'elastic'])
    ap.add_argument('--nu',      type=float, default=0.5)
    ap.add_argument('--E',       type=int, default=None,
                    help='E_candidates override (None → System default)')
    args = ap.parse_args()

    scale = 1.0 if args.p == 300 else float(np.sqrt(args.p / 300.0))

    print(f"\nSweeping dt_factor on hopper N={args.p} (scale={scale:.3f}), "
          f"matched T={args.T_target} simulated time")
    print("Watchdog thresholds: advisory=0.05, strong=0.10, critical=0.20")
    print("Drift stats = ‖x_cm[NEW] − x_cm[REF=dt_factor=0.4]‖ per particle, after T sim-time")
    print("─" * 110)

    # Run only the requested dt_factor and dump x_cm to disk so a separate
    # process can do the comparison. Avoids accumulating JIT compilation
    # state across iterations (which has caused segfaults).
    if args.only is not None:
        f = args.only
        r = run_one(P=args.p, scale=scale, dt_factor=f, T_target=args.T_target,
                    tag=args.tag, Bo=args.Bo, kappa=args.kappa, Oh=args.Oh,
                    ptype=args.ptype, nu=args.nu, E_candidates=args.E)
        if r is None:
            return
        out = OUT_DIR / f"sweep_dt_N{args.p}_{args.tag}_factor{f:.2f}.npz"
        np.savez(out,
                 dt_factor=f, dt=r['dt'], n_steps=r['n_steps'],
                 sim_time=r['sim_time'], elapsed_s=r['elapsed_s'],
                 nan_safe=int(r['nan_safe']),
                 max_disp_ratio_run=r['max_disp_ratio_run'],
                 n_advisory=r['n_advisory'], n_strong=r['n_strong'],
                 n_critical=r['n_critical'],
                 x_cm_final=r['x_cm_final'])
        print(f"saved {out}")
        return

    # Driver mode: spawn a child process per dt_factor
    import subprocess
    factors = [0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0]
    ref_xcm = None
    for f in factors:
        out = OUT_DIR / f"sweep_dt_N{args.p}_{args.tag}_factor{f:.2f}.npz"
        cmd = [sys.executable, __file__, '--p', str(args.p),
               '--T-target', str(args.T_target), '--only', str(f),
               '--tag', args.tag, '--Bo', str(args.Bo),
               '--kappa', str(args.kappa), '--Oh', str(args.Oh),
               '--ptype', args.ptype, '--nu', str(args.nu)]
        if args.E is not None:
            cmd += ['--E', str(args.E)]
        subprocess.run(cmd, check=False)
        if not out.exists():
            print(f"  dt_factor={f}: no output (run failed)")
            continue
        d = np.load(out)
        x_cm = d['x_cm_final']
        nan_safe = bool(d['nan_safe'])
        max_disp = float(d['max_disp_ratio_run'])

        if f == 0.4:
            ref_xcm = x_cm

        # Drift stats
        drift_str = ""
        if ref_xcm is not None and nan_safe:
            diff = x_cm - ref_xcm
            dist = np.linalg.norm(diff, axis=1)
            drift_str = (f"  drift mean={dist.mean():.4f} "
                         f"std={dist.std():.4f} max={dist.max():.4f}")

        flag = ''
        if not nan_safe:
            flag = '  ⚠ NaN'
        elif d['n_critical'] > 0:
            flag = '  CRITICAL'
        elif d['n_strong'] > 0:
            flag = f"  STRONG×{int(d['n_strong'])}"
        elif d['n_advisory'] > 0:
            flag = f"  advisory×{int(d['n_advisory'])}"

        print(f"  dt_factor={f:.2f}  dt={float(d['dt']):.5f}  "
              f"sim_t={float(d['sim_time']):.3f}  steps={int(d['n_steps'])}  "
              f"rho_contact={max_disp:.4f}  "
              f"({float(d['elapsed_s']):.1f}s){flag}{drift_str}")

        if not nan_safe or d['n_critical'] > 0:
            print(f"  Stopping sweep — instability at dt_factor={f}")
            break

    print("─" * 110)


if __name__ == "__main__":
    main()
