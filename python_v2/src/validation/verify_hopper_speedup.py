"""
Phase 2 — verify a code change against the reference test bed.

For each N in {300, 600}:
  1. Build hopper, restore_state from results/profiling/hopper_N{P}_state_50k.npz.
  2. Run 2000 steps; record per-frame (t, x_cm, x_all, theta).
  3. Diff against the reference at results/profiling/hopper_N{P}_ref_2000.npz.
     Report max |Δ| per state component. Pass: < 1e-10.
  4. Re-profile bucket times via profile_hopper.profile_run +
     profile_hopper.profile_tf_buckets, and print the same bucket table
     as the baseline profile so we can compare.

Run after every code change to confirm physics unchanged AND see speedup.
"""
import os
import sys
import time
import json
import argparse
import pathlib

import numpy as np

sys.path.insert(0, os.getcwd())

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.validation.profile_hopper import (
    build_test_bed, profile_run, profile_tf_buckets,
    render_bucket_table, print_table,
)

OUT_DIR = pathlib.Path("results/profiling")


def verify_one(P, scale, label, n_ref=2000):
    print("\n" + "=" * 72)
    print(f"VERIFY — {label} (scale={scale:.3f})")
    print("=" * 72)

    state_path = OUT_DIR / f"hopper_N{P}_state_50k.npz"
    ref_path   = OUT_DIR / f"hopper_N{P}_ref_2000.npz"
    if not state_path.exists() or not ref_path.exists():
        print(f"MISSING reference files: {state_path} / {ref_path}")
        print("Run setup_hopper_reference.py first.")
        return None

    sys_h, info = build_test_bed(P=P, N_nodes=60, scale=scale, verbose=False)
    sys_h.restore_state(state_path)
    print(f"Restored 50k checkpoint from {state_path}")
    print(f"  state shapes: x_cm={tuple(sys_h._state['x_cm'].shape)}  "
          f"x_all={tuple(sys_h._state['x_all'].shape)}")

    ref = np.load(ref_path)
    print(f"Loaded reference frames from {ref_path}: {ref['x_cm'].shape[0]} frames")

    # Run 2000 steps with the same callback as the reference
    new_t   = []
    new_xcm = []
    new_x   = []
    new_th  = []
    def _cb(s):
        snap = s.snapshot()
        new_t.append(float(snap['t']))
        new_xcm.append(snap['x_cm'].copy())
        new_x.append(snap['x_all'].copy())
        new_th.append(snap['theta'].copy())
        return {}

    t0 = time.time()
    sys_h.run(n_ref, sample_every=400, callback=_cb,
              record_initial=True, verbose=False)
    t_run = time.time() - t0
    print(f"\nReplay run: {t_run:.0f}s ({t_run*1000/n_ref:.2f} ms/step)")

    # Diff against reference (frame-by-frame)
    new_xcm = np.array(new_xcm)
    new_x   = np.array(new_x)
    new_th  = np.array(new_th)

    n_min = min(new_xcm.shape[0], ref['x_cm'].shape[0])
    dx_cm = np.abs(new_xcm[:n_min] - ref['x_cm'][:n_min]).max()
    dx_all = np.abs(new_x[:n_min] - ref['x_all'][:n_min]).max()
    d_th = np.abs(new_th[:n_min] - ref['theta'][:n_min]).max()

    pass_ = max(dx_cm, dx_all, d_th) < 1e-10
    print(f"\n── Correctness diff (frames 0..{n_min-1}) ──")
    print(f"  max |Δx_cm|  : {dx_cm:.3e}")
    print(f"  max |Δx_all| : {dx_all:.3e}")
    print(f"  max |Δθ|     : {d_th:.3e}")
    print(f"  {'PASS' if pass_ else 'FAIL'} (threshold 1e-10)")

    # Re-profile to see speedup
    print(f"\nRe-profiling on the same restored state...")
    sys_h2, info2 = build_test_bed(P=P, N_nodes=60, scale=scale, verbose=False)
    sys_h2.restore_state(state_path)
    run_t  = profile_run(sys_h2, n_warmup=200, n_steps=2000, sample_every=400)
    tf_t   = profile_tf_buckets(sys_h2)

    result = dict(label=label, P=P, scale=scale, geometry=info,
                  run=run_t, tf_buckets=tf_t)
    rows, tot, meta = render_bucket_table(result)
    print_table(label, rows, tot, meta)

    return dict(diff_passed=pass_,
                dx_cm=float(dx_cm), dx_all=float(dx_all), dtheta=float(d_th),
                replay_run_ms=t_run * 1000.0,
                replay_ms_per_step=t_run * 1000.0 / n_ref,
                profile=result)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p1', type=int, default=300)
    ap.add_argument('--p2', type=int, default=600)
    ap.add_argument('--tag', type=str, default='unknown',
                    help='Identifier for this iteration, e.g. "stage1_scatter"')
    args = ap.parse_args()

    res1 = verify_one(P=args.p1, scale=1.0, label=f"N={args.p1}")
    res2 = verify_one(P=args.p2, scale=float(np.sqrt(args.p2 / args.p1)),
                      label=f"N={args.p2}")

    out = OUT_DIR / f"verify_{args.tag}.json"
    with open(out, 'w') as f:
        json.dump({'tag': args.tag, args.p1: res1, args.p2: res2},
                  f, indent=2, default=str)
    print(f"\nSaved verify result → {out}")


if __name__ == "__main__":
    main()
