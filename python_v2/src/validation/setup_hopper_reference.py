"""
Phase 2 — one-time test-bed setup.

For each N in {300, 600}:
  1. Build hopper, RSA seed at phi ≈ 0.42 (using outer capsule perimeter
     R0 + r_c for particle area in the phi calculation).
  2. Run 50_000 steps to reach a mid-flow steady state (jamming + raining + outlet exit).
     Save state via System.save_state() → hopper_N{P}_state_50k.npz.
  3. From the 50k checkpoint, run a further 2000 steps and save per-frame
     (t, x_cm, x_all, theta) → hopper_N{P}_ref_2000.npz.
     This is the gold-standard reference for verifying speedup iterations.

Estimated runtime at baseline (RTX 3080, TF 2.15):
  N=300 ~66 ms/step → 50k+2k ≈ 57 min
  N=600 ~124 ms/step → 50k+2k ≈ 1 h 47 min
"""
import os
import sys
import time
import argparse
import pathlib

import numpy as np

sys.path.insert(0, os.getcwd())

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.validation.profile_hopper import build_test_bed, size_hopper

OUT_DIR = pathlib.Path("results/profiling")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def build_reference(P, scale, n_steady=50_000, n_ref=2000, tag="emulsion_Bo005",
                    Bo=0.05, kappa=0.02, Oh=0.15, ptype='emulsion', nu=0.5,
                    E_candidates=None, dt_factor=None):
    label = f"N={P} {tag}"
    print("\n" + "=" * 72)
    print(f"REFERENCE SETUP — {label} (scale={scale:.3f}, dt_factor={dt_factor})")
    print("=" * 72)

    sys_h, info = build_test_bed(P=P, N_nodes=60, scale=scale, verbose=False,
                                 Bo=Bo, kappa=kappa, Oh=Oh, ptype=ptype, nu=nu,
                                 E_candidates=E_candidates)
    if dt_factor is not None:
        # Apply dt_factor BEFORE initialize ran the relax steps:
        # build_test_bed has already initialized at System default, so rescale.
        prev = sys_h._dt_factor
        dt_max = sys_h._dt / prev
        sys_h._dt_factor = float(dt_factor)
        sys_h._dt = float(dt_factor) * dt_max
        import tensorflow as tf
        from src.simulation.tf_sim import NP_DTYPE, DTYPE
        sys_h._dt_tf = tf.constant(NP_DTYPE(sys_h._dt), dtype=DTYPE)
        sys_h._params['_dt_tf'] = sys_h._dt_tf
        print(f"  dt_factor override: {prev} -> {dt_factor}, dt={sys_h._dt:.6f}")
    print(f"Geometry: W_RES={info['W_RES']:.2f} H_RES={info['H_RES']:.2f}"
          f"  LX={info['LX']:.2f} LY={info['LY']:.2f}")
    print(f"RSA: phi_outer = {info['phi_outer']:.4f}")

    # ── 50k steady-state run ─────────────────────────────────────────────
    print(f"\nRunning {n_steady} steps to mid-flow steady state...")
    t0 = time.time()
    chunk = min(5000, max(1, n_steady))
    n_chunks = (n_steady + chunk - 1) // chunk
    steps_done = 0
    for i in range(n_chunks):
        steps_this = min(chunk, n_steady - steps_done)
        sys_h.run(steps_this, sample_every=steps_this, verbose=False)
        steps_done += steps_this
        elapsed = time.time() - t0
        ms_per_step = elapsed * 1000 / steps_done if steps_done > 0 else 0
        eta = (n_steady - steps_done) * ms_per_step / 1000
        print(f"  chunk {i+1}/{n_chunks}: {steps_done}/{n_steady} steps  "
              f"({elapsed:.0f}s elapsed, {ms_per_step:.1f} ms/step, ETA {eta:.0f}s)")

    state_path = OUT_DIR / f"hopper_N{P}_{tag}_state_50k.npz"
    sys_h.save_state(state_path)
    state_size = state_path.stat().st_size / 1024
    print(f"Saved 50k state → {state_path}  ({state_size:.0f} KB)")

    # ── 2000-step reference run with frame logging ───────────────────────
    print(f"\nGenerating reference: {n_ref} steps from 50k checkpoint...")
    frames_t   = []
    frames_xcm = []
    frames_x   = []
    frames_th  = []

    def _ref_cb(s):
        snap = s.snapshot()
        frames_t.append(float(snap['t']))
        frames_xcm.append(snap['x_cm'].copy())
        frames_x.append(snap['x_all'].copy())
        frames_th.append(snap['theta'].copy())
        return {}

    t0 = time.time()
    sys_h.run(n_ref, sample_every=400, callback=_ref_cb,
              record_initial=True, verbose=False)
    print(f"Reference run done: {time.time()-t0:.0f}s, {len(frames_t)} frames")

    ref_path = OUT_DIR / f"hopper_N{P}_{tag}_ref_2000.npz"
    np.savez(ref_path,
             t=np.array(frames_t),
             x_cm=np.array(frames_xcm),
             x_all=np.array(frames_x),
             theta=np.array(frames_th),
             P=np.int64(P),
             scale=np.float64(scale),
             n_steady=np.int64(n_steady),
             n_ref=np.int64(n_ref))
    print(f"Saved 2000-step reference → {ref_path}  "
          f"({ref_path.stat().st_size/1024:.0f} KB)")
    return state_path, ref_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--p',       type=int, default=300)
    ap.add_argument('--n-steady', type=int, default=50_000)
    ap.add_argument('--n-ref',    type=int, default=2000)
    ap.add_argument('--tag',     type=str, default='emulsion_Bo005',
                    help='filename tag, e.g. emulsion_Bo005, emulsion_Bo010, elastic_nu05')
    ap.add_argument('--Bo',      type=float, default=0.05)
    ap.add_argument('--kappa',   type=float, default=0.02)
    ap.add_argument('--Oh',      type=float, default=0.15)
    ap.add_argument('--ptype',   type=str,   default='emulsion',
                    choices=['emulsion', 'elastic'])
    ap.add_argument('--nu',      type=float, default=0.5)
    ap.add_argument('--E',       type=int, default=None)
    ap.add_argument('--dt-factor', type=float, default=None,
                    help='dt_factor override for the setup (None → System default)')
    args = ap.parse_args()

    build_reference(P=args.p, scale=1.0,
                    n_steady=args.n_steady, n_ref=args.n_ref,
                    tag=args.tag, Bo=args.Bo, kappa=args.kappa,
                    Oh=args.Oh, ptype=args.ptype, nu=args.nu,
                    E_candidates=args.E, dt_factor=args.dt_factor)
    print("\nReference test-bed setup complete.")


if __name__ == "__main__":
    main()
