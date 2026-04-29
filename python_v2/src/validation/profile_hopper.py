"""
Phase 2 — bucket profiling on the Test F hopper test bed.

Goal: rank where compute time is spent at N=300 and N=600 droplets.
Test bed = emulsion droplets in a hopper at phi_outer ≈ 0.4, no swell
(RSA seeding only), N_nodes=60, Oh=0.15, kappa=0.02.

Buckets:
  C++ candidacy        — neighbor-list rebuild (instrumented update())
  TF capsule-capsule   — inter_capsule_forces_tf (isolated)
  TF primitive         — primitive_forces_tf (isolated)
  TF node-only forces  — internal_forces_tf (isolated)
  TF k_reg             — k_reg_forces_tf (isolated)
  TF drag              — diff: step_rb_tf with vs without xi_drag
  TF time integration  — step_rb_no_drag - internal_forces - k_reg
  Python bridge        — sample_every callback / snapshot
  Launch + unattributed — total - sum(above)
"""
import os
import sys
import time
import json
import pathlib
import argparse
import numpy as np

sys.path.insert(0, os.getcwd())

import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)

from src.epd.particles import ParticleSpec
from src.epd.objects import HopperRegion
from src.epd.system import System
from src.simulation.tf_sim import (
    internal_forces_tf, k_reg_forces_tf,
    inter_capsule_forces_tf, primitive_forces_tf,
    step_rb_tf,
)

OUT_DIR = pathlib.Path("results/profiling")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ────────────────────────────────────────────────────────────────────────────
# Geometry sizing
# ────────────────────────────────────────────────────────────────────────────
def size_hopper(P, R0_mean=1.0, N_nodes=60, phi_target=0.40,
                W_OUT=4.0, theta_deg=30.0, scale=1.0):
    """Use the OUTER capsule perimeter (R0 + r_c) to set particle area
    when sizing the hopper to a given phi target.
    """
    L0    = 2.0 * np.pi * R0_mean / N_nodes
    r_c   = L0                                      # contact radius
    r_eff = R0_mean + r_c
    eff_area_per_p = np.pi * r_eff**2
    W_RES = 24.0 * scale
    funnel_h = (W_RES - W_OUT) / 2.0 / np.tan(np.radians(theta_deg))
    funnel_area = (W_RES + W_OUT) / 2.0 * funnel_h
    target_area = P * eff_area_per_p / phi_target
    H_RES = max(8.0, (target_area - funnel_area) / W_RES)
    return dict(W_RES=W_RES, W_OUT=W_OUT, THETA_DEG=theta_deg,
                H_RES=H_RES, funnel_h=funnel_h,
                target_area=target_area, funnel_area=funnel_area,
                r_eff=r_eff, eff_area_per_p=eff_area_per_p)


def build_test_bed(P, N_nodes=60, scale=1.0, seed=42, verbose=True,
                   Bo=0.05, kappa=0.02, Oh=0.15, ptype='emulsion', nu=0.5,
                   E_candidates=None):
    """Build the hopper test bed with the given physics knobs.

    ptype='emulsion':  uses kappa, Oh.  Bo = g (with γ=ρ=R0=1).
    ptype='elastic' :  uses nu.         g = Bo, no Oh, no kappa.
    """
    geom = size_hopper(P, N_nodes=N_nodes, scale=scale)
    LX = geom['W_RES'] + 6.0
    Y_BOT = 2.0
    LY = Y_BOT + geom['funnel_h'] + geom['H_RES'] + 8.0  # 8 R0 headroom
    X_C = LX / 2.0

    sys_kwargs = dict(periodic_x=False, periodic_y=False, g=float(Bo))
    if E_candidates is not None:
        sys_kwargs['E_candidates'] = int(E_candidates)
    sys_h = System(LX, LY, **sys_kwargs)
    hopper = HopperRegion(X_C, geom['W_OUT'], geom['W_RES'],
                          geom['THETA_DEG'], geom['H_RES'], Y_BOT)
    sys_h.add_object(hopper)
    if ptype == 'emulsion':
        spec = ParticleSpec(count=P, type='emulsion',
                            kappa=kappa, Oh=Oh, N_nodes=N_nodes)
    elif ptype == 'elastic':
        spec = ParticleSpec(count=P, type='elastic',
                            nu=nu, N_nodes=N_nodes)
    else:
        raise ValueError(f"unknown ptype {ptype!r}")
    sys_h.add_particles(spec)
    sys_h.initialize(phi_target=0.40, seed=seed, verbose=verbose,
                     relax_only=True, n_relax_init=50)
    info = dict(geom)
    info.update(P=P, N_nodes=N_nodes, scale=scale,
                LX=LX, LY=LY, X_C=X_C, Y_BOT=Y_BOT,
                Bo=float(Bo), ptype=ptype,
                kappa=kappa if ptype == 'emulsion' else None,
                Oh=Oh if ptype == 'emulsion' else None,
                nu=nu if ptype == 'elastic' else None,
                phi_outer=float(sys_h.phi_outer))
    return sys_h, info


# ────────────────────────────────────────────────────────────────────────────
# Timing helpers
# ────────────────────────────────────────────────────────────────────────────
def _sync(out):
    if isinstance(out, (tuple, list)):
        for o in out:
            _sync(o)
    elif isinstance(out, dict):
        for o in out.values():
            _sync(o)
    elif hasattr(out, 'numpy'):
        out.numpy()


def time_call(fn, args, n_warmup=3, n_iter=20):
    for _ in range(n_warmup):
        _sync(fn(*args))
    t0 = time.time()
    for _ in range(n_iter):
        _sync(fn(*args))
    return (time.time() - t0) / n_iter * 1000.0


def install_cm_mgr_timer(sys_h):
    """Wrap the CandidacyManager.update method with a timer."""
    cm = sys_h._cm_mgr
    sys_h._cm_update_total_s = 0.0
    sys_h._cm_update_calls   = 0
    orig_update = cm.update

    def _timed_update(*a, **kw):
        t0 = time.time()
        out = orig_update(*a, **kw)
        sys_h._cm_update_total_s += (time.time() - t0)
        sys_h._cm_update_calls   += 1
        return out

    cm.update = _timed_update


# ────────────────────────────────────────────────────────────────────────────
# Bucket profiling
# ────────────────────────────────────────────────────────────────────────────
def profile_run(sys_h, n_warmup=200, n_steps=2000, sample_every=400):
    """Time the full sys.run with bridge + C++ isolation."""
    install_cm_mgr_timer(sys_h)

    # Warmup
    t0 = time.time()
    sys_h.run(n_warmup, sample_every=n_warmup, verbose=False)
    t_warmup = time.time() - t0

    # Reset C++ counters after warmup so measurement isolates the steady-state
    sys_h._cm_update_total_s = 0.0
    sys_h._cm_update_calls = 0

    # Measurement
    bridge_times = []
    def _cb(s):
        tcb0 = time.time()
        s.snapshot()
        bridge_times.append(time.time() - tcb0)
        return {}

    diag_n_before = len(sys_h.diag)
    t0 = time.time()
    sys_h.run(n_steps, sample_every=sample_every, callback=_cb,
              record_initial=False, verbose=False)
    t_run_total = time.time() - t0

    diag_run = sys_h.diag[diag_n_before:]
    return dict(
        n_warmup=n_warmup, n_steps=n_steps, sample_every=sample_every,
        t_warmup_s=t_warmup,
        t_run_total_ms=t_run_total * 1000.0,
        ms_per_step=(t_run_total * 1000.0) / n_steps,
        bridge_total_ms=sum(bridge_times) * 1000.0,
        bridge_n_calls=len(bridge_times),
        cm_update_total_ms=sys_h._cm_update_total_s * 1000.0,
        cm_update_calls=sys_h._cm_update_calls,
        n_chunks=len(diag_run),
        n_cand_checks=sum(d['n_cand_checks']  for d in diag_run),
        n_cand_updates=sum(d['n_cand_updates'] for d in diag_run),
    )


def profile_tf_buckets(sys_h):
    """Isolated timing of each TF function with the current state."""
    state  = sys_h._state
    params = sys_h._params
    g_val  = float(sys_h._g_val)

    x_all = state['x_all']
    P     = int(x_all.shape[0])
    Nn    = int(x_all.shape[1])

    # Flats — same way step_full_tf builds them.
    r_c_flat = tf.repeat(params['r_c_per_p'], Nn)
    k_c_flat = tf.repeat(params['k_c_per_p'], Nn)
    L0_flat  = tf.repeat(params['L0'],        Nn)

    # CapCandidates as TF int32 tensor
    CapCand = tf.constant(sys_h._cm_mgr.CapCandidates, dtype=tf.int32)

    prim_data = sys_h._prim_data
    t_now     = tf.constant(0.0, dtype=tf.float64)

    timings = {}

    # 1. Internal (node-only) forces
    timings['internal_forces'] = time_call(
        internal_forces_tf, (x_all, params, g_val))

    # 2. k_reg
    timings['k_reg'] = time_call(
        k_reg_forces_tf, (x_all, params))

    # 3. Inter-capsule contact
    is_periodic = sys_h._config.get('periodic_x', False) or sys_h._config.get('periodic_y', False)
    if is_periodic:
        box_tf = tf.constant([sys_h.Lx, sys_h.Ly], dtype=tf.float64)
        timings['inter_capsule'] = time_call(
            inter_capsule_forces_tf,
            (x_all, CapCand, r_c_flat, k_c_flat, L0_flat, box_tf))
    else:
        timings['inter_capsule'] = time_call(
            inter_capsule_forces_tf,
            (x_all, CapCand, r_c_flat, k_c_flat, L0_flat))

    # 4. Primitive (wall) forces
    timings['primitive'] = time_call(
        primitive_forces_tf, (x_all, t_now, prim_data, r_c_flat, k_c_flat, L0_flat))

    # 5. step_rb_tf with drag (current Oh)
    f_contact = tf.zeros_like(x_all)
    dt_f = tf.constant(sys_h._dt, dtype=tf.float64)
    alpha_f = tf.constant(0.0, dtype=tf.float64)
    g_f = tf.constant(g_val, dtype=tf.float64)
    timings['step_rb_with_drag'] = time_call(
        step_rb_tf, (state, f_contact, dt_f, alpha_f, g_f, params))

    # 6. step_rb_tf with drag disabled (Oh=0)
    saved_xi = params['xi_drag_per_p']
    params_no_drag = dict(params)
    params_no_drag['xi_drag_per_p'] = tf.zeros_like(saved_xi)
    timings['step_rb_no_drag'] = time_call(
        step_rb_tf, (state, f_contact, dt_f, alpha_f, g_f, params_no_drag))

    # Derived buckets
    timings['drag'] = max(0.0, timings['step_rb_with_drag']
                                - timings['step_rb_no_drag'])
    timings['integration_only'] = max(0.0, timings['step_rb_no_drag']
                                          - timings['internal_forces']
                                          - timings['k_reg'])
    return timings


def profile_one_size(P, scale, label, n_steps=2000, sample_every=400):
    print("\n" + "=" * 72)
    print(f"PROFILING — {label}: P={P}, scale={scale:.3f}")
    print("=" * 72)

    sys_h, info = build_test_bed(P=P, N_nodes=60, scale=scale)
    print(f"\nGeometry:")
    print(f"  W_RES={info['W_RES']:.2f}  H_RES={info['H_RES']:.2f}"
          f"  funnel_h={info['funnel_h']:.2f}")
    print(f"  LX={info['LX']:.2f}  LY={info['LY']:.2f}")
    print(f"  RSA result: phi_outer = {info['phi_outer']:.4f}  (target 0.40)")

    print(f"\nRun timing...")
    run_t = profile_run(sys_h, n_warmup=200, n_steps=n_steps, sample_every=sample_every)
    print(f"  warmup       : {run_t['t_warmup_s']:.1f} s")
    print(f"  measurement  : {run_t['t_run_total_ms']/1000:.2f} s "
          f"({run_t['ms_per_step']:.3f} ms/step)")
    print(f"  bridge       : {run_t['bridge_total_ms']:.1f} ms "
          f"({run_t['bridge_n_calls']} calls)")
    print(f"  C++ update   : {run_t['cm_update_total_ms']:.1f} ms "
          f"({run_t['cm_update_calls']} updates)")
    print(f"  diag: {run_t['n_cand_checks']} checks, {run_t['n_cand_updates']} updates")

    print(f"\nIsolated TF function timing (per-call ms)...")
    tf_t = profile_tf_buckets(sys_h)
    for name, ms in tf_t.items():
        print(f"  {name:24s} : {ms:7.3f} ms")

    return dict(label=label, P=P, scale=scale, geometry=info,
                run=run_t, tf_buckets=tf_t)


def render_bucket_table(result):
    """Build a ranked bucket table. TF bucket times are isolated-call
    measurements; we rescale them to fit (total - cpp - bridge) so the
    column sums to total wall time. The rescaling preserves ranking but
    fixes the percentage interpretation. Raw isolated ms is also shown.
    """
    n_steps   = result['run']['n_steps']
    total_ms  = result['run']['t_run_total_ms']
    bridge_ms = result['run']['bridge_total_ms']
    cpp_ms    = result['run']['cm_update_total_ms']

    tf = result['tf_buckets']
    raw_per_step = {
        'TF capsule-capsule'      : tf['inter_capsule'],
        'TF primitive forces'     : tf['primitive'],
        'TF node-only (internal)' : tf['internal_forces'],
        'TF k_reg'                : tf['k_reg'],
        'TF drag'                 : tf['drag'],
        'TF time integration'     : tf['integration_only'],
    }
    raw_total_per_step = sum(raw_per_step.values())
    isolated_sum_ms    = raw_total_per_step * n_steps

    tf_envelope_ms = max(0.0, total_ms - cpp_ms - bridge_ms)
    # Rescale isolated TF times to the measured envelope.
    if isolated_sum_ms > 1e-9:
        scale = tf_envelope_ms / isolated_sum_ms
    else:
        scale = 0.0

    rows = []
    for k, v_per_step in raw_per_step.items():
        raw_ms     = v_per_step * n_steps
        scaled_ms  = raw_ms * scale
        rows.append((k, scaled_ms, scaled_ms / n_steps,
                     scaled_ms / total_ms * 100.0, raw_ms))

    rows.append(('C++ candidacy mgr', cpp_ms,
                 cpp_ms / n_steps, cpp_ms / total_ms * 100.0, cpp_ms))
    rows.append(('Python bridge', bridge_ms,
                 bridge_ms / n_steps, bridge_ms / total_ms * 100.0, bridge_ms))

    rows.sort(key=lambda r: r[1], reverse=True)
    return rows, total_ms, dict(scale=scale, isolated_sum_ms=isolated_sum_ms,
                                tf_envelope_ms=tf_envelope_ms)


def print_table(label, rows, total_ms, meta):
    print(f"\n── {label} bucket ranking " + "─" * 40)
    print(f"  {'bucket':<28s} {'ms':>10s} {'ms/step':>10s} {'%':>7s} {'raw ms':>10s}")
    print("  " + "-" * 72)
    for name, ms, ms_per_step, pct, raw in rows:
        print(f"  {name:<28s} {ms:10.1f} {ms_per_step:10.4f} {pct:6.1f}% {raw:10.1f}")
    print("  " + "-" * 72)
    print(f"  {'TOTAL':<28s} {total_ms:10.1f}")
    print(f"  TF envelope = total - cpp - bridge = {meta['tf_envelope_ms']:.1f} ms; "
          f"isolated TF sum = {meta['isolated_sum_ms']:.1f} ms; "
          f"rescale factor = {meta['scale']:.3f}")
    print(f"  (raw column = isolated per-call × n_steps; ms column = rescaled "
          f"to fit the envelope so percentages add to 100)")


def render_comparison(res_a, res_b):
    rows_a, tot_a, _ = render_bucket_table(res_a)
    rows_b, tot_b, _ = render_bucket_table(res_b)
    map_a = {r[0]: r for r in rows_a}
    map_b = {r[0]: r for r in rows_b}
    p_a   = res_a['P']
    p_b   = res_b['P']

    print(f"\n\n══════════ N={p_a} vs N={p_b} bucket comparison " + "═" * 18)
    print(f"  {'bucket':<28s} {f'{p_a} ms':>10s} {f'{p_b} ms':>10s} {'B/A':>10s}")
    print("  " + "-" * 64)
    all_keys = sorted(set(map_a) | set(map_b),
                      key=lambda k: -(map_a.get(k, (k,0,0,0))[1] + map_b.get(k, (k,0,0,0))[1]))
    for k in all_keys:
        ms_a = map_a.get(k, (k,0,0,0))[1]
        ms_b = map_b.get(k, (k,0,0,0))[1]
        ratio = ms_b / ms_a if ms_a > 0 else float('nan')
        print(f"  {k:<28s} {ms_a:10.1f} {ms_b:10.1f} {ratio:10.2f}x")
    print("  " + "-" * 64)
    print(f"  {'TOTAL':<28s} {tot_a:10.1f} {tot_b:10.1f} {tot_b/tot_a:10.2f}x")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n-steps',      type=int, default=2000)
    ap.add_argument('--sample-every', type=int, default=400)
    ap.add_argument('--p1', type=int, default=300)
    ap.add_argument('--p2', type=int, default=600)
    args = ap.parse_args()

    res1 = profile_one_size(P=args.p1, scale=1.0, label=f"N={args.p1}",
                            n_steps=args.n_steps, sample_every=args.sample_every)
    rows, tot, meta = render_bucket_table(res1)
    print_table(f"N={args.p1}", rows, tot, meta)
    res1['table'] = [{'bucket': r[0], 'ms': r[1], 'ms_per_step': r[2], 'pct': r[3], 'raw_ms': r[4]} for r in rows]

    scale2 = float(np.sqrt(args.p2 / args.p1))
    res2 = profile_one_size(P=args.p2, scale=scale2, label=f"N={args.p2}",
                            n_steps=args.n_steps, sample_every=args.sample_every)
    rows2, tot2, meta2 = render_bucket_table(res2)
    print_table(f"N={args.p2}", rows2, tot2, meta2)
    res2['table'] = [{'bucket': r[0], 'ms': r[1], 'ms_per_step': r[2], 'pct': r[3], 'raw_ms': r[4]} for r in rows2]

    render_comparison(res1, res2)

    out_p1 = OUT_DIR / f"hopper_N{args.p1}.json"
    out_p2 = OUT_DIR / f"hopper_N{args.p2}.json"
    with open(out_p1, 'w') as f:
        json.dump({k: v for k, v in res1.items() if not k.startswith('_')}, f, indent=2, default=str)
    with open(out_p2, 'w') as f:
        json.dump({k: v for k, v in res2.items() if not k.startswith('_')}, f, indent=2, default=str)
    print(f"\nSaved {out_p1} and {out_p2}")


if __name__ == "__main__":
    main()
