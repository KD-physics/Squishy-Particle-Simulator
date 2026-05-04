"""
initializer.py — RSA seeder and adaptive swell for the EPD user API.

Public API
----------
    rsa_seed(specs, objects, Lx, Ly, seed=42, max_attempts=50000)
        → particles, positions

    adaptive_swell(state, params, cm_mgr, prim_data,
                   phi_target, Lx, Ly, ...)
        → state, Lx, Ly

    compute_phi_outer(state, Lx, Ly, r_c_arr)
        → float  (shoelace packing fraction)
"""

import numpy as np
import tempfile, os

# ── exclusion geometry ────────────────────────────────────────────────────────

def _point_excluded(x, y, objects):
    """
    Return True if point (x, y) is in an excluded region defined by objects.
    An excluded region is determined by each object's .contains() method.
    A point is EXCLUDED if any container says it's outside or any obstacle says inside.
    """
    for obj in objects:
        if not obj.contains([x, y], t=0.0):
            return True
    return False


# ── RSA seeder ────────────────────────────────────────────────────────────────

def rsa_seed(specs, objects, Lx, Ly, seed=42, max_attempts=50000,
             verbose=True):
    """
    Place particles from specs via Random Sequential Adsorption.

    Parameters
    ----------
    specs        : list[ParticleSpec] — defines counts, N_nodes, poly_dist
    objects      : list[SimulationObject] — defines domain exclusions
    Lx, Ly       : float — box dimensions
    seed         : int   — RNG seed
    max_attempts : int   — max placement attempts per particle before giving up
    verbose      : bool

    Returns
    -------
    particles    : list[CapsuleParticle] in placement order
    positions    : (P_total, 2) ndarray — placed CM positions
    """
    from src.simulation.capsule_shell import CapsuleParticle

    rng = np.random.default_rng(seed)

    # Collect all particles from all specs (with sampled R0)
    all_R0_raw = []
    all_specs  = []
    for spec in specs:
        R0_arr = spec.sample_R0(seed=int(rng.integers(0, 2**32)))
        all_R0_raw.append(R0_arr)
        all_specs.extend([(spec, R0_arr[i]) for i in range(spec.count)])

    # Shuffle placement order for unbiased RSA
    order = rng.permutation(len(all_specs))

    # Compute effective radii for overlap check
    def _r_eff(spec, R0_i):
        N = spec.N_nodes
        r_c = 2.0 * R0_i * np.sin(np.pi / N)   # EPD capsule half-gap
        return R0_i + r_c                         # touching radius

    # Compute bounding box for random draws: use tightest 'exterior' container polygon
    bbox_x_lo, bbox_x_hi = 0.0, Lx
    bbox_y_lo, bbox_y_hi = 0.0, Ly
    for obj in objects:
        if getattr(obj, '_exclusion', None) == 'exterior':
            poly = obj.region_polygon(0.0)
            if poly is not None:
                v = np.asarray(poly['vertices'])
                bbox_x_lo = max(bbox_x_lo, float(v[:, 0].min()))
                bbox_x_hi = min(bbox_x_hi, float(v[:, 0].max()))
                bbox_y_lo = max(bbox_y_lo, float(v[:, 1].min()))
                bbox_y_hi = min(bbox_y_hi, float(v[:, 1].max()))

    P_total    = len(all_specs)
    placed_pos = np.zeros((P_total, 2))
    placed_reff = np.zeros(P_total)
    n_placed   = 0

    for idx in order:
        spec_i, R0_i = all_specs[idx]
        reff_i = _r_eff(spec_i, R0_i)

        placed = False
        for _ in range(max_attempts):
            x = rng.uniform(bbox_x_lo, bbox_x_hi)
            y = rng.uniform(bbox_y_lo, bbox_y_hi)

            # 1. Exclusion check
            if _point_excluded(x, y, objects):
                continue

            # 2. Object clearance check (CM must be at least reff_i from any wall)
            too_close_wall = False
            for obj in objects:
                for pdict in obj.resolved(t=0.0):
                    from src.simulation.contact_primitives import LineSegment, Arc
                    prim = pdict['prim']
                    if isinstance(prim, LineSegment):
                        # Distance from point to segment
                        d = _pt_to_seg_dist([x, y], prim.p0, prim.p1)
                        if d < reff_i:
                            too_close_wall = True
                            break
                    elif isinstance(prim, Arc):
                        cx, cy = prim.center
                        r = np.sqrt((x - cx)**2 + (y - cy)**2)
                        excl = pdict.get('exclusion')
                        if excl == 'interior':
                            # Convex obstacle: particle outside, clearance = r - R
                            if (r - prim.radius) < reff_i:
                                too_close_wall = True
                                break
                        else:
                            # Concave container (exclusion='exterior' or None):
                            # particle inside, clearance = R - r
                            if (prim.radius - r) < reff_i:
                                too_close_wall = True
                                break
                if too_close_wall:
                    break
            if too_close_wall:
                continue

            # 3. Overlap check with already-placed particles (periodic images)
            overlap = False
            for j in range(n_placed):
                dx = x - placed_pos[j, 0]
                dy = y - placed_pos[j, 1]
                # Minimum image for periodic BC
                dx -= Lx * np.round(dx / Lx)
                dy -= Ly * np.round(dy / Ly)
                d_cm = np.sqrt(dx**2 + dy**2)
                if d_cm < reff_i + placed_reff[j]:
                    overlap = True
                    break
            if overlap:
                continue

            # Accept
            placed_pos[n_placed]  = [x, y]
            placed_reff[n_placed] = reff_i
            n_placed += 1
            placed = True
            break

        if not placed:
            raise RuntimeError(
                f"RSA failed to place particle {idx} after {max_attempts} attempts. "
                f"Try lower phi_init or larger box.")

    if verbose:
        phi_rsa = np.sum(np.pi * placed_reff**2) / (Lx * Ly)
        print(f"  RSA: placed {P_total} particles, φ_eff={phi_rsa:.4f}, "
              f"Lx={Lx:.3f}, Ly={Ly:.3f}")

    # Build particle list in original spec order (not placement order)
    # Dispatch to EmulsionParticle for emulsion specs, CapsuleParticle otherwise.
    placement_map = {int(order[i]): i for i in range(P_total)}
    particles = []
    d_params = [spec_i.derived for spec_i, _ in all_specs]
    for idx in range(P_total):
        pi        = placement_map[idx]
        spec_i, R0_i = all_specs[idx]
        d         = d_params[idx]
        cx, cy    = float(placed_pos[pi, 0]), float(placed_pos[pi, 1])
        R0_i_f    = float(R0_i)
        if spec_i.type == 'emulsion':
            from src.simulation.emulsion_particle import EmulsionParticle
            gamma_eff = d['gamma'] * R0_i_f / spec_i.R0_mean
            p = EmulsionParticle(
                N=spec_i.N_nodes,
                R0=R0_i_f,
                gamma=gamma_eff,
                K_area=d['K_area'],   # constant — κ preserved
                C=d['C'],             # constant — C̃ preserved
                rho_d=1.0,
                center=(cx, cy),
            )
            p._kappa_target = 1.0 / spec_i._q_eff
            p._C_tilde      = d['C'] * spec_i.R0_mean / d['gamma']
            p._Oh_target    = spec_i.Oh
        else:
            p = CapsuleParticle(
                N=spec_i.N_nodes,
                R0=R0_i_f,
                tau=d['tau'],
                S=1.0,
                C=d['C'],
                K_area=d['K_area'],
                center=(cx, cy),
            )
            p._nu_target = spec_i.nu
            p._q_target  = spec_i._q_eff
            p._Oh_target = spec_i.Oh
        p.alpha_damp = spec_i._alpha_eff * spec_i.R0_mean / R0_i_f
        p.xi_drag    = spec_i._xi_drag   * R0_i_f / spec_i.R0_mean
        particles.append(p)

    # positions in same order as particles
    positions = np.array([placed_pos[placement_map[i]] for i in range(P_total)])
    return particles, positions


def _pt_to_seg_dist(pt, p0, p1):
    """Minimum distance from point pt to segment p0-p1."""
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    pt = np.asarray(pt, dtype=float)
    d  = p1 - p0
    t  = np.clip(np.dot(pt - p0, d) / max(np.dot(d, d), 1e-30), 0, 1)
    return float(np.linalg.norm(pt - (p0 + t * d)))


# ── packing fraction measurement ──────────────────────────────────────────────

def compute_phi_outer(state, Lx, Ly, r_c_arr, accessible_area=None):
    """
    Shoelace packing fraction: sum of outer-perimeter enclosed areas / accessible_area.

    r_c_arr        : (P,) — per-particle half-gap radius
    accessible_area: float or None — if provided, replaces Lx*Ly in denominator.
                     Use Lx*Ly - sum(obstacle.exclusion_area()) to exclude obstacles.
    """
    x_all = state['x_all'].numpy() if hasattr(state['x_all'], 'numpy') else np.asarray(state['x_all'])
    x_cm  = state['x_cm'].numpy()  if hasattr(state['x_cm'],  'numpy') else np.asarray(state['x_cm'])
    P, N, _ = x_all.shape
    total = 0.0
    for i in range(P):
        xy  = x_all[i]           # (N, 2)
        xcm = x_cm[i]            # (2,)
        # Outer perimeter: shift each node radially outward by r_c_arr[i]
        n_hat = xy - xcm[None, :]
        lens  = np.linalg.norm(n_hat, axis=1, keepdims=True)
        n_hat = n_hat / np.maximum(lens, 1e-12)
        outer = xy + float(r_c_arr[i]) * n_hat   # (N, 2)
        xs, ys = outer[:, 0], outer[:, 1]
        total += 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(np.roll(xs, -1), ys))
    denom = float(accessible_area) if accessible_area is not None else (Lx * Ly)
    return total / denom


# ── compression / wrapping ────────────────────────────────────────────────────

def _compress(state, Lx, Ly, phi_now, phi_new):
    """Affine box compression: scale CMs, shift nodes rigidly with CM."""
    import tensorflow as tf
    DTYPE = tf.float64
    scale    = np.sqrt(phi_now / phi_new)
    Lx_n     = Lx * scale
    Ly_n     = Ly * scale
    x_cm_np  = state['x_cm'].numpy() * scale
    x_all_np = state['x_all'].numpy() + (x_cm_np - state['x_cm'].numpy())[:, None, :]
    new_state = {**state,
                 'x_cm':  tf.constant(x_cm_np,  dtype=DTYPE),
                 'x_all': tf.constant(x_all_np, dtype=DTYPE)}
    return new_state, Lx_n, Ly_n


def _wrap(state, Lx, Ly):
    """Wrap particle CMs (and nodes rigidly) into [0,Lx) × [0,Ly)."""
    import tensorflow as tf
    DTYPE = tf.float64
    x_all = state['x_all'].numpy()
    x_cm  = state['x_cm'].numpy()
    dx = x_cm[:, 0] % Lx - x_cm[:, 0]
    dy = x_cm[:, 1] % Ly - x_cm[:, 1]
    mx = np.abs(dx) > 1e-9;  my = np.abs(dy) > 1e-9
    if mx.any():
        x_cm[mx, 0]    += dx[mx]
        x_all[mx, :, 0] += dx[mx, None]
    if my.any():
        x_cm[my, 1]    += dy[my]
        x_all[my, :, 1] += dy[my, None]
    return {**state,
            'x_all': tf.constant(x_all, dtype=DTYPE),
            'x_cm':  tf.constant(x_cm,  dtype=DTYPE)}


def _compute_max_f(state, cm_mgr, Lx, Ly, L0_arr):
    """Max fractional capsule overlap (dimensionless, 0=no overlap, 1=fully overlapping)."""
    x_all  = state['x_all'].numpy()
    P, N, _ = x_all.shape
    x_flat = x_all.reshape(P * N, 2)
    caps   = cm_mgr.CapCandidates   # (K, E) int32
    x_pad  = np.vstack([np.zeros((1, 2)), x_flat])
    E      = caps.shape[1]
    cand_flat = caps.ravel()
    xa = np.repeat(x_flat, E, axis=0)
    xb = x_pad[cand_flat]
    active = cand_flat > 0
    if not active.any():
        return 0.0
    dx = xa[:, 0] - xb[:, 0]
    dx -= Lx * np.round(dx / Lx)
    dy = xa[:, 1] - xb[:, 1]
    dy -= Ly * np.round(dy / Ly)
    dist = np.sqrt(dx**2 + dy**2)
    node_idx_a  = np.arange(P * N)
    L0_a_rep    = np.repeat(np.repeat(L0_arr, N), E)
    node_idx_b  = np.maximum(cand_flat - 1, 0)
    L0_b_rep    = np.repeat(L0_arr, N)[node_idx_b]
    L0_contact  = 0.5 * (L0_a_rep + L0_b_rep)
    f_vals = np.where(active, np.maximum(0.0, (2*L0_contact - dist) / (2*L0_contact)), 0.0)
    return float(f_vals.max())


def _compute_max_v(state):
    """Max perimeter node speed."""
    v_cm  = state['v_cm'].numpy()
    omega = state['omega'].numpy()
    u_dot = state['u_dot'].numpy()
    x_all = state['x_all'].numpy()
    x_cm  = state['x_cm'].numpy()
    x_rel = x_all - x_cm[:, None, :]
    v_rot = np.stack([-omega[:, None] * x_rel[:, :, 1],
                       omega[:, None] * x_rel[:, :, 0]], axis=-1)
    v_node = v_cm[:, None, :] + v_rot + u_dot
    return float(np.max(np.linalg.norm(v_node, axis=-1)))


def _compute_min_circularity(state):
    """Min 4πA/L² across all particles (1=circle, 0=degenerate)."""
    x_all = state['x_all'].numpy()
    P, N, _ = x_all.shape
    min_c = 1.0
    for i in range(P):
        xy = x_all[i]
        xs, ys = xy[:, 0], xy[:, 1]
        A = 0.5 * abs(np.dot(xs, np.roll(ys, -1)) - np.dot(np.roll(xs, -1), ys))
        edges = np.roll(xy, -1, axis=0) - xy
        L = np.sum(np.linalg.norm(edges, axis=1))
        c = 4.0 * np.pi * A / (L**2) if L > 1e-12 else 0.0
        if c < min_c:
            min_c = c
    return min_c


# ── adaptive swell ────────────────────────────────────────────────────────────

def adaptive_swell(state, params, cm_mgr, prim_data,
                   phi_target, Lx, Ly,
                   dphi_init=0.002, dphi_max=0.018, dphi_min=0.0001,
                   n_relax=1500, max_extra_relax=3000,
                   f_warn=0.45, f_crit=0.65,
                   circ_crit=0.30, circ_warn=0.50,
                   V_CRIT_FRAC=50.0, max_restores=6,
                   swell_alpha=10.0,
                   objects=None,
                   verbose=True,
                   skin=None, R0_max=None,
                   cand_check_interval=10,
                   use_tf_function=False):
    """
    Adaptively swell a particle packing from its current state to phi_target.

    Compression-relaxation loop:
    - Compress box by dphi
    - Relax for n_relax steps
    - If f > f_crit: restore checkpoint and halve dphi
    - If f > f_warn: extra relaxation until settled
    - If f < f_warn: advance phi, double dphi (up to dphi_max)
    - Stop when phi_outer >= phi_target

    Sets state['t'] = 0.0 on completion.

    Returns
    -------
    state, Lx, Ly  (updated)
    """
    import tensorflow as tf
    from src.simulation.tf_sim import (step_full_tf, run_simulation_tf,
                                        set_periodic_box, DTYPE, NP_DTYPE)

    # Extract per-particle arrays from params
    r_c_arr = params['r_c_per_p'].numpy() if hasattr(params['r_c_per_p'], 'numpy') \
              else np.asarray(params['r_c_per_p'])
    L0_arr  = params['L0'].numpy() if hasattr(params['L0'], 'numpy') \
              else np.asarray(params['L0'])

    # dt, alpha, g scalars — use swell_alpha (higher damping kills ringing faster)
    from src.simulation.capsule_shell import CapsuleParticle
    dt_tf    = params.get('_dt_tf', tf.constant(NP_DTYPE(0.01)))
    dt_val   = float(dt_tf.numpy())
    alpha_tf = tf.constant(NP_DTYPE(swell_alpha))
    g_tf     = tf.constant(NP_DTYPE(0.0))

    # Velocity scale for v_crit
    v_scale  = float(np.mean(r_c_arr)) / dt_val
    v_crit   = V_CRIT_FRAC * v_scale

    # Skin / R0_max for run_simulation_tf's candidacy threshold check.
    # If caller didn't pass them, fall back to safe defaults from cm_mgr
    # (mean-based; for polydisperse the worst case is just a slightly looser
    # check, so candidacy may refresh more often than strictly necessary).
    if skin is None:
        skin = float(getattr(cm_mgr, 'skin', np.mean(r_c_arr)))
    if R0_max is None:
        R0_max = float(getattr(cm_mgr, 'R0', 1.0))

    # Mutable prim_data container so object rescaling is seen by _relax_block
    _pd = [prim_data]

    def _relax_block(st, n_steps):
        """Run n_steps of relax via run_simulation_tf (one tf.while_loop call).
        Returns the wrapped state.

        NOTE: use_tf_function is forced to False inside the swell regardless of
        the caller's preference. The current wrap rebuilds a tf.function on
        every call, so adaptive_swell's many short relax_blocks would trigger a
        5–10s retrace each — net negative for swell. The wrap is only a win for
        long single-chunk calls (sys.run with large n_steps). Followup task:
        cache the compiled runner so the wrap composes with repeated short
        calls; until then, swell uses cand_check_interval as the only lever.

        For the swell, alpha is forced to swell_alpha by overriding the TF
        scalar inside `params` for the duration of the call, then restored.
        Same trick for g (swell uses g=0). dt is unchanged.
        """
        if n_steps <= 0:
            return st
        phys = {k: v for k, v in st.items() if k not in ('t', 'step')}
        # Inject swell-specific damping/g into params, remember originals.
        # Critical: kernel uses tf.where(per_p > 0, per_p, scalar) so we MUST
        # override the per-particle alpha too — otherwise per_p (~R0_mean ≈ 1)
        # wins and contacts ring under-damped during the swell.
        _alpha_save    = params.get('_alpha_tf')
        _alpha_pp_save = params.get('alpha_damp_per_p')
        _g_save        = params.get('_g_tf')
        params['_alpha_tf'] = alpha_tf
        params['_g_tf']     = g_tf
        if _alpha_pp_save is not None:
            P_count = int(_alpha_pp_save.shape[0])
            params['alpha_damp_per_p'] = tf.constant(
                np.full(P_count, swell_alpha, dtype=np.float64),
                dtype=_alpha_pp_save.dtype)
        try:
            new_phys = run_simulation_tf(
                phys, dt_val, swell_alpha, 0.0,
                params, n_steps, cm_mgr,
                skin=skin, prim_data=_pd[0], R0_max=R0_max,
                cand_check_interval=cand_check_interval,
                diagnostics=False,
                step_offset=0,
                use_tf_function=False,    # forced — see docstring above
            )
        finally:
            if _alpha_save is not None:
                params['_alpha_tf'] = _alpha_save
            if _alpha_pp_save is not None:
                params['alpha_damp_per_p'] = _alpha_pp_save
            if _g_save is not None:
                params['_g_tf'] = _g_save
        # Re-attach t/step if the original state had them (they're not stepped
        # by run_simulation_tf — adaptive_swell resets t at the end anyway).
        if 't' in st:
            new_phys['t'] = st['t']
        if 'step' in st:
            new_phys['step'] = st['step']
        return _wrap(new_phys, Lx, Ly)

    def _step(st):
        """Single-step variant retained for callers that need fine-grained
        Python control between steps (currently unused after the relax-loop
        rewrite, but kept for back-compat / future use)."""
        if cm_mgr.needs_update(st['x_cm'].numpy(), st['theta'].numpy()):
            cm_mgr.update(st['x_cm'].numpy(), st['theta'].numpy(),
                           x_all=st['x_all'].numpy())
        caps    = tf.constant(cm_mgr.CapCandidates, dtype=tf.int32)
        st, _   = step_full_tf(st, caps, dt_tf, alpha_tf, g_tf,
                                params, t=tf.constant(NP_DTYPE(0.0)),
                                prim_data=_pd[0])
        return _wrap(st, Lx, Ly)   # note: Lx/Ly updated via closure below

    def _rescale_objects(scale):
        """Scale all registered objects and rebuild prim_data."""
        if objects is None:
            return
        from src.simulation.tf_sim import make_prim_data as _make_prim_data
        for obj in objects:
            obj.rescale(scale)
        prim_list = []
        for obj in objects:
            prim_list.extend(obj.to_make_prim_list(t=0.0))
        _pd[0] = _make_prim_data(prim_list)

    # Checkpoint (in-memory)
    ckpt = {'state': state, 'Lx': Lx, 'Ly': Ly}

    def save_ckpt():
        ckpt['state'] = {k: tf.identity(v) for k, v in state.items()}
        ckpt['Lx'] = Lx; ckpt['Ly'] = Ly

    def load_ckpt():
        return (
            {k: tf.identity(v) for k, v in ckpt['state'].items()},
            float(ckpt['Lx']),
            float(ckpt['Ly']),
        )

    # Initial relax
    n_init = min(n_relax, 3000)
    if verbose:
        print(f"  Initial relax: {n_init} steps "
              f"(cand_check_interval={cand_check_interval}) ...")
    state = _relax_block(state, n_init)

    # Update cm_mgr with current box
    cm_mgr.Lx = Lx; cm_mgr.Ly = Ly

    phi_now   = compute_phi_outer(state, Lx, Ly, r_c_arr)
    phi_act   = phi_now
    dphi      = dphi_init
    total_steps = n_init
    consecutive_restores = 0

    if verbose:
        print(f"  Start swell: phi_outer={phi_now:.4f} → target={phi_target:.4f} "
              f"dphi_max={dphi_max:.4f}")

    while phi_act < phi_target - 1e-6:
        save_ckpt()

        phi_new = phi_now + dphi
        scale   = np.sqrt(phi_now / phi_new)   # compression factor (< 1)
        state, Lx, Ly = _compress(state, Lx, Ly, phi_now, phi_new)
        state = _wrap(state, Lx, Ly)
        cm_mgr.Lx = Lx; cm_mgr.Ly = Ly
        set_periodic_box(params, Lx, Ly)
        # Compress is an affine rescale of every CM — invalidates L1/L2 ranking
        # AND Q (n0, m0) basis. Force a full PRCM rebuild from scratch and
        # reset cadence; subsequent in-relax PRCM calls then start clean.
        # Use force_full_rebuild on PRCM if available, else fall back to the
        # cadenced update (production CandidacyManager).
        x_cm_np  = state['x_cm'].numpy()
        x_all_np = state['x_all'].numpy()
        theta_np = state['theta'].numpy()
        if hasattr(cm_mgr, '_cpp') and hasattr(cm_mgr._cpp, 'force_full_rebuild'):
            cm_mgr._cpp.force_full_rebuild(x_all_np, x_cm_np, theta_np)
            # Sync the Python-side CapCandidates buffer to the freshly-built C
            cm_mgr.CapCandidates[:] = cm_mgr._cpp.CapCandidates
        else:
            cm_mgr.update(x_cm_np, theta_np, x_all=x_all_np)
        _rescale_objects(scale)

        # Quick assessment after compression
        f_immed = _compute_max_f(state, cm_mgr, Lx, Ly, L0_arr)
        v_immed = _compute_max_v(state)

        if f_immed >= f_crit or v_immed >= v_crit:
            reason = 'f_crit' if f_immed >= f_crit else 'v_crit'
            consecutive_restores += 1
            state, Lx, Ly = load_ckpt()
            cm_mgr.Lx = Lx; cm_mgr.Ly = Ly
            set_periodic_box(params, Lx, Ly)
            dphi = max(dphi * 0.5, dphi_min)
            if verbose:
                print(f"    restore ({reason}) → dphi={dphi:.5f}  consec={consecutive_restores}")
            if consecutive_restores >= max_restores:
                print(f"  *** Stuck at phi={phi_now:.4f} after {max_restores} restores ***")
                break
            continue

        consecutive_restores = 0

        # Relax
        state = _relax_block(state, n_relax)
        total_steps += n_relax

        # Post-relax metrics
        f_post    = _compute_max_f(state, cm_mgr, Lx, Ly, L0_arr)
        v_post    = _compute_max_v(state)
        circ_post = _compute_min_circularity(state)

        # Circularity collapse → treat as critical, restore immediately
        if circ_post <= circ_crit:
            consecutive_restores += 1
            if verbose:
                print(f"    restore (circ={circ_post:.3f}≤{circ_crit}) "
                      f"→ dphi={dphi*0.25:.5f}  consec={consecutive_restores}")
            state, Lx, Ly = load_ckpt()
            cm_mgr.Lx = Lx; cm_mgr.Ly = Ly
            set_periodic_box(params, Lx, Ly)
            dphi = max(dphi * 0.25, dphi_min)   # more aggressive backoff for shape collapse
            if consecutive_restores >= max_restores:
                print(f"  *** Stuck (circ) at phi={phi_now:.4f} after {max_restores} restores ***")
                break
            continue

        # Extra relax if still overlapping — chunk into 200-step blocks
        # so we can still poll metrics between chunks (same cadence as before)
        if f_post >= f_warn:
            extra = 0
            while f_post >= f_warn and extra < max_extra_relax:
                step_chunk = min(200, max_extra_relax - extra)
                state = _relax_block(state, step_chunk)
                extra += step_chunk
                f_post    = _compute_max_f(state, cm_mgr, Lx, Ly, L0_arr)
                circ_post = _compute_min_circularity(state)
                if circ_post <= circ_crit:
                    break   # shape collapsed during extra relax → exit and restore next iter
            total_steps += extra
            f_post    = _compute_max_f(state, cm_mgr, Lx, Ly, L0_arr)
            circ_post = _compute_min_circularity(state)

        phi_now = phi_new
        phi_act = compute_phi_outer(state, Lx, Ly, r_c_arr)

        # Rate adaptation
        if f_post < f_warn * 0.5 and circ_post > circ_warn:
            dphi = min(dphi * 1.2, dphi_max)
        elif circ_post <= circ_warn:
            dphi = max(dphi * 0.7, dphi_min)   # slow down near shape limit

        if verbose and total_steps % (n_relax * 5) < n_relax:
            print(f"    phi_box={phi_now:.4f}  phi_outer={phi_act:.4f}  "
                  f"f={f_post:.3f}  circ_min={circ_post:.3f}  steps={total_steps}")

        if phi_act >= phi_target - 1e-6:
            break

    # Reset time
    state = {**state, 't': tf.constant(NP_DTYPE(0.0))} if 't' in state else state

    if verbose:
        phi_final = compute_phi_outer(state, Lx, Ly, r_c_arr)
        print(f"  Swell complete: phi_outer={phi_final:.4f}  steps={total_steps}")

    return state, Lx, Ly
