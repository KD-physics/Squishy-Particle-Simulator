"""
tf_sim.py — TensorFlow-accelerated EPD/Emulsion simulation engine.

Architecture
------------
  C++ CandidacyManager  →  CapCandidates (K×E int32)  →  TF step graph
       ↑                                                        │
       └──────── x_cm (P,2),  theta (P,)  ────────────────────┘

All force and integration arithmetic lives in one @tf.function per step.
The C++ manager owns neighbor topology; Python orchestrates.

Data structures (fixed shape after init — required for static TF graph)
-----------------------------------------------------------------------
  P   — number of particles (uniform N nodes each)
  N   — nodes per particle
  E   — max candidates per edge (skin-buffered)
  K   = P * N   — total edges; global index k = p*N + e, k∈[1..K]; 0 = ghost

  CapCandidates : (K, E)    int32  — candidate edge indices; 0 = ghost
  x_all         : (P, N, 2) DTYPE
  x_cm          : (P, 2)    DTYPE
  v_cm          : (P, 2)    DTYPE
  theta         : (P,)      DTYPE
  omega         : (P,)      DTYPE
  u             : (P, N, 2) DTYPE  — elastic displacement, body frame
  u_dot         : (P, N, 2) DTYPE

DTYPE: tf.float32 by default; set DTYPE = tf.float64 for full precision.
C++ buffers use the matching C type — no silent truncation.
"""

import numpy as np
import tensorflow as tf

# ── dtype flag ────────────────────────────────────────────────────────────────
DTYPE    = tf.float32
NP_DTYPE = np.float32

def set_dtype(dtype):
    """Switch global DTYPE (call before constructing any TFSim)."""
    global DTYPE, NP_DTYPE
    DTYPE    = dtype
    NP_DTYPE = np.float64 if dtype == tf.float64 else np.float32


# ── geometry helpers (batched over P) ─────────────────────────────────────────

def _shoelace_area(x):
    """Signed shoelace area.  x: (P, N, 2) → (P,)"""
    xr = tf.roll(x, shift=-1, axis=1)
    return 0.5 * tf.reduce_sum(
        x[:, :, 0] * xr[:, :, 1] - xr[:, :, 0] * x[:, :, 1], axis=1)


def _edge_tangents_normals(x):
    """
    Returns t_hat (P,N,2), n_hat (P,N,2), L (P,N).
    n_hat is the 90°-CCW rotation of t_hat (INWARD for CCW polygon).
    """
    xnext = tf.roll(x, shift=-1, axis=1)
    e     = xnext - x                              # (P, N, 2)
    L     = tf.norm(e, axis=2)                    # (P, N)
    t_hat = e / tf.maximum(L[:, :, None], 1e-30)
    n_hat = tf.stack([-t_hat[:, :, 1], t_hat[:, :, 0]], axis=2)
    return t_hat, n_hat, L


def _node_normals(n_hat):
    """Average of adjacent edge normals, normalised.  (P,N,2)→(P,N,2)"""
    n_prev = tf.roll(n_hat, shift=1, axis=1)
    avg    = n_hat + n_prev
    norms  = tf.maximum(tf.norm(avg, axis=2, keepdims=True), 1e-30)
    return avg / norms


# ── internal forces kernel ────────────────────────────────────────────────────

@tf.function(jit_compile=True)
def internal_forces_tf(x, params, g):
    """
    Internal forces per particle (pressure + edge + bending + optional gravity).

    x      : (P, N, 2)   DTYPE — current node positions
    params : dict of (P,) DTYPE tensors: L0, A0, K_area, El_t, EI, rho_d, theta0
    g      : scalar DTYPE

    Returns f : (P, N, 2) DTYPE
    """
    t_hat, n_hat, L = _edge_tangents_normals(x)
    n_node          = _node_normals(n_hat)

    # Scalar (P,) params — expand to (P,1,1) for (P,N,2) broadcast or (P,1) for (P,N)
    L0_s     = params['L0']         # (P,)
    A0_s     = params['A0']
    K_area_s = params['K_area']
    El_t_s   = params['El_t']
    EI_s     = params['EI']

    # ── 1. Fluid pressure ─────────────────────────────────────────────────────
    A   = tf.abs(_shoelace_area(x))                       # (P,)
    P_s = K_area_s * (A0_s - A) / A0_s                   # (P,) positive when compressed
    # n_node is INWARD; subtract to push outward
    Fp  = (P_s * L0_s)[:, None, None] * n_node           # (P,1,1)*(P,N,2) = (P,N,2)
    f   = -Fp

    # ── 2. Hydrostatic gravity ────────────────────────────────────────────────
    if g != 0.0:
        rho_d  = params['rho_d']                          # (P,)
        x_next = tf.roll(x, shift=-1, axis=1)
        y_top  = tf.reduce_max(x[:, :, 1], axis=1)       # (P,)
        y_mid  = 0.5 * (x[:, :, 1] + x_next[:, :, 1])   # (P, N)
        P_mid  = rho_d[:, None] * g * (y_top[:, None] - y_mid)   # (P, N)
        # outward normal = -n_hat; distribute half to each endpoint
        dF_edge = P_mid[:, :, None] * L[:, :, None] * (-n_hat)   # (P, N, 2)
        # Enforce Σ Fx = 0 per particle
        fx_mean = tf.reduce_mean(dF_edge[:, :, 0], axis=1, keepdims=True)  # (P,1)
        dF_x    = dF_edge[:, :, 0] - fx_mean             # (P, N)
        dF_edge = tf.stack([dF_x, dF_edge[:, :, 1]], axis=2)      # (P, N, 2)
        f = f + 0.5 * dF_edge + 0.5 * tf.roll(dF_edge, shift=1, axis=1)

    # ── 3. Edge elasticity ────────────────────────────────────────────────────
    eps    = (L - L0_s[:, None]) / L0_s[:, None]         # (P, N)
    F_edge = El_t_s[:, None, None] * eps[:, :, None] * t_hat   # (P, N, 2)
    f      = f + F_edge - tf.roll(F_edge, shift=1, axis=1)

    # ── 4. Bending (3-node hinge) ─────────────────────────────────────────────
    t_prev = tf.roll(t_hat, shift=1, axis=1)
    n_prev = tf.roll(n_hat, shift=1, axis=1)

    cross = t_prev[:, :, 0]*t_hat[:, :, 1] - t_prev[:, :, 1]*t_hat[:, :, 0]
    dot_  = t_prev[:, :, 0]*t_hat[:, :, 0] + t_prev[:, :, 1]*t_hat[:, :, 1]
    theta_cur = tf.math.atan2(cross, dot_)                # (P, N)

    M  = (EI_s[:, None] / L0_s[:, None]) * (theta_cur - params['theta0'])  # (P,N)
    Mb = M / L0_s[:, None]                                # (P, N)

    dF_im1 = -Mb[:, :, None] * n_prev
    dF_i   =  Mb[:, :, None] * (n_prev + n_hat)
    dF_ip1 = -Mb[:, :, None] * n_hat

    f = (f
         + tf.roll(dF_im1, shift=-1, axis=1)
         + dF_i
         + tf.roll(dF_ip1, shift=1, axis=1))

    # ── 5. Line tension (emulsion droplets; zero for elastic capsules) ─────────
    gamma_lt_s = params['gamma_lt']                      # (P,)
    t_hat_prev = tf.roll(t_hat, shift=1, axis=1)         # (P, N, 2)
    F_lt = gamma_lt_s[:, None, None] * (t_hat - t_hat_prev)   # net tangential tension
    f = f + F_lt

    return f


# ── rigid-body decomposition + integration ────────────────────────────────────

@tf.function(jit_compile=True)
def k_reg_forces_tf(x, params):
    """Tangential regularization pseudo-force that keeps nodes evenly spaced.
    Net sum != 0 for deformed shapes, so it must be excluded from RB sums."""
    k_reg   = params['k_reg_per_p']                        # (P,)
    t_hat, _, L = _edge_tangents_normals(x)                # (P,N,2), _, (P,N)
    L_plus  = L                                            # edge k:   x[k+1]-x[k]
    L_minus = tf.roll(L, shift=1, axis=1)                  # edge k-1: x[k]-x[k-1]
    t_prev  = tf.roll(t_hat, shift=1, axis=1)              # t̂ of edge k-1
    t_node  = 0.5 * (t_hat + t_prev)                       # bisector tangent at node k
    t_norm  = tf.maximum(tf.norm(t_node, axis=2, keepdims=True), 1e-30)
    t_node  = t_node / t_norm
    return k_reg[:, None, None] * 0.5 * (L_plus - L_minus)[:, :, None] * t_node


# ── background-flow presets (TF-native, runs inside tf.while_loop) ────────────

def _eval_U_bg(x_nodes, t, params):
    """
    Evaluate background flow U_bg at every node position.

    x_nodes : (P, N, 2) — node positions
    t       : scalar DTYPE — current time (unused by current presets)
    params  : must contain 'U_bg_type' (scalar int32) and 'U_bg_params' (4,) float64

    Returns U : (P, N, 2) — background velocity at each node

    Preset index:
        0 — zero (quiescent)
        1 — constant:       params[0]=Ux, params[1]=Uy
        2 — shear:          U = (rate*y, 0)       params[0]=rate
        3 — parabolic:      U = (U_max*(1-(y/H)^2), 0)  params[0]=U_max, params[1]=H
        4 — extensional:    U = (rate*x, -rate*y) params[0]=rate
        5 — poiseuille_v:   U = (0, -U_max*(1-((x-x_c)/H)^2))
                            params[0]=U_max, params[1]=H, params[2]=x_c
    """
    bg_type   = params.get('U_bg_type',   tf.constant(0, dtype=tf.int32))
    bg_params = params.get('U_bg_params', tf.zeros([4], dtype=DTYPE))

    y = x_nodes[:, :, 1]   # (P, N)
    x = x_nodes[:, :, 0]   # (P, N)
    zeros = tf.zeros_like(y)

    def _zero():
        return tf.zeros_like(x_nodes)

    def _constant():
        Ux = tf.fill(tf.shape(y), bg_params[0])
        Uy = tf.fill(tf.shape(y), bg_params[1])
        return tf.stack([Ux, Uy], axis=2)

    def _shear():
        rate = bg_params[0]
        return tf.stack([rate * y, zeros], axis=2)

    def _parabolic():
        U_max = bg_params[0]
        H     = bg_params[1]
        Ux    = U_max * (1.0 - (y / H) ** 2)
        return tf.stack([Ux, zeros], axis=2)

    def _extensional():
        rate = bg_params[0]
        return tf.stack([rate * x, -rate * y], axis=2)

    def _poiseuille_v():
        U_max = bg_params[0]
        H     = bg_params[1]
        x_c   = bg_params[2]
        Uy    = -U_max * (1.0 - ((x - x_c) / H) ** 2)
        return tf.stack([zeros, Uy], axis=2)

    return tf.switch_case(bg_type, branch_fns=[_zero, _constant, _shear,
                                               _parabolic, _extensional,
                                               _poiseuille_v])


def step_rb_tf(state, f_contact, dt, alpha_damp, g, params, t=None):
    """
    One rigid-body decomposition step.

    state       : dict with x_all, x_cm, v_cm, theta, omega, u, u_dot, X_ref
    f_contact   : (P, N, 2) DTYPE — nodal contact forces
    dt, alpha_damp, g : scalar DTYPE
    t           : scalar DTYPE — current time (for driven-particle trajectories)

    params may include (all default to zero/False if absent):
      driven_mask  (P,)    — 1 = CM velocity prescribed, 0 = force-integrated
      shape_frozen (P,)    — 1 = nodes locked to rigid-body frame (elastic DOFs frozen)
      traj         (P, 18) — DC+AC trajectory:
                             [vx_dc, vx_ac, vx_freq, vx_phase,
                              vy_dc, vy_ac, vy_freq, vy_phase,
                              ws_dc, ws_ac, ws_freq, ws_phase,   (spin)
                              wo_dc, wo_ac, wo_freq, wo_phase,   (orbit)
                              rx_ref, ry_ref]

    Returns new_state dict, metrics dict.
    """
    x     = state['x_all']    # (P, N, 2)
    x_cm  = state['x_cm']     # (P, 2)
    v_cm  = state['v_cm']
    theta = state['theta']     # (P,)
    omega = state['omega']
    u     = state['u']         # (P, N, 2)
    u_dot = state['u_dot']
    X_ref = state['X_ref']

    M_disk = params['M_disk']  # (P,)
    I_disk = params['I_disk']
    m_node = params['m_node']

    # ── Stokes drag (node-level, arc-length weighted) ─────────────────────────
    xi_drag = params.get('xi_drag_per_p',    tf.zeros_like(M_disk))  # (P,)
    # Per-particle alpha: use params entry when available; fall back to scalar arg.
    _alpha_pp = params.get('alpha_damp_per_p', tf.zeros_like(M_disk))  # (P,)
    _alpha_any = tf.reduce_any(_alpha_pp > 0.0)
    alpha_p = tf.where(_alpha_any,
                       _alpha_pp,
                       tf.fill(tf.shape(M_disk), alpha_damp))  # (P,)
    xi_any  = tf.reduce_any(xi_drag > 0.0)

    def _compute_drag():
        # Node velocity in lab frame: v_cm + omega×r + R(theta)·u_dot
        r      = x - x_cm[:, None, :]                             # (P, N, 2)
        v_rot  = omega[:, None, None] * tf.stack(
                     [-r[:, :, 1], r[:, :, 0]], axis=2)           # (P, N, 2)
        c_th   = tf.cos(theta)[:, None]
        s_th   = tf.sin(theta)[:, None]
        v_el_x = c_th * u_dot[:, :, 0] - s_th * u_dot[:, :, 1]
        v_el_y = s_th * u_dot[:, :, 0] + c_th * u_dot[:, :, 1]
        v_el   = tf.stack([v_el_x, v_el_y], axis=2)               # (P, N, 2)
        v_node = v_cm[:, None, :] + v_rot + v_el                   # (P, N, 2)

        # Background flow at each node
        t_val  = tf.constant(0.0, dtype=DTYPE) if t is None else t
        U_bg   = _eval_U_bg(x, t_val, params)                      # (P, N, 2)

        # Arc-length weights: dL_i = |x_{i+1} - x_{i-1}| / 2
        x_next = tf.roll(x, shift=-1, axis=1)
        x_prev = tf.roll(x,  shift=1, axis=1)
        dL     = 0.5 * tf.norm(x_next - x_prev, axis=2)           # (P, N)

        # f_drag,i = -xi * (v_node - U_bg) * dL_i
        return -xi_drag[:, None, None] * (v_node - U_bg) * dL[:, :, None]

    f_drag  = tf.cond(xi_any, _compute_drag, lambda: tf.zeros_like(x))

    # ── elastic + regularization forces ──────────────────────────────────────
    f_el    = internal_forces_tf(x, params, g)     # (P, N, 2)
    f_reg   = k_reg_forces_tf(x, params)           # (P, N, 2) — pseudo-force, excluded from RB
    f_total = f_el + f_reg + f_contact + f_drag

    # ── rigid-body resultants from contact forces only ────────────────────────
    # Subtract both elastic and regularization sums (neither drives RB motion)
    F_el    = tf.reduce_sum(f_el,    axis=1)        # (P, 2)
    F_reg   = tf.reduce_sum(f_reg,   axis=1)        # (P, 2)
    F_total = tf.reduce_sum(f_total, axis=1)
    F       = F_total - F_el - F_reg                # (P, 2) = Σf_contact only

    # Re-add gravity as rigid-body body force
    if g != 0.0:
        F = F + tf.stack([tf.zeros_like(M_disk), -M_disk * g], axis=1)

    r       = x - x_cm[:, None, :]                                      # (P, N, 2)
    T_total = tf.reduce_sum(r[:,:,0]*f_total[:,:,1] - r[:,:,1]*f_total[:,:,0], axis=1)
    T_el    = tf.reduce_sum(r[:,:,0]*f_el[:,:,1]    - r[:,:,1]*f_el[:,:,0],    axis=1)
    T_reg   = tf.reduce_sum(r[:,:,0]*f_reg[:,:,1]   - r[:,:,1]*f_reg[:,:,0],   axis=1)
    T       = T_total - T_el - T_reg                                     # (P,)

    # ── rigid-body integration with optional dissipation ─────────────────────
    # Semi-implicit Euler: v_new = (v_old + F/M·dt) / (1 + β·dt)
    # β has units 1/time and damps the rigid-body modes (v_cm and omega).
    # Distinct from xi_drag (per-node Stokes drag, physical). beta_rb is a
    # numerical knob for ad-hoc energy vacuum during initialization / settle.
    # Default zero → no effect (backwards-compatible Newtonian integration).
    # Per-particle override: beta_rb_per_p; scalar fallback: _beta_rb_tf.
    v_cm_free  = v_cm  + (F / M_disk[:, None]) * dt
    omega_free = omega + (T / I_disk) * dt
    _beta_pp   = params.get('beta_rb_per_p', tf.zeros_like(M_disk))
    _beta_any  = tf.reduce_any(_beta_pp > 0.0)
    _beta_scal = params.get('_beta_rb_tf', tf.constant(0.0, dtype=DTYPE))
    beta_p     = tf.where(_beta_any,
                          _beta_pp,
                          tf.fill(tf.shape(M_disk), _beta_scal))     # (P,)
    inv_factor = 1.0 / (1.0 + beta_p * dt)                            # (P,)
    v_cm_free  = v_cm_free  * inv_factor[:, None]
    omega_free = omega_free * inv_factor

    # ── driven-particle override ──────────────────────────────────────────────
    if 'driven_mask' in params and 'traj' in params:
        traj = params['traj']           # (P, 18)
        t_val = tf.constant(0.0, dtype=DTYPE) if t is None else t
        # Evaluate DC+AC trajectory for all particles
        vx  = traj[:,0] + traj[:,1] * tf.cos(traj[:,2]  * t_val + traj[:,3])
        vy  = traj[:,4] + traj[:,5] * tf.cos(traj[:,6]  * t_val + traj[:,7])
        ws  = traj[:,8] + traj[:,9] * tf.cos(traj[:,10] * t_val + traj[:,11])
        wo  = traj[:,12] + traj[:,13] * tf.cos(traj[:,14] * t_val + traj[:,15])
        r_ref = traj[:, 16:18]          # (P, 2)
        # Orbital contribution: v_orbit = ω_orbit * ⊥(x_cm - r_ref)
        dx = x_cm - r_ref               # (P, 2)
        v_orb = wo[:, None] * tf.stack([-dx[:, 1], dx[:, 0]], axis=1)
        v_prescribed  = tf.stack([vx, vy], axis=1) + v_orb  # (P, 2)
        omega_prescribed = ws                                 # (P,)
        mask   = params['driven_mask'][:, None]  # (P, 1)
        mask_s = params['driven_mask']           # (P,)
        v_cm_new  = v_cm_free  * (1.0 - mask)   + v_prescribed * mask
        omega_new = omega_free * (1.0 - mask_s)  + omega_prescribed * mask_s
    else:
        v_cm_new  = v_cm_free
        omega_new = omega_free

    x_cm_new  = x_cm  + v_cm_new * dt
    theta_new = theta + omega_new * dt

    # ── deviatoric force → elastic deformation ────────────────────────────────
    # Frozen branch: f_dev = f_phys - F_phys/N  where F_phys = F_total (sum of ALL forces)
    N_f   = tf.cast(tf.shape(x)[1], DTYPE)
    f_dev = f_total - F_total[:, None, :] / N_f                         # (P, N, 2)

    # Rotate f_dev to body frame: R(-θ) = [[c, s], [-s, c]]
    c = tf.cos(theta_new)[:, None]
    s = tf.sin(theta_new)[:, None]
    f_bx  =  c * f_dev[:,:,0] + s * f_dev[:,:,1]
    f_by  = -s * f_dev[:,:,0] + c * f_dev[:,:,1]
    f_body = tf.stack([f_bx, f_by], axis=2)                             # (P, N, 2)

    accel_el  = f_body / m_node[:, None, None] - alpha_p[:, None, None] * u_dot
    u_dot_el  = u_dot + accel_el * dt
    u_el      = u     + u_dot_el * dt

    # Shape-freeze: lock elastic DOFs to frozen reference shape (not necessarily zero)
    if 'shape_frozen' in params:
        frz   = params['shape_frozen'][:, None, None]   # (P, 1, 1)
        u_ref = params.get('u_frozen', tf.zeros_like(u_el))  # (P, N, 2)
        u_new     = u_el * (1.0 - frz) + u_ref * frz
        u_dot_new = u_dot_el * (1.0 - frz)
    else:
        u_new     = u_el
        u_dot_new = u_dot_el

    # ── perimeter reconstruct ─────────────────────────────────────────────────
    body  = X_ref + u_new                                                # (P, N, 2)
    c2    = tf.cos(theta_new)[:, None]
    s2    = tf.sin(theta_new)[:, None]
    x_new = tf.stack([
        c2 * body[:,:,0] - s2 * body[:,:,1],
        s2 * body[:,:,0] + c2 * body[:,:,1],
    ], axis=2) + x_cm_new[:, None, :]

    # ── metrics ───────────────────────────────────────────────────────────────
    A_new  = tf.abs(_shoelace_area(x_new))
    dA_rel = tf.abs(A_new - params['A0']) / params['A0']
    KE_rb  = (0.5 * M_disk * tf.reduce_sum(v_cm_new**2, axis=1)
             + 0.5 * I_disk * omega_new**2)
    KE_el  = 0.5 * m_node * tf.reduce_sum(
                 tf.reduce_sum(u_dot_new**2, axis=2), axis=1)

    new_state = dict(
        x_all = x_new,
        x_cm  = x_cm_new,
        v_cm  = v_cm_new,
        theta = theta_new,
        omega = omega_new,
        u     = u_new,
        u_dot = u_dot_new,
        X_ref = X_ref,
    )
    metrics = dict(dA_rel=dA_rel, KE_rb=KE_rb, KE_el=KE_el)
    return new_state, metrics


# ── inter-capsule contact forces (TF, one-directional, polydisperse) ─────────

@tf.function(jit_compile=True)
def inter_capsule_forces_tf(x_all, CapCandidates, r_c_flat, k_c_flat, L0_flat,
                             box=None):
    """
    TF inter-capsule contact forces from CapCandidates index matrix.

    One-directional (pA < pB rows only).  Source Gauss points on A;
    Newton 3rd scatter to candidate (B) nodes via tensor_scatter_nd_add.

    x_all         : (P, N, 2) DTYPE
    CapCandidates : (K, E)    int32  — 0 = ghost
    r_c_flat      : (K,) DTYPE — capsule radius per source edge
    k_c_flat      : (K,) DTYPE — contact stiffness per source edge
    L0_flat       : (K,) DTYPE — edge reference length per source edge
    box           : (2,) DTYPE [Lx, Ly] or None — periodic box; None = non-periodic

    Returns f_contact : (P, N, 2) DTYPE
    """
    dtype = x_all.dtype
    P = x_all.shape[0]
    N = x_all.shape[1]
    K = P * N
    E = CapCandidates.shape[1]

    # Ghost-prepended node arrays (index 0 = ghost at infinity)
    ghost  = tf.fill([1, 2], tf.constant(1e9, dtype=dtype))
    x_flat      = tf.reshape(x_all, [K, 2])
    x_flat_next = tf.reshape(tf.roll(x_all, shift=-1, axis=1), [K, 2])
    x_src  = tf.concat([ghost, x_flat],      axis=0)   # (K+1, 2)
    x_nxt  = tf.concat([ghost, x_flat_next], axis=0)   # (K+1, 2)

    # Source edge endpoints
    a0 = x_flat       # (K, 2) — start of source edge k
    a1 = x_flat_next  # (K, 2) — end   of source edge k

    # Global node index of a1 for each source edge k
    k_arr     = tf.range(K, dtype=tf.int32)
    e_arr     = k_arr % N
    p_arr     = k_arr // N
    a1_global = p_arr * N + (e_arr + 1) % N              # (K,)

    # Gather candidate edge endpoints from CapCandidates (1-indexed; 0=ghost)
    cand      = CapCandidates                             # (K, E) int32
    cand_flat = tf.reshape(cand, [K * E])
    b0 = tf.reshape(tf.gather(x_src, cand_flat), [K, E, 2])
    b1 = tf.reshape(tf.gather(x_nxt, cand_flat), [K, E, 2])

    # Periodic min-image: shift candidate edge endpoints (b0, b1) into the
    # nearest image of the source-edge anchor a0. b0 and b1 are on the same
    # particle (consecutive nodes), so they share the same shift — applying
    # it preserves the candidate edge vector ab = b1 - b0 unchanged.
    # box encoding: box[i] > 0 means periodic on axis i (length box[i]);
    # box[i] == 0 means non-periodic on that axis (no shift). For non-periodic
    # boundaries (open / hard walls / hopper) box=None and this block is skipped
    # entirely, so there is zero overhead in the non-periodic path.
    if box is not None:
        a0_e = tf.expand_dims(a0, axis=1)                    # (K, 1, 2)
        delta = b0 - a0_e                                    # (K, E, 2)
        box_b = tf.reshape(box, [1, 1, 2])
        # Per-axis: shift = -box * round(delta/box) when box>0, else 0
        periodic_mask = tf.cast(box_b > 0, dtype)            # (1,1,2) ∈ {0,1}
        box_safe = tf.where(box_b > 0, box_b, tf.ones_like(box_b))
        shift = -periodic_mask * box_safe * tf.round(delta / box_safe)
        b0 = b0 + shift
        b1 = b1 + shift

    # Safe 0-indexed candidate node addresses (ghost → 0, clamped)
    cand_0safe = tf.maximum(cand, 1) - 1                 # (K, E) 0-indexed, ghost→0
    e_b        = cand_0safe % N
    p_b        = cand_0safe // N
    b1_global  = p_b * N + (e_b + 1) % N                 # (K, E)

    # Active mask (non-ghost candidates)
    active = tf.cast(tf.not_equal(cand, 0), dtype)       # (K, E)

    # Candidate edge vectors
    ab  = b1 - b0                                         # (K, E, 2)
    ab2 = tf.reduce_sum(ab * ab, axis=2)                  # (K, E)

    # 2-pt Gauss quadrature params — tf.constant preserves float64 precision
    sg0 = tf.constant((1.0 - 1.0 / 3.0 ** 0.5) / 2.0, dtype=dtype)
    sg1 = tf.constant((1.0 + 1.0 / 3.0 ** 0.5) / 2.0, dtype=dtype)
    wg  = tf.constant(0.5, dtype=dtype)

    # Polydisperse contact radius per (source_edge, candidate_edge):
    # contact_r[k, e] = r_c_flat[k] + r_c_flat[cand_0safe[k, e]]
    r_c_cand = tf.gather(r_c_flat, tf.reshape(cand_0safe, [-1]))   # (K*E,)
    r_c_cand = tf.reshape(r_c_cand, [K, E])                         # (K, E)
    contact_r = r_c_flat[:, None] + r_c_cand                        # (K, E)

    _zero  = tf.constant(0.0,   dtype=dtype)
    _half  = tf.constant(0.5,   dtype=dtype)
    _one   = tf.constant(1.0,   dtype=dtype)
    _tiny  = tf.constant(1e-30, dtype=dtype)
    _eps   = tf.constant(1e-15, dtype=dtype)

    # Stage 1 scatter restructure: instead of 6 scatters (2× a1 source + 4×
    # Newton-3rd b0/b1 across both Gauss points), accumulate into:
    #   - f_a0     : source-side a0 contributions (direct add)
    #   - f_a1_src : source-side a1 contributions, indexed at the SOURCE edge,
    #                then rolled by +1 along the N axis to land at a1's node
    #                (no scatter needed — same as the a0 direct add pattern)
    #   - scatter_*_list : flat (index, contribution) pairs for the Newton-3rd
    #                deposits onto candidate edge endpoints b0 and b1, across
    #                BOTH Gauss points → ONE combined tensor_scatter_nd_add.
    # Net: 6 scatter_nd_add calls → 1.
    f_a0     = tf.zeros([K, 2], dtype=dtype)
    f_a1_src = tf.zeros([K, 2], dtype=dtype)
    scatter_indices_list = []
    scatter_contribs_list = []

    for sg, w0 in [(sg0, _one - sg0), (sg1, _one - sg1)]:
        xq     = w0 * a0 + sg * a1                             # (K, 2)
        xq_exp = xq[:, None, :]                                # (K, 1, 2)

        diff_b0 = xq_exp - b0                                  # (K, E, 2)
        # Periodic minimum-image: wrap diff_b0 so we measure through nearest image
        if box is not None:
            diff_b0_x = diff_b0[:, :, 0] - box[0] * tf.math.round(diff_b0[:, :, 0] / box[0])
            diff_b0_y = diff_b0[:, :, 1] - box[1] * tf.math.round(diff_b0[:, :, 1] / box[1])
            diff_b0 = tf.stack([diff_b0_x, diff_b0_y], axis=2)
        dot_num = tf.reduce_sum(diff_b0 * ab, axis=2)          # (K, E)
        t_B = tf.clip_by_value(
            tf.where(ab2 > _tiny,
                     dot_num / tf.maximum(ab2, _tiny),
                     tf.fill([K, E], _half)),
            _zero, _one)                                        # (K, E)

        # diff = vector from image-closest-point to xq (already in minimum-image frame)
        diff = diff_b0 - t_B[:, :, None] * ab                 # (K, E, 2)
        d    = tf.sqrt(tf.reduce_sum(diff * diff, axis=2))     # (K, E)
        gap  = d - contact_r                                    # (K, E)

        in_contact = tf.cast(gap < _zero, dtype) * active
        d_safe = tf.maximum(d, _eps)
        n_hat  = diff / d_safe[:, :, None]                     # (K, E, 2)

        # Polydisperse force magnitude: k_c comes from source edge's particle
        F_mag  = k_c_flat[:, None] * (-gap) * L0_flat[:, None] * wg   # (K, E)
        F_mag  = F_mag * in_contact                            # (K, E), ghost/far → 0
        F_vec  = F_mag[:, :, None] * n_hat                    # (K, E, 2)

        # Source (A) side: reduce_sum over E, weighted to edge endpoints
        F_src = tf.reduce_sum(F_vec, axis=1)                   # (K, 2)
        f_a0     = f_a0     + w0 * F_src                       # a0: direct add
        f_a1_src = f_a1_src + sg * F_src                       # a1: roll later

        # Newton 3rd — collect (index, contribution) for one combined scatter
        F_vec_flat = tf.reshape(F_vec, [K * E, 2])
        t_B_flat   = tf.reshape(t_B,   [K * E, 1])
        scatter_indices_list.append(tf.reshape(cand_0safe, [K * E, 1]))
        scatter_contribs_list.append(-(_one - t_B_flat) * F_vec_flat)
        scatter_indices_list.append(tf.reshape(b1_global,  [K * E, 1]))
        scatter_contribs_list.append(-t_B_flat * F_vec_flat)

    # Source a1 contribution: roll +1 along the N axis within each particle.
    # f_a1_src is indexed by source edge k = (p, e); a1's target node is
    # (p, (e+1) % N). After roll(+1, axis=N), entry [p, e] holds f_a1_src[p, e-1],
    # i.e. the contribution from edge (p, e-1) whose a1 IS node (p, e). Direct add.
    f_a1_src_3d = tf.reshape(f_a1_src, [P, N, 2])
    f_a1_3d     = tf.roll(f_a1_src_3d, shift=1, axis=1)
    f_a1        = tf.reshape(f_a1_3d, [K, 2])

    # Newton-3rd scatter REMOVED — PRCM lists each contact pair in both rows
    # (bidirectional), so the per-row reduce_sum on f_a0/f_a1 already captures
    # forces on both sides. The earlier scatter was double-counting forces.
    # NOTE: this assumes the candidacy provider lists pairs bidirectionally
    # (PRCM does). Production CandidacyManager uses pA<pB convention which
    # would now under-count — only PRCM is supported with this kernel.
    # all_indices  = tf.concat(scatter_indices_list,  axis=0)
    # all_contribs = tf.concat(scatter_contribs_list, axis=0)

    f_out = f_a0 + f_a1
    # f_out = tf.tensor_scatter_nd_add(f_out, all_indices, all_contribs)

    return tf.reshape(f_out, [P, N, 2])


# ── primitive boundary forces ─────────────────────────────────────────────────

def make_prim_data(prim_list, dtype=None):
    """
    Build prim_data dict from a list of (primitive_obj, k_pen_mult, vel) tuples.

    primitive_obj : LineSegment, Arc, Polygon, Box, or Hopper from contact_primitives
    k_pen_mult    : scalar — multiplier on particle's k_c (1.0 = use particle stiffness)
    vel           : (2,) array — constant translation velocity (or None for static)

    The force formula for wall contact is:
      F = k_c_particle * k_pen_mult * max(0, -(gap - r_c_particle - r_c_wall)) * L0 * wg
    matching the frozen capsule_shell.py::apply_primitive_forces exactly when k_pen_mult=1.0.

    Returns dict of float64 tensors for use with primitive_forces_tf().
    """
    from src.simulation.contact_primitives import LineSegment, Arc, Polygon, Hopper

    if dtype is None:
        dtype = NP_DTYPE

    seg_p0_list    = []
    seg_p1_list    = []
    seg_n_list     = []
    seg_k_pen_list = []  # k_pen multiplier (applied to k_c_flat)
    seg_r_c_list   = []  # wall capsule radius (usually 0)
    seg_vel_list   = []
    seg_omega_list = []  # angular velocity about r_ref (default 0)
    seg_r_ref_list = []  # reference point for rotation (default [0,0])
    seg_group_list = []  # -1 = standalone; >=0 = polygon group id

    arc_c_list     = []
    arc_R_list     = []
    arc_sgn_list   = []
    arc_k_pen_list = []
    arc_r_c_list   = []
    arc_vel_list   = []
    arc_omega_list = []  # angular velocity of arc centre about r_ref (default 0)
    arc_r_ref_list = []  # reference point for arc rotation (default [0,0])

    group_id = 0

    seg_osc_A_list     = []
    seg_osc_omega_list = []
    seg_osc_sign_list  = []

    def _add_seg(p0, p1, normal, k_pen, r_c, vel, omega, r_ref, group=-1,
                 osc_A=0.0, osc_omega=0.0, osc_sign=0.0):
        seg_p0_list.append(np.asarray(p0, dtype=float))
        seg_p1_list.append(np.asarray(p1, dtype=float))
        seg_n_list.append(np.asarray(normal, dtype=float))
        seg_k_pen_list.append(float(k_pen))
        seg_r_c_list.append(float(r_c))
        seg_vel_list.append(np.asarray(vel, dtype=float))
        seg_omega_list.append(float(omega))
        seg_r_ref_list.append(np.asarray(r_ref, dtype=float))
        seg_group_list.append(int(group))
        seg_osc_A_list.append(float(osc_A))
        seg_osc_omega_list.append(float(osc_omega))
        seg_osc_sign_list.append(float(osc_sign))

    for entry in prim_list:
        # Tuple formats (all backward-compatible):
        #   3: (obj, k_pen, vel)
        #   4: (obj, k_pen, vel, r_c_wall)
        #   5: (obj, k_pen, vel, omega, r_ref)
        #   6: (obj, k_pen, vel, r_c_wall, omega, r_ref)
        #   7: (obj, k_pen, vel, r_c_wall, omega, r_ref, osc)
        #      osc = (A, omega_osc, sign) for sinusoidal oscillation, or None
        osc_A_v = osc_omega_v = osc_sign_v = 0.0
        if len(entry) == 3:
            prim_obj, k_pen, vel = entry
            r_c_wall, omega, r_ref = 0.0, 0.0, np.zeros(2)
        elif len(entry) == 4:
            prim_obj, k_pen, vel, r_c_wall = entry
            omega, r_ref = 0.0, np.zeros(2)
        elif len(entry) == 5:
            prim_obj, k_pen, vel, omega, r_ref = entry
            r_c_wall = 0.0
        elif len(entry) == 6:
            prim_obj, k_pen, vel, r_c_wall, omega, r_ref = entry
        else:
            prim_obj, k_pen, vel, r_c_wall, omega, r_ref, osc = entry
            if osc is not None:
                osc_A_v, osc_omega_v, osc_sign_v = float(osc[0]), float(osc[1]), float(osc[2])

        vel   = np.zeros(2)  if vel   is None else np.asarray(vel,   dtype=float)
        r_ref = np.zeros(2)  if r_ref is None else np.asarray(r_ref, dtype=float)
        omega = float(omega) if omega is not None else 0.0

        if isinstance(prim_obj, Arc):
            sgn = 1.0 if prim_obj.convex else -1.0
            arc_c_list.append(prim_obj.center.copy())
            arc_R_list.append(prim_obj.radius)
            arc_sgn_list.append(sgn)
            arc_k_pen_list.append(float(k_pen))
            arc_r_c_list.append(float(r_c_wall))
            arc_vel_list.append(vel.copy())
            arc_omega_list.append(omega)
            arc_r_ref_list.append(r_ref.copy())

        elif isinstance(prim_obj, LineSegment):
            _add_seg(prim_obj.p0, prim_obj.p1, prim_obj.normal,
                     k_pen, r_c_wall, vel, omega, r_ref, group=-1,
                     osc_A=osc_A_v, osc_omega=osc_omega_v, osc_sign=osc_sign_v)

        elif isinstance(prim_obj, Hopper):
            _add_seg(prim_obj._left.p0,  prim_obj._left.p1,  prim_obj._left.normal,
                     k_pen, r_c_wall, vel, omega, r_ref, group=group_id,
                     osc_A=osc_A_v, osc_omega=osc_omega_v, osc_sign=osc_sign_v)
            _add_seg(prim_obj._right.p0, prim_obj._right.p1, prim_obj._right.normal,
                     k_pen, r_c_wall, vel, omega, r_ref, group=group_id,
                     osc_A=osc_A_v, osc_omega=osc_omega_v, osc_sign=osc_sign_v)
            group_id += 1

        elif isinstance(prim_obj, Polygon):
            for (p0, p1, seg_normal) in prim_obj._segments:
                _add_seg(p0, p1, seg_normal, k_pen, r_c_wall, vel, omega, r_ref,
                         group=group_id)
            group_id += 1

        else:
            raise ValueError(f"Unknown primitive type: {type(prim_obj)}")

    Ng = group_id   # number of polygon groups

    # Build polygon grouping tables
    if Ng > 0:
        group_segs = [[] for _ in range(Ng)]
        for s_idx, gid in enumerate(seg_group_list):
            if gid >= 0:
                group_segs[gid].append(s_idx)
        Ms = max(len(g) for g in group_segs)

        grp_seg  = np.full((Ng, Ms), -1, dtype=np.int32)
        grp_mask = np.zeros((Ng, Ms), dtype=bool)
        for g, segs in enumerate(group_segs):
            for slot, s_idx in enumerate(segs):
                grp_seg[g, slot] = s_idx
                grp_mask[g, slot] = True
    else:
        grp_seg  = np.zeros((0, 1), dtype=np.int32)
        grp_mask = np.zeros((0, 1), dtype=bool)
        Ms = 1

    # Precompute inverse mapping: seg_grp_flat_idx[s] = g*Ms + k  (or -1 for standalone)
    Ls = len(seg_p0_list)
    seg_grp_flat_idx = np.full(Ls, -1, dtype=np.int32)
    if Ng > 0:
        for g in range(Ng):
            for k in range(Ms):
                s = grp_seg[g, k]
                if s >= 0:
                    seg_grp_flat_idx[s] = g * Ms + k

    # Use dummy entries so shapes are never empty (avoids TF trace issues)
    if Ls == 0:
        seg_p0_list    = [np.zeros(2)];       seg_p1_list    = [np.zeros(2)]
        seg_n_list     = [np.array([0.,1.])]; seg_k_pen_list = [0.0]
        seg_r_c_list   = [0.0];               seg_vel_list   = [np.zeros(2)]
        seg_omega_list = [0.0];               seg_r_ref_list = [np.zeros(2)]
        seg_group_list = [-1];                seg_grp_flat_idx = np.array([-1], dtype=np.int32)
        seg_osc_A_list = [0.0]; seg_osc_omega_list = [0.0]; seg_osc_sign_list = [0.0]
        Ls = 1

    La = len(arc_c_list)
    if La == 0:
        arc_c_list     = [np.zeros(2)];  arc_R_list     = [1.0]
        arc_sgn_list   = [1.0];          arc_k_pen_list = [0.0]
        arc_r_c_list   = [0.0];          arc_vel_list   = [np.zeros(2)]
        arc_omega_list = [0.0];          arc_r_ref_list = [np.zeros(2)]
        La = 1

    def _f(lst): return np.array(lst, dtype=dtype)

    # Watchdog support: precompute max wall velocity per primitive (over a chunk).
    # Max |v| for a segment: translation + osc-peak + |omega|·max_arm_length.
    seg_p0_arr = np.array(seg_p0_list, dtype=np.float64)
    seg_p1_arr = np.array(seg_p1_list, dtype=np.float64)
    seg_r_ref_arr = np.array(seg_r_ref_list, dtype=np.float64)
    seg_vel_arr = np.array(seg_vel_list, dtype=np.float64)
    seg_omega_arr = np.array(seg_omega_list, dtype=np.float64)
    seg_osc_A_arr = np.array(seg_osc_A_list, dtype=np.float64)
    seg_osc_omega_arr = np.array(seg_osc_omega_list, dtype=np.float64)
    seg_arm_max = np.maximum(
        np.linalg.norm(seg_p0_arr - seg_r_ref_arr, axis=1),
        np.linalg.norm(seg_p1_arr - seg_r_ref_arr, axis=1))                  # (Ls,)
    seg_v_max = (np.linalg.norm(seg_vel_arr, axis=1)
                 + np.abs(seg_omega_arr) * seg_arm_max
                 + np.abs(seg_osc_A_arr * seg_osc_omega_arr))                # (Ls,)

    arc_R_arr     = np.array(arc_R_list, dtype=np.float64)
    arc_vel_arr   = np.array(arc_vel_list, dtype=np.float64)
    arc_omega_arr = np.array(arc_omega_list, dtype=np.float64)
    arc_v_max = (np.linalg.norm(arc_vel_arr, axis=1)
                 + np.abs(arc_omega_arr) * arc_R_arr)                        # (La,)

    prim_data = dict(
        # Line segments
        seg_p0          = tf.constant(_f(seg_p0_list)),           # (Ls, 2)
        seg_p1          = tf.constant(_f(seg_p1_list)),           # (Ls, 2)
        seg_n           = tf.constant(_f(seg_n_list)),            # (Ls, 2) reference normals
        seg_k_pen       = tf.constant(_f(seg_k_pen_list)),        # (Ls,)
        seg_r_c         = tf.constant(_f(seg_r_c_list)),          # (Ls,)
        seg_vel         = tf.constant(_f(seg_vel_list)),          # (Ls, 2) translation vel
        seg_omega       = tf.constant(_f(seg_omega_list)),        # (Ls,)  angular vel
        seg_r_ref       = tf.constant(_f(seg_r_ref_list)),        # (Ls, 2) rotation pivot
        seg_group       = tf.constant(np.array(seg_group_list, dtype=np.int32)),
        seg_grp_flat_idx= tf.constant(seg_grp_flat_idx),
        # Sinusoidal oscillation: y(t) = y0 + osc_A * osc_sign * sin(osc_omega * t)
        seg_osc_A       = tf.constant(_f(seg_osc_A_list)),        # (Ls,) amplitude
        seg_osc_omega   = tf.constant(_f(seg_osc_omega_list)),    # (Ls,) angular freq
        seg_osc_sign    = tf.constant(_f(seg_osc_sign_list)),     # (Ls,) ±1
        # Polygon groups
        grp_seg         = tf.constant(grp_seg),
        grp_mask        = tf.constant(grp_mask),
        # Arcs
        arc_c           = tf.constant(_f(arc_c_list)),            # (La, 2) reference centres
        arc_R           = tf.constant(_f(arc_R_list)),            # (La,)
        arc_sgn         = tf.constant(_f(arc_sgn_list)),          # (La,)
        arc_k_pen       = tf.constant(_f(arc_k_pen_list)),        # (La,)
        arc_r_c         = tf.constant(_f(arc_r_c_list)),          # (La,)
        arc_vel         = tf.constant(_f(arc_vel_list)),          # (La, 2) translation vel
        arc_omega       = tf.constant(_f(arc_omega_list)),        # (La,)  angular vel
        arc_r_ref       = tf.constant(_f(arc_r_ref_list)),        # (La, 2) rotation pivot
        # Watchdog: peak wall velocity per primitive (precomputed; static)
        seg_v_max       = tf.constant(_f(seg_v_max.tolist())),    # (Ls,)
        arc_v_max       = tf.constant(_f(arc_v_max.tolist())),    # (La,)
    )
    return prim_data


@tf.function(jit_compile=True)
def primitive_forces_tf(x_all, t, prim_data, r_c_flat, k_c_flat, L0_flat,
                          box=None):
    """
    Nodal forces from rigid boundary primitives (line segments, polygons, arcs).

    Matches frozen capsule_shell.py::apply_primitive_forces exactly when
    prim_data seg_k_pen = 1.0 and seg_r_c = 0.0.

    Uses 2-point Gauss quadrature over source edges (same as inter-capsule forces).

    x_all    : (P, N, 2) DTYPE
    t        : scalar DTYPE — current simulation time
    prim_data: dict from make_prim_data()
    r_c_flat : (K,) DTYPE — per-edge capsule radius
    k_c_flat : (K,) DTYPE — per-edge contact stiffness
    L0_flat  : (K,) DTYPE — per-edge reference length
    box      : (2,) DTYPE [Lx, Ly] or None — periodic box; per-axis encoding
               (box[i] > 0 means periodic on axis i; 0 means non-periodic).
               When box is None, this function behaves identically to the
               non-periodic version (zero overhead — Python guard skips the
               min-image block).

    Returns f : (P, N, 2) DTYPE
    """
    dtype = x_all.dtype
    P = x_all.shape[0]
    N = x_all.shape[1]
    K = P * N

    # Source edge endpoints
    a0_flat = tf.reshape(x_all, [K, 2])                              # (K, 2)
    a1_flat = tf.reshape(tf.roll(x_all, shift=-1, axis=1), [K, 2])  # (K, 2)

    # a1 global node index for scatter
    k_arr = tf.range(K, dtype=tf.int32)
    e_arr = k_arr % N
    p_arr = k_arr // N
    a1_global = p_arr * N + (e_arr + 1) % N   # (K,)

    f_out = tf.zeros([K, 2], dtype=dtype)

    _zero = tf.constant(0.0,   dtype=dtype)
    _one  = tf.constant(1.0,   dtype=dtype)
    _tiny = tf.constant(1e-30, dtype=dtype)
    _eps  = tf.constant(1e-15, dtype=dtype)
    _INF  = tf.constant(1e30,  dtype=dtype)
    _wg   = tf.constant(0.5,   dtype=dtype)
    sg0   = tf.constant((1.0 - 1.0 / 3.0**0.5) / 2.0, dtype=dtype)
    sg1   = tf.constant((1.0 + 1.0 / 3.0**0.5) / 2.0, dtype=dtype)

    # Current segment endpoint positions (translation + rotation about r_ref)
    seg_omega = prim_data['seg_omega']   # (Ls,)
    seg_r_ref = prim_data['seg_r_ref']   # (Ls, 2)
    cos_s = tf.cos(seg_omega * t)        # (Ls,)
    sin_s = tf.sin(seg_omega * t)        # (Ls,)
    # Rotate p0 and p1 about r_ref, then translate
    def _rotate_pts(pts_ref):            # pts_ref: (Ls, 2)
        dp = pts_ref - seg_r_ref         # (Ls, 2)
        rot_x = cos_s * dp[:,0] - sin_s * dp[:,1]
        rot_y = sin_s * dp[:,0] + cos_s * dp[:,1]
        return seg_r_ref + tf.stack([rot_x, rot_y], axis=1) + prim_data['seg_vel'] * t
    seg_p0 = _rotate_pts(prim_data['seg_p0'])   # (Ls, 2)
    seg_p1 = _rotate_pts(prim_data['seg_p1'])   # (Ls, 2)
    # Normal also rotates
    dn = prim_data['seg_n']                       # (Ls, 2)
    seg_n = tf.stack([cos_s * dn[:,0] - sin_s * dn[:,1],
                      sin_s * dn[:,0] + cos_s * dn[:,1]], axis=1)  # (Ls, 2)

    # Sinusoidal oscillation: shift y-coordinate by A*sign*sin(omega*t)
    # Non-oscillating walls have osc_A=0, so this is a no-op for them.
    osc_dy = (prim_data['seg_osc_A'] * prim_data['seg_osc_sign']
              * tf.sin(prim_data['seg_osc_omega'] * t))             # (Ls,)
    zeros  = tf.zeros_like(osc_dy)
    seg_p0 = seg_p0 + tf.stack([zeros, osc_dy], axis=1)
    seg_p1 = seg_p1 + tf.stack([zeros, osc_dy], axis=1)

    Ls     = seg_p0.shape[0]

    # Segment vectors and squared lengths
    ab_seg  = seg_p1 - seg_p0           # (Ls, 2)
    ab2_seg = tf.reduce_sum(ab_seg * ab_seg, axis=1)  # (Ls,)

    # Current arc centres (translation + rotation of centre about r_ref)
    arc_omega = prim_data['arc_omega']   # (La,)
    arc_r_ref = prim_data['arc_r_ref']   # (La, 2)
    cos_a = tf.cos(arc_omega * t)
    sin_a = tf.sin(arc_omega * t)
    dc = prim_data['arc_c'] - arc_r_ref  # (La, 2)
    arc_c_rot = arc_r_ref + tf.stack([cos_a * dc[:,0] - sin_a * dc[:,1],
                                       sin_a * dc[:,0] + cos_a * dc[:,1]], axis=1)
    arc_c = arc_c_rot + prim_data['arc_vel'] * t   # (La, 2)
    La    = arc_c.shape[0]

    # Polygon group tables
    Ng = prim_data['grp_seg'].shape[0]
    Ms = prim_data['grp_seg'].shape[1]

    for sg in [sg0, sg1]:
        w0 = _one - sg   # Gauss weight for a0 node

        # Quadrature point for each source edge
        xq = w0 * a0_flat + sg * a1_flat   # (K, 2)

        # ── Line segments ──────────────────────────────────────────────────────
        xq_exp = xq[:, None, :]            # (K, 1, 2)
        p0_exp = seg_p0[None, :, :]        # (1, Ls, 2)
        ab_exp = ab_seg[None, :, :]        # (1, Ls, 2)

        # Periodic min-image for source-particle Gauss point relative to each
        # segment endpoint p0. Shifts xq into the nearest image of p0 per
        # (k, s) pair. p0, p1 stay raw — particle's image moves to meet the wall.
        if box is not None:
            box_b   = tf.reshape(box, [1, 1, 2])
            pmask   = tf.cast(box_b > 0, dtype)
            box_safe = tf.where(box_b > 0, box_b, tf.ones_like(box_b))
            delta_pre = xq_exp - p0_exp                                      # (K, Ls, 2)
            shift_xq  = -pmask * box_safe * tf.round(delta_pre / box_safe)   # (K, Ls, 2)
            xq_exp = xq_exp + shift_xq                                       # (K, Ls, 2)

        diff0  = xq_exp - p0_exp           # (K, Ls, 2)
        dot_n  = tf.reduce_sum(diff0 * ab_exp, axis=2)           # (K, Ls)
        t_raw  = dot_n / tf.maximum(ab2_seg[None, :], _tiny)     # (K, Ls) unclamped
        t_seg  = tf.clip_by_value(t_raw, _zero, _one)            # (K, Ls)

        # Interior-projection mask: only fire contact when the perpendicular foot
        # lands on the segment (t_raw ∈ [0,1]).  Prevents endpoint clamping from
        # generating spurious long-range forces on particles that are far off the
        # end of a finite wall.  Applied to standalone segments only; polygon-group
        # segments use their own argmax-based active-face selection instead.
        t_in_seg   = tf.cast((t_raw >= _zero) & (t_raw <= _one), dtype)  # (K, Ls)
        is_standalone = tf.cast(
            prim_data['seg_grp_flat_idx'] < 0, dtype)[None, :]           # (1, Ls)
        seg_in_seg = _one - is_standalone + is_standalone * t_in_seg     # (K, Ls)

        cp_seg = p0_exp + t_seg[:, :, None] * ab_exp             # (K, Ls, 2)
        diff_s = xq_exp - cp_seg                                  # (K, Ls, 2)
        dist_s = tf.sqrt(tf.reduce_sum(diff_s * diff_s, axis=2)) # (K, Ls)
        n_exp  = seg_n[None, :, :]                               # (1, Ls, 2)
        gap_s  = tf.reduce_sum(diff_s * n_exp, axis=2)           # (K, Ls)

        # ── Polygon group masking ──────────────────────────────────────────────
        # For a convex polygon obstacle, the correct active segment is the one
        # with the MAXIMUM gap (the polygon SDF = max gap over all half-planes).
        # Using argmax(gap) instead of argmin(dist) prevents spurious contact
        # forces on outside nodes near corners: at a corner, the nearest-by-
        # distance wall may have a negative gap even for geometrically-outside
        # nodes, while the correct wall (whose half-plane the node is truly
        # outside of) has a positive gap.
        #   gap > 0  →  node outside that half-plane (SDF > 0 = outside polygon)
        #   gap ≤ 0  →  node inside all half-planes (SDF ≤ 0 = inside polygon)
        seg_active = tf.ones([K, Ls], dtype=dtype)

        if Ng > 0:
            # Gather gap for each (source_edge, polygon_group, slot)
            grp_seg_flat = tf.reshape(prim_data['grp_seg'], [-1])       # (Ng*Ms,)
            grp_seg_safe = tf.maximum(grp_seg_flat, 0)                  # clamp -1→0
            gap_grp = tf.gather(gap_s, grp_seg_safe, axis=1)           # (K, Ng*Ms)
            gap_grp = tf.reshape(gap_grp, [K, Ng, Ms])                 # (K, Ng, Ms)

            # Mask inactive slots (padded -1 entries) with -INF
            mask_exp = prim_data['grp_mask'][None, :, :]               # (1, Ng, Ms)
            gap_grp_masked = tf.where(mask_exp,
                                      gap_grp,
                                      tf.fill(tf.shape(gap_grp), -_INF))

            # Argmax(gap) within each polygon group → one-hot active slot
            argmax_poly  = tf.argmax(gap_grp_masked, axis=2, output_type=tf.int32)  # (K, Ng)
            active_slots = (tf.one_hot(argmax_poly, Ms, dtype=dtype)
                            * tf.cast(mask_exp, dtype))                # (K, Ng, Ms)
            active_flat  = tf.reshape(active_slots, [K, Ng * Ms])     # (K, Ng*Ms)

            # Map each segment to its (group, slot) flat index
            seg_grp_idx  = prim_data['seg_grp_flat_idx']              # (Ls,) int32
            is_poly      = tf.cast(seg_grp_idx >= 0, dtype)           # (Ls,)
            grp_safe     = tf.maximum(seg_grp_idx, 0)                 # clamp -1→0
            poly_act     = tf.gather(active_flat, grp_safe, axis=1)   # (K, Ls)

            # Standalone segments always active; polygon segments use poly_act
            seg_active = (_one - is_poly[None, :]) + is_poly[None, :] * poly_act

        # ── Apply segment forces ───────────────────────────────────────────────
        # gap_capsule = gap_s - r_c_flat[k] - seg_r_c[s]
        r_c_eff  = r_c_flat[:, None] + prim_data['seg_r_c'][None, :]   # (K, Ls)
        gap_eff  = gap_s - r_c_eff                                      # (K, Ls)
        in_c     = tf.cast(gap_eff < _zero, dtype) * seg_active * seg_in_seg  # (K, Ls)

        # F per Gauss pt: k_c_particle * k_pen_mult * max(0,-gap_eff) * L0 * wg
        k_eff    = k_c_flat[:, None] * prim_data['seg_k_pen'][None, :] # (K, Ls)
        F_mag    = k_eff * tf.maximum(-gap_eff, _zero) * L0_flat[:, None] * _wg  # (K, Ls)
        F_mag    = F_mag * in_c                                         # (K, Ls)

        # Normal pointing from primitive surface toward particle
        F_vec    = F_mag[:, :, None] * n_exp                            # (K, Ls, 2)
        F_seg    = tf.reduce_sum(F_vec, axis=1)                         # (K, 2)

        # ── Arcs ──────────────────────────────────────────────────────────────
        c_exp   = arc_c[None, :, :]       # (1, La, 2)
        xq_arc  = xq[:, None, :]          # (K, 1, 2)
        # Periodic min-image: shift xq into nearest image of each arc center.
        if box is not None:
            box_b   = tf.reshape(box, [1, 1, 2])
            pmask   = tf.cast(box_b > 0, dtype)
            box_safe = tf.where(box_b > 0, box_b, tf.ones_like(box_b))
            delta_pre = xq_arc - c_exp                                       # (K, La, 2)
            shift_xq  = -pmask * box_safe * tf.round(delta_pre / box_safe)   # (K, La, 2)
            xq_arc = xq_arc + shift_xq
        diff_a  = xq_arc - c_exp          # (K, La, 2)
        dist_a  = tf.sqrt(tf.reduce_sum(diff_a * diff_a, axis=2))       # (K, La)

        # gap: convex (+1): dist - R; concave (-1): R - dist
        # Both: gap = arc_sgn * (dist - R) ... verified:
        #   convex sgn=+1: gap = dist-R (positive outside = safe) ✓
        #   concave sgn=-1: gap = -(dist-R) = R-dist (positive inside = safe) ✓
        arc_sgn  = prim_data['arc_sgn']
        gap_a    = arc_sgn[None, :] * (dist_a - prim_data['arc_R'][None, :])  # (K, La)

        # Normal: convex → radial outward; concave → radial inward
        d_safe  = tf.maximum(dist_a, _eps)
        radial  = diff_a / d_safe[:, :, None]             # (K, La, 2)
        n_arc   = arc_sgn[None, :, None] * radial          # (K, La, 2)

        r_c_a   = r_c_flat[:, None] + prim_data['arc_r_c'][None, :]    # (K, La)
        gap_a_eff = gap_a - r_c_a                                       # (K, La)
        in_ca   = tf.cast(gap_a_eff < _zero, dtype)                    # (K, La)

        k_eff_a = k_c_flat[:, None] * prim_data['arc_k_pen'][None, :]  # (K, La)
        F_mag_a = k_eff_a * tf.maximum(-gap_a_eff, _zero) * L0_flat[:, None] * _wg  # (K, La)
        F_mag_a = F_mag_a * in_ca

        F_arc   = tf.reduce_sum(F_mag_a[:, :, None] * n_arc, axis=1)  # (K, 2)

        # ── Accumulate (Gauss point contribution) ────────────────────────────
        F_total_gp = F_seg + F_arc   # (K, 2)

        # Distribute to edge nodes: a0 gets w0 fraction, a1 gets sg fraction
        f_out = f_out + w0 * F_total_gp
        f_out = tf.tensor_scatter_nd_add(
            f_out, tf.expand_dims(a1_global, 1), sg * F_total_gp)

    return tf.reshape(f_out, [P, N, 2])


# ── force evaluation (no integration step) ───────────────────────────────────

def eval_forces_tf(state, CapCandidates, alpha_damp, g, params,
                   t=None, prim_data=None):
    """
    Compute per-node force breakdown without advancing the simulation.

    Returns a dict of (P, N, 2) float64 arrays:
      f_elastic  — internal elastic + area + line-tension forces
      f_contact  — inter-capsule and wall contact forces
      f_drag     — Stokes drag forces (zero if xi_drag not set)
      f_reg      — tangential regularisation pseudo-force (emulsion only)
      f_total    — sum of all four
    """
    x        = state['x_all']
    x_cm     = state['x_cm']
    v_cm     = state['v_cm']
    theta    = state['theta']
    omega    = state['omega']
    u_dot    = state['u_dot']
    M_disk   = params['M_disk']
    N_nodes  = x.shape[1]

    r_c_flat = tf.repeat(params['r_c_per_p'], N_nodes)
    k_c_flat = tf.repeat(params['k_c_per_p'], N_nodes)
    L0_flat  = tf.repeat(params['L0'],        N_nodes)

    if t is None:
        t = tf.constant(0.0, dtype=DTYPE)

    # contact
    f_contact = inter_capsule_forces_tf(
        x, CapCandidates, r_c_flat, k_c_flat, L0_flat,
        box=params.get('box', None))
    if prim_data is not None:
        f_contact = f_contact + primitive_forces_tf(
            x, t, prim_data, r_c_flat, k_c_flat, L0_flat,
            box=params.get('box', None))

    # elastic + regularisation
    f_elastic = internal_forces_tf(x, params, g)
    f_reg     = k_reg_forces_tf(x, params)

    # drag
    xi_drag = params.get('xi_drag_per_p', tf.zeros_like(M_disk))
    xi_any  = tf.reduce_any(xi_drag > 0.0)
    _alpha_pp  = params.get('alpha_damp_per_p', tf.zeros_like(M_disk))
    _alpha_any = tf.reduce_any(_alpha_pp > 0.0)
    alpha_p    = tf.where(_alpha_any, _alpha_pp,
                          tf.fill(tf.shape(M_disk), tf.cast(alpha_damp, DTYPE)))

    def _drag():
        r      = x - x_cm[:, None, :]
        v_rot  = omega[:, None, None] * tf.stack([-r[:,:,1], r[:,:,0]], axis=2)
        c_th   = tf.cos(theta)[:, None]
        s_th   = tf.sin(theta)[:, None]
        v_el_x = c_th * u_dot[:,:,0] - s_th * u_dot[:,:,1]
        v_el_y = s_th * u_dot[:,:,0] + c_th * u_dot[:,:,1]
        v_node = v_cm[:, None, :] + v_rot + tf.stack([v_el_x, v_el_y], axis=2)
        U_bg   = _eval_U_bg(x, t, params)
        x_next = tf.roll(x, shift=-1, axis=1)
        x_prev = tf.roll(x,  shift=1, axis=1)
        dL     = 0.5 * tf.norm(x_next - x_prev, axis=2)
        return -xi_drag[:, None, None] * (v_node - U_bg) * dL[:,:,None]

    f_drag  = tf.cond(xi_any, _drag, lambda: tf.zeros_like(x))
    f_total = f_elastic + f_reg + f_contact + f_drag

    return {
        'f_elastic': f_elastic.numpy(),
        'f_contact': f_contact.numpy(),
        'f_drag':    f_drag.numpy(),
        'f_reg':     f_reg.numpy(),
        'f_total':   f_total.numpy(),
    }


# ── fused full step (contact + primitives + internal + RB integration) ────────

@tf.function(jit_compile=True)
def step_full_tf(state, CapCandidates, dt, alpha_damp, g, params,
                 t=None, prim_data=None):
    """
    Full step: inter-capsule + primitive forces, then step_rb_tf.

    state         : dict with x_all, x_cm, v_cm, theta, omega, u, u_dot, X_ref
    CapCandidates : (K, E) int32
    dt, alpha_damp, g : scalar DTYPE
    params        : dict from make_state() — includes r_c_per_p, k_c_per_p
    t             : scalar DTYPE — current time (for moving primitives); None → 0
    prim_data     : dict from make_prim_data() — None → no wall forces

    Returns new_state, metrics.
    """
    N_nodes  = state['x_all'].shape[1]

    # Per-edge flat arrays (polydisperse support)
    r_c_flat = tf.repeat(params['r_c_per_p'], N_nodes)   # (K,)
    k_c_flat = tf.repeat(params['k_c_per_p'], N_nodes)   # (K,)
    L0_flat  = tf.repeat(params['L0'],        N_nodes)   # (K,)

    if t is None:
        t = tf.constant(0.0, dtype=state['x_all'].dtype)

    box = params.get('box', None)
    f_contact = inter_capsule_forces_tf(
        state['x_all'], CapCandidates, r_c_flat, k_c_flat, L0_flat, box=box)

    if prim_data is not None:
        f_prim    = primitive_forces_tf(
            state['x_all'], t, prim_data, r_c_flat, k_c_flat, L0_flat,
            box=box)
        f_contact = f_contact + f_prim

    new_state, metrics = step_rb_tf(state, f_contact, dt, alpha_damp, g, params, t=t)

    # ──────────────────────────────────────────────────────────────────────
    # Closing-rate watchdog (replaces absolute-displacement criterion).
    # Symmetric definition: per contact pair, the change in separation per
    # step, normalised by the pair's combined contact radius.  Captures both
    # particle-particle and particle-wall contacts; insensitive to bulk
    # translation (e.g. all particles falling together).
    # ──────────────────────────────────────────────────────────────────────
    dtype = state['x_all'].dtype
    P_  = state['x_all'].shape[0]
    N_  = N_nodes
    K_  = P_ * N_
    E_  = CapCandidates.shape[1]
    disp_step = new_state['x_all'] - state['x_all']               # (P, N, 2)
    disp_flat = tf.reshape(disp_step, [K_, 2])                     # (K, 2)

    # Particle-particle: per (k, e), |Δdisp_a - Δdisp_b| / (r_c_a + r_c_b)
    cand        = CapCandidates                                    # (K, E)
    cand_0safe  = tf.maximum(cand, 1) - 1                          # (K, E)
    disp_b_flat = tf.gather(disp_flat, tf.reshape(cand_0safe, [-1]))   # (K*E, 2)
    disp_b      = tf.reshape(disp_b_flat, [K_, E_, 2])
    rel_disp    = disp_b - disp_flat[:, None, :]                   # (K, E, 2)
    rel_mag     = tf.sqrt(tf.reduce_sum(rel_disp * rel_disp, axis=2))   # (K, E)
    r_c_cand    = tf.gather(r_c_flat, tf.reshape(cand_0safe, [-1]))
    r_c_cand    = tf.reshape(r_c_cand, [K_, E_])
    r_c_pair    = r_c_flat[:, None] + r_c_cand                     # (K, E)
    active      = tf.cast(tf.not_equal(cand, 0), dtype)            # (K, E)
    ratio_pp    = tf.reduce_max(rel_mag / tf.maximum(r_c_pair, tf.constant(1e-30, dtype=dtype))
                                * active)                          # scalar

    # Particle-wall (segments + arcs): conservative upper bound on closing rate
    # per (node k, primitive s):
    #   numerator   = |disp_node[k]| + |v_wall[s]|·dt
    #   denominator = r_c_particle[k] + r_c_wall[s]
    # (wall r_c is usually 0, so denominator is dominated by particle's r_c).
    if prim_data is not None:
        disp_mag_flat = tf.sqrt(tf.reduce_sum(disp_flat * disp_flat, axis=1))  # (K,)
        tiny = tf.constant(1e-30, dtype=dtype)

        seg_v_max = prim_data['seg_v_max']           # (Ls,)
        seg_r_c   = prim_data['seg_r_c']             # (Ls,)
        seg_num   = disp_mag_flat[:, None] + seg_v_max[None, :] * dt   # (K, Ls)
        seg_den   = r_c_flat[:, None] + seg_r_c[None, :]               # (K, Ls)
        seg_active = tf.cast(prim_data['seg_k_pen'] > 0, dtype)[None, :]
        seg_ratio_max = tf.reduce_max(seg_num / tf.maximum(seg_den, tiny) * seg_active)

        arc_v_max = prim_data['arc_v_max']           # (La,)
        arc_r_c   = prim_data['arc_r_c']             # (La,)
        arc_num   = disp_mag_flat[:, None] + arc_v_max[None, :] * dt   # (K, La)
        arc_den   = r_c_flat[:, None] + arc_r_c[None, :]               # (K, La)
        arc_active = tf.cast(prim_data['arc_k_pen'] > 0, dtype)[None, :]
        arc_ratio_max = tf.reduce_max(arc_num / tf.maximum(arc_den, tiny) * arc_active)

        ratio_wall = tf.maximum(seg_ratio_max, arc_ratio_max)
    else:
        ratio_wall = tf.constant(0.0, dtype=dtype)

    metrics['max_closing_ratio'] = tf.maximum(ratio_pp, ratio_wall)
    return new_state, metrics


# ── state construction from CapsuleParticle objects ───────────────────────────

def make_state(particles):
    """
    Build TF state + params from a list of CapsuleParticle objects.
    Supports polydisperse particles (different R0 per particle).
    All tensors cast to DTYPE.
    """
    def _t1(vals):   # scalar per particle → (P,)
        return tf.constant(np.array(vals, dtype=NP_DTYPE))
    def _t2(arrs):   # (N,2) per particle → (P,N,2)
        return tf.constant(np.stack(arrs).astype(NP_DTYPE))
    def _t1n(arrs):  # (N,) per particle → (P,N)
        return tf.constant(np.stack(arrs).astype(NP_DTYPE))

    state = dict(
        x_all = _t2([p.x     for p in particles]),
        x_cm  = tf.constant(np.stack([p.x_cm  for p in particles]).astype(NP_DTYPE)),
        v_cm  = tf.constant(np.stack([p.v_cm  for p in particles]).astype(NP_DTYPE)),
        theta = _t1([p.theta for p in particles]),
        omega = _t1([p.omega for p in particles]),
        u     = _t2([p.u     for p in particles]),
        u_dot = _t2([p.u_dot for p in particles]),
        X_ref = _t2([p.X_ref for p in particles]),
    )
    P = len(particles)
    params = dict(
        L0       = _t1([p.L0     for p in particles]),
        A0       = _t1([p.A0     for p in particles]),
        K_area   = _t1([p.K_area for p in particles]),
        El_t     = _t1([p.El_t   for p in particles]),
        EI       = _t1([p.EI     for p in particles]),
        m_node   = _t1([p.m_node for p in particles]),
        M_disk   = _t1([p.M_disk for p in particles]),
        I_disk   = _t1([p.I_disk for p in particles]),
        rho_d    = _t1([p.rho_d  for p in particles]),
        theta0   = _t1n([p.theta0 for p in particles]),
        # Polydisperse per-particle arrays
        r_c_per_p = _t1([p.r_c for p in particles]),   # (P,)
        k_c_per_p = _t1([p.k_c for p in particles]),   # (P,)
        # Line tension per particle (emulsion; 0 for elastic capsules)
        gamma_lt  = _t1([float(getattr(p, 'gamma', 0.0)) for p in particles]),
        k_reg_per_p = _t1([float(getattr(p, 'k_reg', 0.0)) for p in particles]),
        # Driven particle / periodic BC — all off by default (backward compatible)
        driven_mask  = tf.zeros([P], dtype=DTYPE),
        shape_frozen = tf.zeros([P], dtype=DTYPE),
        traj         = tf.zeros([P, 18], dtype=DTYPE),
        box          = tf.constant([1e9, 1e9], dtype=DTYPE),  # non-periodic
        # Stokes drag — off by default (backward compatible)
        xi_drag_per_p    = tf.zeros([P], dtype=DTYPE),        # (P,) drag per arc length
        # Per-particle internal elastic damping (∝ R0^{-1}); 0 = use scalar fallback
        alpha_damp_per_p = tf.zeros([P], dtype=DTYPE),        # (P,)
        U_bg_type        = tf.constant(0, dtype=tf.int32),    # 0=zero preset
        U_bg_params      = tf.zeros([4], dtype=DTYPE),        # preset parameters
    )
    return state, params


# ── driven-particle helpers ───────────────────────────────────────────────────

def make_traj(v_dc=(0.0, 0.0), v_ac=(0.0, 0.0), freq=(0.0, 0.0), phase=(0.0, 0.0),
              omega_spin_dc=0.0, omega_spin_ac=0.0, omega_spin_freq=0.0, omega_spin_phase=0.0,
              omega_orbit_dc=0.0, omega_orbit_ac=0.0, omega_orbit_freq=0.0, omega_orbit_phase=0.0,
              r_ref=(0.0, 0.0)):
    """
    Build a (18,) trajectory parameter array for one driven particle.

    Velocity components: v(t) = dc + ac * cos(freq * t + phase)
    - v_dc, v_ac, freq, phase : (x, y) tuples
    - omega_spin_*  : self-rotation about own CM
    - omega_orbit_* : orbital rotation about r_ref
    - r_ref         : (x, y) reference point for orbital motion

    Example — constant rightward:
        make_traj(v_dc=(1.0, 0.0))

    Example — oscillatory shear:
        make_traj(v_ac=(A, 0), freq=(omega, 0))

    Example — Couette orbit at angular velocity Omega:
        make_traj(omega_orbit_dc=Omega, r_ref=(cx, cy))
    """
    return np.array([
        v_dc[0],        v_ac[0],        freq[0],        phase[0],
        v_dc[1],        v_ac[1],        freq[1],        phase[1],
        omega_spin_dc,  omega_spin_ac,  omega_spin_freq, omega_spin_phase,
        omega_orbit_dc, omega_orbit_ac, omega_orbit_freq, omega_orbit_phase,
        r_ref[0],       r_ref[1],
    ], dtype=NP_DTYPE)


def set_driven(params, indices, traj_rows, frozen=False):
    """
    Mark particles as driven and set their trajectories.

    Modifies params dicts in place (replaces driven_mask, shape_frozen, traj tensors).

    Parameters
    ----------
    params     : dict from make_state() — will be mutated
    indices    : list of int particle indices to drive
    traj_rows  : list of (18,) arrays from make_traj(), one per index
    frozen     : bool or list of bool — freeze elastic shape for each driven particle
    """
    P = params['driven_mask'].shape[0]

    dm  = params['driven_mask'].numpy().copy()
    sf  = params['shape_frozen'].numpy().copy()
    tr  = params['traj'].numpy().copy()

    if isinstance(frozen, bool):
        frozen = [frozen] * len(indices)

    for i, (idx, trow, frz) in enumerate(zip(indices, traj_rows, frozen)):
        dm[idx]     = 1.0
        sf[idx]     = 1.0 if frz else 0.0
        tr[idx, :]  = np.asarray(trow, dtype=NP_DTYPE)

    params['driven_mask']  = tf.constant(dm)
    params['shape_frozen'] = tf.constant(sf)
    params['traj']         = tf.constant(tr)


def set_periodic_box(params, Lx, Ly, periodic_x=True, periodic_y=True):
    """Set the periodic box dimensions in params (activates periodic BC in forces).

    For mixed-periodicity setups (e.g. periodic_x=True with hard top/bottom walls
    in a Couette-like geometry), pass periodic_x / periodic_y individually.
    Non-periodic axes are encoded as 0 in the box tensor; the force kernel
    skips min-image on those axes (open / hard-wall behavior).
    """
    bx = NP_DTYPE(Lx) if periodic_x else NP_DTYPE(0.0)
    by = NP_DTYPE(Ly) if periodic_y else NP_DTYPE(0.0)
    params['box'] = tf.constant([bx, by])


def state_to_numpy(state):
    return {k: v.numpy() for k, v in state.items()}


def update_particles_from_state(particles, state):
    """Write TF state back into CapsuleParticle objects for frozen-branch interop."""
    s = state_to_numpy(state)
    for i, p in enumerate(particles):
        p.x     = s['x_all'][i].astype(np.float64)
        p.x_cm  = s['x_cm'][i].astype(np.float64)
        p.v_cm  = s['v_cm'][i].astype(np.float64)
        p.theta = float(s['theta'][i])
        p.omega = float(s['omega'][i])
        p.u     = s['u'][i].astype(np.float64)
        p.u_dot = s['u_dot'][i].astype(np.float64)


# ── tf.while_loop simulation runner ──────────────────────────────────────────

def run_simulation_tf(state0, dt, alpha_damp, g, params, n_steps, mgr_cpp,
                      skin, prim_data=None, R0_max=1.0,
                      cand_check_interval=10, diagnostics=False,
                      step_offset=0, use_tf_function=False):
    """
    Run n_steps using tf.while_loop with C++ candidacy via tf.py_function.

    prim_data is treated as permanently static — oscillating walls encode their
    motion as (seg_osc_A, seg_osc_omega, seg_osc_sign) parameters computed
    inside the TF graph from the current time, so prim_data never needs rebuilding.

    Parameters
    ----------
    state0              : dict — initial state from make_state()
    dt                  : float
    alpha_damp          : float
    g                   : float
    params              : dict — from make_state()
    n_steps             : int
    mgr_cpp             : CandidacyManager
    skin                : float
    prim_data           : dict or None — built ONCE via make_prim_data(); never rebuilt
    R0_max              : float
    cand_check_interval : int — call tf.py_function every this many steps (default 10)
                          Reduces Python callbacks from n_steps to n_steps/interval.
    diagnostics         : bool — if True, return (final_state, diag_dict)
    step_offset         : int — global step index of the first step in this call.
                          Used when running in chunks: pass chunk_i * steps_per_chunk
                          so that t = (step_offset + local_idx) * dt is correct
                          for moving primitives across chunk boundaries.
    use_tf_function     : bool — if True, wrap the tf.while_loop in a tf.function
                          so the entire loop runs as a single graph (eliminates
                          per-iteration eager-dispatch overhead). Costs one ~5–10s
                          trace on the first call per call-site. Default False
                          preserves the original eager dispatch behaviour.

    Returns
    -------
    final_state              — if diagnostics=False
    (final_state, diag_dict) — if diagnostics=True
      diag_dict keys: n_cand_checks, n_cand_updates, n_retraces_step_full, n_steps
    """
    dtype = DTYPE

    # Initial candidacy
    xc0 = np.ascontiguousarray(state0['x_cm'].numpy(), dtype=np.float64)
    th0 = np.ascontiguousarray(state0['theta'].numpy(), dtype=np.float64)
    xa0 = np.ascontiguousarray(state0['x_all'].numpy(), dtype=np.float64)
    mgr_cpp.update(xc0, th0, x_all=xa0)

    K = mgr_cpp.P * mgr_cpp.N
    E = mgr_cpp.E
    CapCand_init = tf.constant(mgr_cpp.CapCandidates)

    # Use stable tensor objects stored in params so that loop_body closures share
    # the SAME tensor references across chunks → step_full_tf is traced only once.
    dt_tf    = params.get('_dt_tf',    tf.constant(NP_DTYPE(dt)))
    alpha_tf = params.get('_alpha_tf', tf.constant(NP_DTYPE(alpha_damp)))
    g_tf     = params.get('_g_tf',     tf.constant(NP_DTYPE(g)))

    skin_thresh = float(skin * 0.5)

    last_xc       = [xc0.copy()]
    last_th       = [th0.copy()]
    n_cand_checks  = [0]
    n_cand_updates = [0]

    # Managers can opt out of tf_sim's skin gating by exposing `always_update`
    # = True (PRCM does this — it has its own per-particle internal triggers
    # and cascade, so it must be called every candidacy_step to keep L1/L2
    # current. Production CM keeps the skin gating to avoid expensive rebuilds
    # when nothing has moved).
    mgr_always_update = bool(getattr(mgr_cpp, 'always_update', False))

    def candidacy_step(x_cm, theta, x_all):
        """Python callback: update C++ candidacy.

        x_all is plumbed through for managers that need per-node positions
        (e.g., PRCM). Production CM ignores it.
        """
        n_cand_checks[0] += 1
        xc = x_cm.numpy()
        th = theta.numpy()
        xa = x_all.numpy()
        if mgr_always_update:
            should_update = True
        else:
            dx  = float(np.max(np.linalg.norm(xc - last_xc[0], axis=1)))
            dth = float(np.max(np.abs(th - last_th[0])))
            should_update = (dx + R0_max * dth) > skin_thresh
        if should_update:
            mgr_cpp.update(np.ascontiguousarray(xc, dtype=np.float64),
                           np.ascontiguousarray(th, dtype=np.float64),
                           x_all=np.ascontiguousarray(xa, dtype=np.float64))
            last_xc[0] = xc.copy()
            last_th[0] = th.copy()
            n_cand_updates[0] += 1
        return mgr_cpp.CapCandidates.copy()

    cand_interval_tf = tf.constant(cand_check_interval, dtype=tf.int32)
    step_end         = tf.constant(step_offset + n_steps, dtype=tf.int32)

    def loop_body(state, CapCand, step_idx, max_closing_ratio):
        # Call Python candidacy check every cand_check_interval steps only.
        new_cand = tf.cond(
            tf.equal(step_idx % cand_interval_tf, 0),
            true_fn=lambda: tf.ensure_shape(
                tf.py_function(candidacy_step,
                               [state['x_cm'], state['theta'], state['x_all']],
                               Tout=tf.int32),
                [K, E]),
            false_fn=lambda: CapCand,
        )

        # t is global simulation time — correct across chunk boundaries
        t = tf.cast(step_idx, dtype) * dt_tf

        new_state, metrics = step_full_tf(
            state, new_cand, dt_tf, alpha_tf, g_tf, params,
            t=t, prim_data=prim_data)

        # Track running max of contact-closing ratio across the chunk
        # for the System-level stability watchdog (covers both particle-particle
        # and particle-wall contacts; symmetric in who is moving).
        new_max = tf.maximum(max_closing_ratio, metrics['max_closing_ratio'])
        return new_state, new_cand, step_idx + 1, new_max

    def loop_cond(state, CapCand, step_idx, max_closing_ratio):
        return step_idx < step_end

    step_init = tf.constant(step_offset, dtype=tf.int32)
    max_init  = tf.constant(0.0, dtype=dtype)

    def _exec_loop():
        return tf.while_loop(
            loop_cond,
            loop_body,
            loop_vars=(state0, CapCand_init, step_init, max_init),
            parallel_iterations=1)

    if use_tf_function:
        # Graph-compile the while_loop so all n_steps iterations dispatch
        # as a single graph execution (no per-iteration Python/runtime tick).
        # tf.py_function inside the loop body still works in graph mode
        # (it's only XLA jit_compile=True that forbids py_function).
        # Pays a ~5–10s trace cost on first call per closure signature.
        final_state, _, _, max_closing_ratio_chunk = tf.function(_exec_loop)()
    else:
        final_state, _, _, max_closing_ratio_chunk = _exec_loop()

    if diagnostics:
        diag = dict(
            n_cand_checks   = n_cand_checks[0],
            n_cand_updates  = n_cand_updates[0],
            n_retraces_step_full = step_full_tf.experimental_get_tracing_count(),
            n_steps         = n_steps,
            max_closing_ratio_chunk = float(max_closing_ratio_chunk.numpy()),
        )
        return final_state, diag
    return final_state
