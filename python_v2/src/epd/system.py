"""
system.py — High-level EPD simulation orchestrator.

Typical usage
-------------
    from src.epd.system import System
    from src.epd.particles import ParticleSpec

    sys = System(Lx=12, Ly=12)
    sys.add_particles(ParticleSpec(count=20, nu=0.5, poly_dist=0.05))
    sys.initialize(phi_target=0.80)
    sys.run(1000, record_every=50, output_dir='results/run1/')
    sys.save('results/run1/checkpoint/')
"""

import os
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class System:
    """
    High-level EPD simulation system.

    Parameters
    ----------
    Lx, Ly          : float — simulation box dimensions
    periodic_x/y    : bool  — periodic boundary conditions (default True)
    dt_factor       : float — multiplier on auto-stable dt (default 0.4)
    alpha_damp      : float — damping coefficient; None → auto from first spec
    E_candidates    : int   — candidacy table width (default 128)
    skin            : float — candidacy skin factor (multiplies mean R0; default 1.0)
    """

    def __init__(self,
                 Lx, Ly,
                 periodic_x=True, periodic_y=True,
                 dt_factor=1.5,
                 alpha_damp=None,
                 E_candidates=32,
                 skin=1.0,
                 g=0.0,
                 # Stability watchdog: per-step max nodal displacement / r_c thresholds.
                 # Crossing these fires warnings during sys.run(). Configurable for users
                 # who know their problem; defaults are conservative.
                 disp_advisory=0.05,
                 disp_strong=0.10,
                 disp_critical=0.20,
                 use_tf_function=False,
                 candidacy_kind='production'):
        self.Lx = float(Lx)
        self.Ly = float(Ly)
        self._config = {
            'periodic_x': bool(periodic_x),
            'periodic_y': bool(periodic_y),
        }
        self._dt_factor     = float(dt_factor) if dt_factor is not None else 0.4
        self._alpha_damp    = float(alpha_damp) if alpha_damp is not None else None
        self._E_candidates  = int(E_candidates)
        self._skin          = float(skin)
        self._g_val         = float(g)
        self._candidacy_kind = str(candidacy_kind)
        if self._candidacy_kind not in ('production', 'prcm'):
            raise ValueError(f"candidacy_kind must be 'production' or 'prcm', "
                              f"got {self._candidacy_kind!r}")
        self._disp_advisory = float(disp_advisory)
        self._disp_strong   = float(disp_strong)
        self._disp_critical = float(disp_critical)
        self._disp_advisory_fired = False     # one-shot per run.clear_recording()
        self._max_disp_ratio_run  = 0.0       # running max for end-of-run summary
        # Default dispatch mode for run/run_fast. When True, the inner
        # tf.while_loop runs as a single graph (one ~5–10s trace per call site,
        # then no per-iteration Python/runtime tick). Set once here or via the
        # `use_tf_function` property; per-call kwarg overrides this default.
        self._use_tf_function = bool(use_tf_function)

        # Populated by add_particles / add_object
        self._specs    = []    # list[ParticleSpec]
        self._objects  = []    # list[SimulationObject]

        # Set after initialize()
        self._particles  = []
        self._state      = None
        self._params     = None
        self._cm_mgr     = None
        self._prim_data  = None
        self._dt         = None
        self._dt_max     = None     # CFL upper bound on dt (set during initialize)
        self._dt_tf      = None
        self._alpha_tf   = None
        self._g_tf       = None

        # Background flow (None = quiescent; set via self.U_background = ...)
        self._U_background    = None

        # Recording — populated by run(); wiped by clear_recording()
        self.frames           = []   # snapshot dicts from snapshot()
        self.callback_data    = []   # callback() return dicts
        self.diag             = []   # per-chunk plumbing dicts
        self._particle_colors = []   # auto-assigned or user-set
        self._pending_palette = None # (name, seed) deferred from set_color_palette

    # ── dt / dt_factor sync API ───────────────────────────────────────────────
    # Two equivalent ways to control the integration step.  Setting either
    # property updates the other and propagates through every TF tensor that
    # carries the value (self._dt_tf, params['_dt_tf']).
    #
    #   sys.dt_factor = 1.5   # multiplicative scale on the CFL bound dt_max
    #   sys.dt        = 0.002 # absolute dt; dt_factor = dt / dt_max recomputed
    #
    # Before initialize() runs, dt_max is unknown; assignments are stored and
    # applied at initialize time.
    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        v = float(value)
        self._dt = v
        if self._dt_max is not None and self._dt_max > 0:
            self._dt_factor = v / self._dt_max
        self._sync_dt_tensors()

    @property
    def dt_factor(self):
        return self._dt_factor

    @dt_factor.setter
    def dt_factor(self, value):
        v = float(value)
        self._dt_factor = v
        if self._dt_max is not None and self._dt_max > 0:
            self._dt = v * self._dt_max
        self._sync_dt_tensors()

    @property
    def use_tf_function(self):
        """Default dispatch mode for run/run_fast. When True, the inner
        tf.while_loop is wrapped in tf.function and executes as a single
        graph (eliminates per-iteration eager-dispatch overhead). Per-call
        kwarg on run/run_fast overrides this default when explicitly set."""
        return self._use_tf_function

    @use_tf_function.setter
    def use_tf_function(self, value):
        self._use_tf_function = bool(value)

    def _sync_dt_tensors(self):
        """Propagate self._dt to TF tensors used by run_simulation_tf and step_full_tf."""
        if self._dt is None:
            return
        import tensorflow as tf
        from src.simulation.tf_sim import NP_DTYPE, DTYPE
        self._dt_tf = tf.constant(NP_DTYPE(self._dt), dtype=DTYPE)
        if self._params is not None:
            self._params['_dt_tf'] = self._dt_tf

    # ── registration ──────────────────────────────────────────────────────────

    def add_particles(self, spec):
        """Register a ParticleSpec. Returns self for chaining."""
        spec._system = self   # back-reference so setters can push live params
        self._specs.append(spec)
        return self

    def add_object(self, obj):
        """Register a SimulationObject (wall, obstacle, etc.). Returns self."""
        self._objects.append(obj)
        return self

    # ── recording helpers ─────────────────────────────────────────────────────

    def clear_recording(self):
        """Wipe frames/callback_data/diag without touching simulation state."""
        self.frames        = []
        self.callback_data = []
        self.diag          = []
        # Reset the per-run stability watchdog so a fresh run gets fresh warnings
        self._disp_advisory_fired = False
        self._max_disp_ratio_run  = 0.0
        return self

    def _build_cm(self, P, N_val, R0_arr, particles):
        """Construct the candidacy manager — production sliding-fill or PRCM
        depending on self._candidacy_kind.
        """
        from src.simulation.candidacy_manager import CandidacyManager
        if self._candidacy_kind == 'prcm':
            from src.simulation.pair_registration_cm import PairRegistrationCM
            r_c_arr = np.array([p.r_c for p in particles])
            return PairRegistrationCM(
                P=P, N=N_val, R0_arr=R0_arr, r_c_per_p=r_c_arr,
                Lx=self.Lx, Ly=self.Ly,
                M1=20, M2=10,
                delta=(self._E_candidates - 1) // 2,
                periodic_x=self._config['periodic_x'],
                periodic_y=self._config['periodic_y'])
        return CandidacyManager(
            P=P, N=N_val,
            R0=float(np.mean(R0_arr)),
            E=self._E_candidates,
            skin=self._skin * float(np.mean(R0_arr)),
            periodic=self._config['periodic_x'] and self._config['periodic_y'],
            periodic_x=self._config['periodic_x'],
            periodic_y=self._config['periodic_y'],
            Lx=self.Lx, Ly=self.Ly,
            R0_arr=R0_arr,
        )

    def _assign_particle_colors(self):
        """Auto-assign per-particle colors from palette (or p.color if set).

        Honors user-set palettes:
          - If a palette name was set via set_color_palette() before particles
            were realized, apply it now.
          - If _particle_colors is already populated (matches len(_particles)),
            assume the user assigned them and skip the default.
        """
        if getattr(self, '_pending_palette', None) is not None:
            from src.epd.palettes import apply_palette
            name, seed = self._pending_palette
            apply_palette(self, name, seed=seed)
            self._pending_palette = None
            return
        if len(self._particle_colors) == len(self._particles):
            return
        from src.simulation.emulsion_particle import EmulsionParticle
        palette_e  = ['cornflowerblue', 'salmon', 'mediumseagreen', 'orchid',
                      'goldenrod', 'lightcoral', 'steelblue', 'sandybrown']
        palette_em = ['steelblue', 'tomato', 'mediumturquoise', 'mediumpurple',
                      'coral', 'mediumaquamarine']
        self._particle_colors = []
        ei = 0; emi = 0
        for p in self._particles:
            if hasattr(p, 'color') and p.color is not None:
                self._particle_colors.append(p.color)
            elif isinstance(p, EmulsionParticle):
                self._particle_colors.append(palette_em[emi % len(palette_em)])
                emi += 1
            else:
                self._particle_colors.append(palette_e[ei % len(palette_e)])
                ei += 1

    def set_color_palette(self, name, seed=None):
        """Randomly recolor every particle from the named palette in
        src.epd.palettes.PALETTES (e.g. 'palette1' through 'palette12').

        Safe to call anywhere in the build sequence:
          - Before initialize(): defers; palette is applied automatically
            once particles are realized.
          - After initialize(): applies immediately.
        Either way, the palette survives any later auto-color call inside
        initialize(). Per-particle overrides (e.g. driven-ring tinting in
        set_driven_particles handlers) should still come last."""
        from src.epd.palettes import apply_palette, PALETTES
        if name not in PALETTES:
            raise ValueError(f"unknown palette {name!r}; "
                              f"available: {sorted(PALETTES)}")
        if len(self._particles) == 0:
            self._pending_palette = (name, seed)
            return None
        return apply_palette(self, name, seed=seed)

    # ── initialization ────────────────────────────────────────────────────────

    def initialize(self, phi_target=0.80, seed=42, verbose=True,
                   relax_only=False, n_relax_init=200,
                   cand_check_interval=10, **swell_kwargs):
        """
        RSA seed + adaptive swell to phi_target.
        Sets state['t']=0 and state['step']=0 after completion.

        cand_check_interval : int — Python candidacy callback every N steps
                              during initial relax / adaptive_swell. Default 10.
                              Larger values reduce CPU↔GPU sync count (helpful
                              on virtualized PCIe like Colab) at the cost of
                              looser candidacy freshness within each block.

        Inherits self.use_tf_function for the inner relax loops in both
        branches (relax_only and adaptive_swell).

        Returns self.
        """
        import tensorflow as tf
        import src.simulation.tf_sim as tf_sim_mod
        tf_sim_mod.set_dtype(tf.float64)
        from src.simulation.tf_sim import (make_state, set_periodic_box,
                                           make_prim_data, DTYPE, NP_DTYPE)
        from src.simulation.candidacy_manager import CandidacyManager
        from src.epd.initializer import (rsa_seed, compute_phi_outer,
                                         adaptive_swell)

        # 1. Compute RSA initial box
        #    Sample sizes to know total area; use phi_rsa_init=0.25 for placement.
        phi_rsa_init = 0.25
        _rng_tmp = np.random.default_rng(seed)
        _R0_preview = []
        for sp in self._specs:
            _R0_preview.extend(sp.sample_R0(seed=int(_rng_tmp.integers(0, 2**32))))
        _N_nodes_preview = [sp.N_nodes for sp in self._specs
                            for _ in range(sp.count)]
        _R_eff = [_R0_preview[i] + 2 * _R0_preview[i] * np.sin(np.pi / _N_nodes_preview[i])
                  for i in range(len(_R0_preview))]
        _A_total = np.sum(np.pi * np.array(_R_eff)**2)
        if relax_only:
            # For relax-only (no swell), never auto-expand: RSA draws from
            # the container bounding box and rejects out-of-region samples.
            Lx_rsa, Ly_rsa = self.Lx, self.Ly
        else:
            phi_at_Lx = _A_total / (self.Lx * self.Ly)
            if phi_at_Lx > phi_rsa_init:
                # box too small for RSA; auto-expand preserving aspect ratio
                aspect = self.Lx / self.Ly if self.Ly > 0 else 1.0
                A_rsa  = _A_total / phi_rsa_init
                Lx_rsa = np.sqrt(A_rsa * aspect)
                Ly_rsa = A_rsa / Lx_rsa
                if verbose:
                    print(f"  [System] auto-expanding RSA box: "
                          f"{self.Lx:.2f}×{self.Ly:.2f} → {Lx_rsa:.2f}×{Ly_rsa:.2f} "
                          f"(phi@Lx={phi_at_Lx:.3f} > phi_rsa={phi_rsa_init})")
            else:
                Lx_rsa, Ly_rsa = self.Lx, self.Ly

        # 1b. RSA placement in initial box
        particles, positions = rsa_seed(
            self._specs, self._objects, Lx_rsa, Ly_rsa,
            seed=seed, verbose=verbose)
        self._particles = particles
        self.Lx = Lx_rsa
        self.Ly = Ly_rsa

        # 2. Build TF state + params
        state, params = make_state(particles)

        if self._config['periodic_x'] or self._config['periodic_y']:
            set_periodic_box(params, self.Lx, self.Ly,
                              periodic_x=self._config['periodic_x'],
                              periodic_y=self._config['periodic_y'])

        # 2b. Wire per-particle xi_drag and alpha_damp from particle attributes.
        #     particles[i].xi_drag and .alpha_damp were set by ParticleSpec.build()
        #     with correct R0-scaling (xi ∝ R0^{+1}, alpha ∝ R0^{-1}).
        xi_arr    = np.array([p.xi_drag    for p in particles], dtype=np.float64)
        alpha_arr = np.array([p.alpha_damp for p in particles], dtype=np.float64)
        params['xi_drag_per_p']    = tf.constant(xi_arr,    dtype=DTYPE)
        params['alpha_damp_per_p'] = tf.constant(alpha_arr, dtype=DTYPE)

        # 3. Set shape_frozen from specs
        _frozen = np.zeros(len(particles))
        idx = 0
        for sp in self._specs:
            if sp.frozen_shape:
                _frozen[idx:idx+sp.count] = 1.0
            idx += sp.count
        params['shape_frozen'] = tf.constant(_frozen, dtype=DTYPE)

        # 4. Determine dt and alpha
        # For polydisperse systems: use the smallest-R0 particle (highest frequency,
        # most constraining dt). T_wave ∝ R0^(3/2), so smallest R0 → smallest T_wave.
        R0_arr_all = np.array([p.R0 for p in particles])
        R0_mean    = float(np.mean(R0_arr_all))
        if verbose and len(particles) > 1:
            R0_std = float(np.std(R0_arr_all))
            print(f"  [System] <R0>={R0_mean:.6f}  std={R0_std:.4f}  "
                  f"min={R0_arr_all.min():.4f}  max={R0_arr_all.max():.4f}")
        assert abs(R0_mean - 1.0) < 0.01, \
            f"<R0>={R0_mean:.6f} deviates from 1.0 by more than 1% — " \
            f"check ParticleSpec.sample_R0() normalisation"

        p0 = min(particles, key=lambda p: p.R0)   # most constrained (smallest R0)
        if p0.El_t > 0:
            # Elastic capsule: use edge-spring wave speed
            om_edge   = p0.N * np.sqrt(p0.El_t / (p0.rho_d * p0.L0)) / p0.R0
            T_wave    = 2.0 * np.pi * p0.R0 / np.sqrt(p0.El_t * p0.R0 / (p0.rho_d * p0.L0))
        else:
            # Emulsion particle: use surface-tension capillary wave speed
            gamma_p   = float(getattr(p0, 'gamma', 1.0))
            c_cap     = np.sqrt(gamma_p / (p0.rho_d * p0.R0))
            om_edge   = c_cap * p0.N / p0.R0
            T_wave    = 2.0 * np.pi / np.sqrt(6.0) / np.sqrt(gamma_p / (p0.rho_d * p0.R0**3))

        # Constrain dt by contact stiffness (prevents instability for soft/emulsion particles)
        if p0.m_node > 0:
            om_contact = np.sqrt(p0.k_c / p0.m_node)
            om_max     = max(om_edge, om_contact)
        else:
            om_max = om_edge
        dt_max_cfl = float(2.0 / om_max)            # CFL upper bound (factor=1.0)
        dt_val     = float(self._dt_factor * dt_max_cfl)

        if self._alpha_damp is not None:
            alpha_val = self._alpha_damp
        else:
            alpha_val = 2.0 / T_wave

        self._dt        = dt_val
        self._dt_max    = dt_max_cfl                 # remembered for dt/dt_factor sync API
        self._alpha_damp = alpha_val

        # Membrane wave speed reference (for ringing metric in SimMonitor)
        v_mem = R0_mean * (2.0 * np.pi / T_wave)   # R0_mean × ω_mem

        params['_dt_tf']    = tf.constant(NP_DTYPE(dt_val))
        params['_alpha_tf'] = tf.constant(NP_DTYPE(alpha_val))
        params['_g_tf']     = tf.constant(NP_DTYPE(self._g_val))   # stable ref for retrace fix
        params['_v_mem']    = tf.constant(NP_DTYPE(v_mem))

        # 5. Candidacy manager (production sliding-fill or PRCM, opt-in)
        P      = len(particles)
        N_val  = particles[0].N
        R0_arr = np.array([p.R0 for p in particles])
        cm_mgr = self._build_cm(P, N_val, R0_arr, particles)
        cm_mgr.update(state['x_cm'].numpy(), state['theta'].numpy(),
                       x_all=state['x_all'].numpy())

        # 6. Build prim_data
        t_now    = 0.0
        prim_list = []
        for obj in self._objects:
            prim_list.extend(obj.to_make_prim_list(t=t_now))
        prim_data = make_prim_data(prim_list)

        # 7. Adaptive swell (or relax-only for fixed-wall geometries)
        r_c_arr = params['r_c_per_p'].numpy()
        if relax_only:
            # No compression: settle in-place using the same tf.while_loop path
            # as run_fast — one compiled op, not 2000 individual Python calls.
            from src.simulation.tf_sim import run_simulation_tf, DTYPE as _DTYPE, NP_DTYPE as _NPDTYPE
            if verbose:
                print(f"  [System] relax_only: running {n_relax_init} settle steps ...")
            if n_relax_init > 0:
                skin_abs    = self._skin * float(np.mean(R0_arr))
                R0_max_val  = float(np.max(R0_arr))
                phys_state  = {k: v for k, v in state.items() if k not in ('t', 'step')}
                new_phys = run_simulation_tf(
                    phys_state,
                    float(params['_dt_tf']),
                    float(params['_alpha_tf']),
                    0.0,                      # g=0 during relax
                    params, n_relax_init, cm_mgr,
                    skin=skin_abs,
                    prim_data=prim_data,
                    R0_max=R0_max_val,
                    cand_check_interval=cand_check_interval,
                    diagnostics=False,
                    step_offset=0,
                    use_tf_function=self._use_tf_function,
                )
                state = {
                    **new_phys,
                    't':    tf.constant(_NPDTYPE(0.0), dtype=_DTYPE),
                    'step': tf.constant(0, dtype=tf.int64),
                }
        else:
            # Convert user's phi_target (fraction of accessible area) to the
            # box-fraction used internally by adaptive_swell (which divides by
            # Lx*Ly).  The ratio acc/box stays constant throughout swell, so a
            # single pre-swell correction is exact.
            acc_area = self._accessible_area()
            box_area = self.Lx * self.Ly
            phi_target_swell = phi_target * (acc_area / box_area)
            skin_abs   = self._skin * float(np.mean(R0_arr_all))
            R0_max_val = float(np.max(R0_arr_all))
            state, Lx_new, Ly_new = adaptive_swell(
                state, params, cm_mgr, prim_data,
                phi_target=phi_target_swell,
                Lx=self.Lx, Ly=self.Ly,
                verbose=verbose,
                objects=self._objects,
                skin=skin_abs, R0_max=R0_max_val,
                cand_check_interval=cand_check_interval,
                use_tf_function=self._use_tf_function,
                **swell_kwargs,
            )
            self.Lx = Lx_new
            self.Ly = Ly_new

        # Rebuild prim_data from scaled objects after swell
        prim_list = []
        for obj in self._objects:
            prim_list.extend(obj.to_make_prim_list(t=0.0))
        prim_data = make_prim_data(prim_list)

        # Update box in params after swell
        if self._config['periodic_x'] or self._config['periodic_y']:
            set_periodic_box(params, self.Lx, self.Ly,
                              periodic_x=self._config['periodic_x'],
                              periodic_y=self._config['periodic_y'])

        # 8. Wire driven particles from ParticleSpec.motion
        from src.simulation.tf_sim import set_driven
        _driven_idxs = []
        _traj_rows   = []
        _frozen_flags = []
        particle_offset = 0
        for sp in self._specs:
            if sp.motion is not None:
                for k in range(sp.count):
                    _driven_idxs.append(particle_offset + k)
                    _traj_rows.append(sp.motion.to_traj_row())
                    _frozen_flags.append(sp.frozen_shape)
            particle_offset += sp.count
        if _driven_idxs:
            set_driven(params, _driven_idxs, _traj_rows, frozen=_frozen_flags)

        # 9. Finalize state
        state['t']    = tf.constant(0.0, dtype=DTYPE)
        state['step'] = tf.constant(0,   dtype=tf.int64)

        self._state    = state
        self._params   = params
        self._cm_mgr   = cm_mgr
        self._prim_data = prim_data

        # Store TF scalars
        self._dt_tf    = tf.constant(NP_DTYPE(dt_val))
        self._alpha_tf = tf.constant(NP_DTYPE(alpha_val))
        self._g_tf     = tf.constant(NP_DTYPE(self._g_val))

        # Rebuild candidacy at final positions
        cm_mgr.update(state['x_cm'].numpy(), state['theta'].numpy(),
                       x_all=state['x_all'].numpy())

        if verbose:
            phi = compute_phi_outer(state, self.Lx, self.Ly, r_c_arr,
                                    accessible_area=self._accessible_area())
            print(f"  [System] initialized: P={P}, phi_outer={phi:.4f}, "
                  f"dt={dt_val:.4e}, alpha={alpha_val:.4f}")

        self._assign_particle_colors()
        self.clear_recording()
        return self

    def initialize_from_particles(self, particles, alpha_damp=None, dt=None):
        """
        Initialize System from pre-built, pre-positioned particle objects.
        Bypasses RSA seeding and adaptive swell — use when exact initial
        geometry is required (e.g. Test A two-particle squeeze).

        particles  : list of CapsuleParticle / EmulsionParticle already at
                     their desired positions
        alpha_damp : override damping coefficient (None → auto from T_wave)
        dt         : override time step (None → auto from stability estimate)

        Returns self.
        """
        import tensorflow as tf
        import src.simulation.tf_sim as tf_sim_mod
        tf_sim_mod.set_dtype(tf.float64)
        from src.simulation.tf_sim import (make_state, make_prim_data,
                                           DTYPE, NP_DTYPE)
        from src.simulation.candidacy_manager import CandidacyManager

        self._particles = list(particles)
        P  = len(particles)
        p0 = particles[0]

        # 1. TF state + params
        state, params = make_state(particles)

        # 2. dt — prefer explicit; otherwise derive from wave-speed stability
        if dt is not None:
            dt_val = float(dt)
        else:
            if p0.El_t > 0:
                om_edge = p0.N * np.sqrt(p0.El_t / (p0.rho_d * p0.L0)) / p0.R0
            else:
                gamma_p = float(getattr(p0, 'gamma', 1.0))
                om_edge = np.sqrt(gamma_p / (p0.rho_d * p0.R0)) * p0.N / p0.R0
            om_max = max(om_edge,
                         np.sqrt(p0.k_c / p0.m_node) if p0.m_node > 0 else 0.0)
            dt_max_cfl = float(2.0 / om_max)
            dt_val = float(self._dt_factor * dt_max_cfl)
            self._dt_max = dt_max_cfl

        # 3. alpha_damp
        if alpha_damp is not None:
            alpha_val = float(alpha_damp)
        elif self._alpha_damp is not None:
            alpha_val = float(self._alpha_damp)
        else:
            if p0.El_t > 0:
                T_w = 2.0 * np.pi * p0.R0 / np.sqrt(p0.El_t * p0.R0 /
                                                       (p0.rho_d * p0.L0))
            else:
                gamma_p = float(getattr(p0, 'gamma', 1.0))
                T_w = 2.0 * np.pi / np.sqrt(6.0 * gamma_p /
                                              (p0.rho_d * p0.R0**3))
            alpha_val = 2.0 / T_w

        self._dt         = dt_val
        self._alpha_damp = alpha_val

        # 4. CandidacyManager — identical parameters to TF test defaults
        R0_arr = np.array([p.R0 for p in particles])
        N_val  = p0.N
        cm_mgr = self._build_cm(P, N_val, R0_arr, particles)
        cm_mgr.update(state['x_cm'].numpy(), state['theta'].numpy(),
                       x_all=state['x_all'].numpy())

        # 5. prim_data at t=0 from registered objects
        prim_list = []
        for obj in self._objects:
            prim_list.extend(obj.to_make_prim_list(t=0.0))
        prim_data = make_prim_data(prim_list)

        # 6. Finalise state and wire TF scalars
        state['t']    = tf.constant(0.0, dtype=DTYPE)
        state['step'] = tf.constant(0,   dtype=tf.int64)

        self._state     = state
        self._params    = params
        self._cm_mgr    = cm_mgr
        self._prim_data = prim_data

        self._dt_tf    = tf.constant(NP_DTYPE(dt_val),    dtype=DTYPE)
        self._alpha_tf = tf.constant(NP_DTYPE(alpha_val), dtype=DTYPE)
        self._g_tf     = tf.constant(NP_DTYPE(self._g_val), dtype=DTYPE)

        # Store stable tensor refs in params for retrace fix
        params['_dt_tf']    = self._dt_tf
        params['_alpha_tf'] = self._alpha_tf
        params['_g_tf']     = self._g_tf

        # Per-particle alpha_damp and xi_drag (use particle attrs if set, else zeros)
        alpha_arr = np.array([p.alpha_damp for p in particles], dtype=np.float64)
        xi_arr    = np.array([p.xi_drag    for p in particles], dtype=np.float64)
        params['alpha_damp_per_p'] = tf.constant(alpha_arr, dtype=DTYPE)
        params['xi_drag_per_p']    = tf.constant(xi_arr,    dtype=DTYPE)

        self._assign_particle_colors()
        self.clear_recording()
        return self

    # ── parameter update API ─────────────────────────────────────────────────

    @property
    def U_background(self):
        """Background flow preset. None = quiescent."""
        return self._U_background

    @U_background.setter
    def U_background(self, value):
        """
        Set background flow preset.

        value : None
                ('zero',)
                ('constant',    {'U': (Ux, Uy)})
                ('shear',       {'rate': gamma_dot})
                ('parabolic',   {'U_max': Um, 'H': H})
                ('extensional', {'rate': eps_dot})
        """
        import tensorflow as tf
        from src.simulation.tf_sim import DTYPE

        _presets = {'zero': 0, 'constant': 1, 'shear': 2,
                    'parabolic': 3, 'extensional': 4, 'poiseuille_v': 5}

        if value is None or value == 'zero' or value == ('zero',):
            bg_type   = 0
            bg_params = np.zeros(4)
        else:
            name, kw = value[0], value[1] if len(value) > 1 else {}
            if name not in _presets:
                raise ValueError(f"Unknown U_background preset '{name}'. "
                                 f"Choose from: {list(_presets.keys())}")
            bg_type   = _presets[name]
            bg_params = np.zeros(4)
            if name == 'constant':
                U = kw.get('U', (0.0, 0.0))
                bg_params[0], bg_params[1] = float(U[0]), float(U[1])
            elif name == 'shear':
                bg_params[0] = float(kw.get('rate', 0.0))
            elif name == 'parabolic':
                bg_params[0] = float(kw.get('U_max', 0.0))
                bg_params[1] = float(kw.get('H', 1.0))
            elif name == 'extensional':
                bg_params[0] = float(kw.get('rate', 0.0))
            elif name == 'poiseuille_v':
                bg_params[0] = float(kw.get('U_max', 0.0))
                bg_params[1] = float(kw.get('H', 1.0))
                bg_params[2] = float(kw.get('x_c', 0.0))

        self._U_background = value

        # Push to live params if initialized
        if self._params is not None:
            self._params['U_bg_type']   = tf.constant(bg_type, dtype=tf.int32)
            self._params['U_bg_params'] = tf.constant(bg_params, dtype=DTYPE)

    def _particle_index_range(self, spec):
        """Return (start, end) indices in self._particles for the given spec."""
        offset = 0
        for sp in self._specs:
            if sp is spec:
                return offset, offset + sp.count
            offset += sp.count
        raise ValueError("spec not registered with this System.")

    def _push_param_for_spec(self, spec, param_name, value):
        """
        Patch params[param_name] for all particles belonging to spec.
        Called by ParticleSpec setters when post-init.
        """
        if self._params is None:
            return   # pre-init — nothing to push yet
        import tensorflow as tf
        from src.simulation.tf_sim import DTYPE
        start, end = self._particle_index_range(spec)
        indices = list(range(start, end))
        arr = self._params[param_name].numpy().copy()
        arr[indices] = float(value)
        self._params[param_name] = tf.constant(arr, dtype=DTYPE)

    def set_param(self, name, value, particles='all'):
        """
        Update a material parameter for specific particles.

        name      : str  — parameter name. Supported:
                    High-level:  'Oh', 'kappa', 'nu'
                    Mid-level:   'xi', 'K_area', 'El_t'
        value     : float or array-like — new value(s)
        particles : 'all', list of int indices, or range object

        Works pre- and post-initialization.
        Post-init: patches live TF params tensors immediately.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import DTYPE

        P = sum(sp.count for sp in self._specs)

        # Resolve particle indices
        if particles == 'all':
            idxs = list(range(P))
        else:
            idxs = list(particles)

        # Build index → spec mapping
        spec_of = {}
        offset = 0
        for sp in self._specs:
            for k in range(sp.count):
                spec_of[offset + k] = (sp, k)
            offset += sp.count

        # Group by spec
        from collections import defaultdict
        spec_groups = defaultdict(list)
        for i in idxs:
            sp, local_k = spec_of[i]
            spec_groups[id(sp)].append((i, sp, local_k))

        # Apply per-spec
        for group in spec_groups.values():
            sp = group[0][1]
            if name == 'Oh':
                sp.set_Oh(float(value))
            elif name == 'kappa':
                sp.set_kappa(float(value))
            elif name == 'nu':
                sp.set_nu(float(value))
            elif name == 'xi':
                sp.set_xi(float(value))
            elif name in ('K_area', 'El_t', 'xi_drag_per_p'):
                # Direct mid-level patch
                sp_val = float(value)
                if hasattr(sp, f'_{name}'):
                    setattr(sp, f'_{name}', sp_val)
                if self._params is not None:
                    arr = self._params[name].numpy().copy()
                    for (i, _, _k) in group:
                        arr[i] = sp_val
                    self._params[name] = tf.constant(arr, dtype=DTYPE)
            elif name == 'U_max':
                # Update U_bg_params[0] (U_max for poiseuille_v and parabolic presets).
                # Does not loop over per-spec groups — applies globally.
                if self._params is not None:
                    bg_p = self._params['U_bg_params'].numpy().copy()
                    bg_p[0] = float(value)
                    self._params['U_bg_params'] = tf.constant(bg_p, dtype=DTYPE)
                break  # single global update; no need to iterate over all groups
            else:
                raise ValueError(f"set_param: unknown parameter '{name}'. "
                                 f"Supported: 'Oh', 'kappa', 'nu', 'xi', 'K_area', 'El_t', "
                                 f"'U_max'.")

    # ── R0-scaling recompute ──────────────────────────────────────────────────

    def adjust_params_for_size(self):
        """
        Recompute per-particle physics parameters from each particle's actual R0.

        The actual R0 is read from the current node geometry (mean chord length
        times N / 2π).  All force constants (El_t, EI, k_c, K_area, gamma, C,
        m_node, M_disk, I_disk, alpha_damp, xi_drag) are updated via the R0
        property setter, which calls _recompute_params() on each particle.

        After particle attributes are updated, the live TF params tensors are
        refreshed so the running simulation sees the new values immediately.

        Call this after any geometry rescaling or after manually changing
        particle positions/sizes.  Safe to call multiple times (idempotent).
        """
        for p in self._particles:
            # Actual R0 from geometry: L0 (mean chord) × N / (2π)
            r0_actual = p.L0 * p.N / (2.0 * np.pi)
            p._R0 = r0_actual          # set backing store directly (no infinite loop)
            p._recompute_params()      # update force constants

            # Rescale per-particle damping using the spec's reference values if
            # the particle carries target info; otherwise leave as-is.
            # alpha ∝ R0^{-1}: find spec's alpha at R0_mean and rescale
            if hasattr(p, '_Oh_target') and p._Oh_target is not None:
                # Back out alpha_ref at R0_mean from current alpha_damp * R0
                alpha_r0 = p.alpha_damp * r0_actual  # = alpha_ref * R0_mean (invariant)
                p.alpha_damp = alpha_r0 / r0_actual  # = same (no-op, already correct)
            # xi ∝ R0^{+1}: xi / R0 is constant
            if p.xi_drag > 0.0:
                xi_r0 = p.xi_drag / r0_actual   # = xi_ref / R0_mean (invariant)
                p.xi_drag = xi_r0 * r0_actual   # = same (no-op, already correct)

        if self._params is not None:
            self._refresh_params()

    def _refresh_params(self):
        """
        Push current per-particle attributes into the live TF params tensors.
        Called after adjust_params_for_size() or any bulk particle modification.
        Invalidates the compiled runner so the next run_fast call retraces once.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import DTYPE

        ps = self._particles
        self._params['El_t']             = tf.constant(np.array([p.El_t    for p in ps]), dtype=DTYPE)
        self._params['EI']               = tf.constant(np.array([p.EI      for p in ps]), dtype=DTYPE)
        self._params['m_node']           = tf.constant(np.array([p.m_node  for p in ps]), dtype=DTYPE)
        self._params['M_disk']           = tf.constant(np.array([p.M_disk  for p in ps]), dtype=DTYPE)
        self._params['I_disk']           = tf.constant(np.array([p.I_disk  for p in ps]), dtype=DTYPE)
        self._params['r_c_per_p']        = tf.constant(np.array([p.r_c     for p in ps]), dtype=DTYPE)
        self._params['k_c_per_p']        = tf.constant(np.array([p.k_c     for p in ps]), dtype=DTYPE)
        self._params['gamma_lt']         = tf.constant(np.array([float(getattr(p, 'gamma', 0.0)) for p in ps]), dtype=DTYPE)
        self._params['k_reg_per_p']      = tf.constant(np.array([float(getattr(p, 'k_reg',  0.0)) for p in ps]), dtype=DTYPE)
        self._params['alpha_damp_per_p'] = tf.constant(np.array([p.alpha_damp for p in ps]), dtype=DTYPE)
        self._params['xi_drag_per_p']    = tf.constant(np.array([p.xi_drag    for p in ps]), dtype=DTYPE)

# ── stepping ──────────────────────────────────────────────────────────────

    def _build_prim_data(self, t):
        """Rebuild prim_data from registered objects at time t."""
        from src.simulation.tf_sim import make_prim_data
        prim_list = []
        for obj in self._objects:
            prim_list.extend(obj.to_make_prim_list(t=float(t)))
        return make_prim_data(prim_list)

    def _has_moving_objects(self):
        """True if any registered object has time-varying geometry."""
        return any(obj.is_time_varying() for obj in self._objects)

    def _wrap_state(self, state):
        """Wrap CMs (and nodes) into periodic box."""
        from src.epd.initializer import _wrap
        if self._config['periodic_x'] or self._config['periodic_y']:
            return _wrap(state, self.Lx, self.Ly)
        return state

    def step(self, N=1):
        """
        Advance simulation by N dynamics steps.
        Updates state in place. Returns self.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import step_full_tf, DTYPE, NP_DTYPE

        has_mov = self._has_moving_objects()

        # Extract mutable t/step scalars (step_full_tf doesn't carry them)
        t_val   = float(self._state['t'])
        step_val = int(self._state.get('step', tf.constant(0, dtype=tf.int64)))
        # Physics state (without t/step — step_full_tf works on these keys)
        phys_state = {k: v for k, v in self._state.items()
                      if k not in ('t', 'step')}

        for _ in range(N):
            # Rebuild prim_data if objects are time-varying
            if has_mov:
                self._prim_data = self._build_prim_data(t_val)

            # Candidacy update
            x_cm_np  = phys_state['x_cm'].numpy()
            theta_np = phys_state['theta'].numpy()
            if self._cm_mgr.needs_update(x_cm_np, theta_np):
                self._cm_mgr.update(x_cm_np, theta_np,
                                     x_all=phys_state['x_all'].numpy())

            # TF step
            t_tf = tf.constant(NP_DTYPE(t_val), dtype=DTYPE)
            phys_state, metrics = step_full_tf(
                phys_state,
                self._cm_mgr.CapCandidates,
                self._dt_tf,
                self._alpha_tf,
                self._g_tf,
                self._params,
                t=t_tf,
                prim_data=self._prim_data,
            )

            # Stability watchdog: per-contact closing rate ρ = |Δ(rel_disp)|/r_c
            ratio = float(metrics['max_closing_ratio'])
            if ratio > self._max_disp_ratio_run:
                self._max_disp_ratio_run = ratio
            if ratio > self._disp_critical:
                rec = max(1.0, ratio / self._disp_advisory)
                print(f"CRITICAL: step {step_val}  ρ_contact = {ratio:.3f}  "
                      f">  critical {self._disp_critical:.2f}.  "
                      f"Likely unstable; reduce dt_factor by ≥ {rec:.1f}× "
                      f"(currently {self._dt_factor:.3f}).")
            elif ratio > self._disp_strong:
                rec = max(1.0, ratio / self._disp_advisory)
                print(f"WARNING: step {step_val}  ρ_contact = {ratio:.3f}  "
                      f">  strong {self._disp_strong:.2f}.  "
                      f"Reduce dt_factor by ~{rec:.1f}× "
                      f"(currently {self._dt_factor:.3f}).")
            elif ratio > self._disp_advisory and not self._disp_advisory_fired:
                rec = max(1.0, ratio / self._disp_advisory)
                print(f"NOTE: step {step_val}  ρ_contact = {ratio:.3f}  "
                      f">  advisory {self._disp_advisory:.2f}.  "
                      f"Approaching the closing-rate stability limit; "
                      f"consider dt_factor ≤ {self._dt_factor / rec:.3f}.  "
                      f"(One-shot — won't fire again until clear_recording().)")
                self._disp_advisory_fired = True

            # Periodic wrap
            phys_state = self._wrap_state(phys_state)

            # Advance time and step counter
            t_val    += self._dt
            step_val += 1

        # Reattach t/step to physics state
        self._state = {
            **phys_state,
            't':    tf.constant(NP_DTYPE(t_val), dtype=DTYPE),
            'step': tf.constant(step_val, dtype=tf.int64),
        }
        return self

    # ── tf-fast chunked run ───────────────────────────────────────────────────

    def run_fast(self, n_steps, cand_check_interval=10, diagnostics=False,
                 use_tf_function=None):
        """
        Run n_steps using run_simulation_tf (tf.while_loop, static prim_data).

        Wall oscillation is encoded in prim_data at initialize time — prim_data
        is never rebuilt. step_offset is inferred from self.step_count so that
        global simulation time t = (step_offset + i) * dt is correct across
        chunk boundaries.

        Parameters
        ----------
        n_steps             : int   — steps to advance
        cand_check_interval : int   — Python candidacy callback every N steps (default 10)
        diagnostics         : bool  — if True, return diag dict; else return self
        use_tf_function     : bool or None — when not None, forwarded to
                              run_simulation_tf. When None (default), falls back
                              to self.use_tf_function (set on the System
                              constructor or via the property). When True, the
                              inner tf.while_loop runs as a single graph
                              (fewer per-iteration CPU↔GPU round-trips, larger
                              gain on virtualised PCIe like Colab).

        Returns
        -------
        self      — if diagnostics=False (fluent API)
        diag dict — if diagnostics=True (for test/validation scripts)
        """
        import tensorflow as tf
        from src.simulation.tf_sim import run_simulation_tf, DTYPE, NP_DTYPE

        if use_tf_function is None:
            use_tf_function = self._use_tf_function

        step_offset = self.step_count
        phys_state  = {k: v for k, v in self._state.items() if k not in ('t', 'step')}

        R0_arr   = np.array([p.R0 for p in self._particles])
        skin_abs = self._skin * float(np.mean(R0_arr))
        R0_max   = float(np.max(R0_arr))

        result = run_simulation_tf(
            phys_state,
            self._dt,
            float(self._alpha_damp),
            self._g_val,
            self._params,
            n_steps,
            self._cm_mgr,
            skin=skin_abs,
            prim_data=self._prim_data,
            R0_max=R0_max,
            cand_check_interval=cand_check_interval,
            diagnostics=diagnostics,
            step_offset=step_offset,
            use_tf_function=use_tf_function,
        )

        if diagnostics:
            new_phys_state, diag = result
        else:
            new_phys_state = result
            diag = None

        t_val    = (step_offset + n_steps) * self._dt
        step_val = step_offset + n_steps
        self._state = {
            **new_phys_state,
            't':    tf.constant(NP_DTYPE(t_val), dtype=DTYPE),
            'step': tf.constant(step_val,        dtype=tf.int64),
        }

        return diag if diagnostics else self

    def snapshot(self):
        """
        Return current state as a numpy dict with CMs/nodes wrapped into [0,Lx)×[0,Ly).
        Wrapping is done in numpy here so the live TF state is never perturbed.
        """
        x_all = self._state['x_all'].numpy()
        x_cm  = self._state['x_cm'].numpy()
        if self._config['periodic_x']:
            dx = x_cm[:, 0] % self.Lx - x_cm[:, 0]
            mx = np.abs(dx) > 1e-9
            if mx.any():
                x_cm[mx, 0]     += dx[mx]
                x_all[mx, :, 0] += dx[mx, None]
        if self._config['periodic_y']:
            dy = x_cm[:, 1] % self.Ly - x_cm[:, 1]
            my = np.abs(dy) > 1e-9
            if my.any():
                x_cm[my, 1]     += dy[my]
                x_all[my, :, 1] += dy[my, None]
        return {
            'x_all': x_all,
            'x_cm':  x_cm,
            'theta': self._state['theta'].numpy(),
            't':     float(self._state['t']),
            'step':  int(self._state['step']),
        }

    def eval_forces(self):
        """
        Compute per-node force breakdown at the current simulation state.

        Returns a dict with keys (P, N, 2) numpy arrays:
          'f_elastic'  — internal elastic + area + line-tension forces
          'f_contact'  — inter-capsule and wall contact forces
          'f_drag'     — Stokes drag (zero if Oh not set)
          'f_reg'      — tangential regularisation pseudo-force (emulsion only)
          'f_total'    — sum of all four

        Typical use inside a callback::

            def my_cb(sys):
                fb = sys.eval_forces()
                return {'F_contact_mean': float(np.linalg.norm(fb['f_contact'], axis=-1).mean())}

        Or post-hoc from recorded frames (re-run one step from a saved state).
        """
        import tensorflow as tf
        from src.simulation.tf_sim import eval_forces_tf
        caps = tf.constant(self._cm_mgr.CapCandidates, dtype=tf.int32)
        zero = tf.constant(0.0, dtype=tf.float64)
        return eval_forces_tf(
            state         = self._state,
            CapCandidates = caps,
            alpha_damp    = self._params.get('_alpha_tf', zero),
            g             = self._params.get('_g_tf',     zero),
            params        = self._params,
            t             = self._state['t'],
            prim_data     = self._prim_data,
        )

    # ── batch run ─────────────────────────────────────────────────────────────

    def run(self, N, sample_every=None, callback=None, record_initial=True,
            cand_check_interval=10, verbose=True, use_tf_function=None):
        """
        Run N steps in chunks of sample_every, recording to self.frames /
        self.callback_data / self.diag.

        sample_every is the single cadence knob: chunk size = sampling =
        callback frequency.  For per-step debugging set sample_every=1.

        Appends to self.frames, self.callback_data, self.diag on every call.
        Call clear_recording() to wipe without touching simulation state.
        Returns self (fluent API).

        Parameters
        ----------
        N                   : int  — total steps to run
        sample_every        : int  — chunk size AND output cadence (default: N)
        callback            : callable(sys) → dict or None
        record_initial      : bool — record snapshot/callback before first step
        cand_check_interval : int  — C++ candidacy poll interval (default 10)
        verbose             : bool — print progress every 20% of chunks (default True)
        use_tf_function     : bool or None — when not None, forwarded to
                              run_fast → run_simulation_tf. When None (default),
                              falls back to self.use_tf_function (set on the
                              System constructor or via the property). When True,
                              each chunk's tf.while_loop runs as a single graph
                              (fewer per-iteration CPU↔GPU round-trips, larger
                              gain on virtualised PCIe like Colab). First chunk
                              pays one ~5–10s trace cost.
        """
        import time as _time
        if use_tf_function is None:
            use_tf_function = self._use_tf_function
        if sample_every is None:
            sample_every = N

        n_full      = N // sample_every
        remain      = N  % sample_every
        n_chunks    = n_full + (1 if remain > 0 else 0)
        _t0         = _time.time()

        def _record(step_i):
            snap = self.snapshot()
            self.frames.append(snap)
            cbd = callback(self) if callback is not None else {}
            self.callback_data.append(cbd if cbd is not None else {})

        if record_initial:
            _record(0)

        for ci in range(n_full):
            d = self.run_fast(sample_every,
                              cand_check_interval=cand_check_interval,
                              diagnostics=True,
                              use_tf_function=use_tf_function)
            self.diag.append(d)
            self._check_disp_watchdog(d.get('max_closing_ratio_chunk', 0.0), ci + 1)
            _record(ci + 1)
            if verbose and n_chunks >= 5 and (ci + 1) % max(1, n_chunks // 5) == 0:
                elapsed = _time.time() - _t0
                print(f"    chunk {ci+1}/{n_chunks}  t={self.t:.3f}  ({elapsed:.0f}s)")

        if remain > 0:
            d = self.run_fast(remain,
                              cand_check_interval=cand_check_interval,
                              diagnostics=True,
                              use_tf_function=use_tf_function)
            self.diag.append(d)
            self._check_disp_watchdog(d.get('max_closing_ratio_chunk', 0.0), n_full + 1)
            _record(n_full + 1)
            if verbose and n_chunks >= 5:
                elapsed = _time.time() - _t0
                print(f"    chunk {n_chunks}/{n_chunks}  t={self.t:.3f}  ({elapsed:.0f}s)")

        return self

    def _check_disp_watchdog(self, ratio, chunk_idx):
        """Stability watchdog — fires when max_closing_ratio_chunk
        crosses configured thresholds. Called once per chunk by run().
        ratio = max per-contact closing rate (|Δ(rel_disp)| / r_c_pair)
        observed during the chunk, including particle-particle and
        particle-wall pairs.
        """
        if ratio <= 0.0:
            return
        if ratio > self._max_disp_ratio_run:
            self._max_disp_ratio_run = ratio
        if ratio > self._disp_critical:
            rec = max(1.0, ratio / self._disp_advisory)
            print(f"CRITICAL: chunk {chunk_idx}  ρ_contact = {ratio:.3f}  "
                  f">  critical {self._disp_critical:.2f}.  "
                  f"Likely unstable; reduce dt_factor by ≥ {rec:.1f}x "
                  f"(currently {self._dt_factor:.3f}).")
        elif ratio > self._disp_strong:
            rec = max(1.0, ratio / self._disp_advisory)
            print(f"WARNING (watchdog): chunk {chunk_idx}  ρ_contact = {ratio:.3f}  "
                  f">  strong {self._disp_strong:.2f}.  "
                  f"Reduce dt_factor by ~{rec:.1f}x "
                  f"(currently {self._dt_factor:.3f}).")
        elif ratio > self._disp_advisory and not self._disp_advisory_fired:
            rec = max(1.0, ratio / self._disp_advisory)
            print(f"NOTE (watchdog): chunk {chunk_idx}  ρ_contact = {ratio:.3f}  "
                  f">  advisory {self._disp_advisory:.2f}.  "
                  f"Approaching closing-rate stability limit; "
                  f"consider dt_factor <= {self._dt_factor / rec:.3f}.  "
                  f"(One-shot — won't fire again until clear_recording().)")
            self._disp_advisory_fired = True

    # ── lightweight state save / restore ──────────────────────────────────────

    def save_state(self, path):
        """
        Save post-swell simulation state to a .npz file.

        Stores all TF state tensors (positions, velocities, internal DOFs),
        box dimensions, dt, and alpha_damp.  Used to skip expensive swells on
        subsequent runs.  Requires initialize() to have been called first.

        Returns pathlib.Path to the saved file.
        """
        import pathlib
        out = pathlib.Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        arrays = {
            'Lx':         np.float64(self.Lx),
            'Ly':         np.float64(self.Ly),
            'dt':         np.float64(self._dt),
            'alpha_damp': np.float64(self._alpha_damp),
        }
        for k, v in self._state.items():
            if hasattr(v, 'numpy'):
                arrays[k] = v.numpy()
        # Persist registered-object geometry exactly so restore_state doesn't
        # have to reconstruct via affine rescale (which can drift relative to
        # particle positions over many compression cycles).
        for oi, obj in enumerate(self._objects):
            cls = obj.__class__.__name__
            arrays[f'obj{oi}_class'] = np.array(cls)
            if hasattr(obj, '_origin'):
                arrays[f'obj{oi}_origin'] = np.asarray(obj._origin, dtype=np.float64)
            if hasattr(obj, '_inner_r'):
                arrays[f'obj{oi}_inner_r'] = np.float64(obj._inner_r)
            if hasattr(obj, '_outer_r'):
                arrays[f'obj{oi}_outer_r'] = np.float64(obj._outer_r)
        np.savez(str(out), **arrays)
        return out

    def restore_state(self, path):
        """
        Overwrite current state from a .npz written by save_state().

        System must already be initialize()-d (positions will be overwritten).
        Rescales registered objects to match checkpoint box dimensions, then
        rebuilds prim_data and CandidacyManager.

        Returns self.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import (make_prim_data, DTYPE, NP_DTYPE,
                                           set_periodic_box)
        from src.simulation.candidacy_manager import CandidacyManager

        ckpt     = np.load(path)
        Lx_ckpt  = float(ckpt['Lx'])
        Ly_ckpt  = float(ckpt['Ly'])

        # Rescale objects from current box to checkpoint box (fallback path)
        f = Lx_ckpt / self.Lx
        for obj in self._objects:
            if hasattr(obj, 'rescale'):
                obj.rescale(f)
        self.Lx = Lx_ckpt
        self.Ly = Ly_ckpt

        # Override with EXACT object geometry if persisted in the checkpoint.
        # This avoids accumulated drift between particle positions and
        # affine-reconstructed object positions over many compression cycles.
        for oi, obj in enumerate(self._objects):
            ok_orig = f'obj{oi}_origin'  in ckpt.files
            ok_inr  = f'obj{oi}_inner_r' in ckpt.files
            ok_outr = f'obj{oi}_outer_r' in ckpt.files
            if not (ok_orig or ok_inr or ok_outr):
                continue
            if ok_orig and hasattr(obj, '_origin'):
                obj._origin = np.array(ckpt[f'obj{oi}_origin'], dtype=np.float64)
            if ok_inr and hasattr(obj, '_inner_r'):
                obj._inner_r = float(ckpt[f'obj{oi}_inner_r'])
            if ok_outr and hasattr(obj, '_outer_r'):
                obj._outer_r = float(ckpt[f'obj{oi}_outer_r'])
            # If this is a CompositeObject, rebuild children to reflect new radii
            if hasattr(obj, '_children') and ok_inr and ok_outr:
                from src.epd.objects import ArcWall
                obj._children = []
                outer = ArcWall((0, 0), obj._outer_r, convex=False)
                outer.set_exclusion('exterior')
                obj.add_primitive(outer)
                inner = ArcWall((0, 0), obj._inner_r, convex=True)
                inner.set_exclusion('interior')
                obj.add_primitive(inner)

        # Restore integration scalars
        self._dt         = float(ckpt['dt'])
        self._alpha_damp = float(ckpt['alpha_damp'])
        self._dt_tf      = tf.constant(NP_DTYPE(self._dt))
        self._alpha_tf   = tf.constant(NP_DTYPE(self._alpha_damp))
        # Keep dt_factor consistent with restored dt (dt_max is unchanged)
        if self._dt_max is not None and self._dt_max > 0:
            self._dt_factor = self._dt / self._dt_max

        # Update periodic box in params
        if self._config['periodic_x'] or self._config['periodic_y']:
            set_periodic_box(self._params, self.Lx, self.Ly,
                              periodic_x=self._config['periodic_x'],
                              periodic_y=self._config['periodic_y'])

        # Overwrite state tensors from checkpoint
        new_state = {}
        for k, v in self._state.items():
            if k == 't' and 't' in ckpt:
                new_state[k] = tf.constant(NP_DTYPE(float(ckpt['t'])), dtype=DTYPE)
            elif k == 'step' and 'step' in ckpt:
                new_state[k] = tf.constant(int(ckpt['step']), dtype=tf.int64)
            elif k in ckpt:
                new_state[k] = tf.constant(ckpt[k], dtype=DTYPE)
            else:
                new_state[k] = v
        self._state = new_state

        # Rebuild prim_data at checkpoint box scale
        prim_list = []
        for obj in self._objects:
            prim_list.extend(obj.to_make_prim_list(t=0.0))
        self._prim_data = make_prim_data(prim_list)

        # Rebuild CandidacyManager at new scale
        R0_arr = np.array([p.R0 for p in self._particles])
        N_val  = self._particles[0].N
        self._cm_mgr = self._build_cm(len(self._particles), N_val, R0_arr,
                                        self._particles)
        self._cm_mgr.update(self._state['x_cm'].numpy(),
                            self._state['theta'].numpy(),
                            x_all=self._state['x_all'].numpy())
        return self

    def set_driven_particles(self, indices, traj_rows, frozen=True):
        """
        Assign prescribed motion to particle indices after initialize().

        When frozen=True the particle's current deformed shape (post-swell u)
        is captured and stored as params['u_frozen'], so the shape retained
        throughout the run is the deformed post-swell shape, not the circular
        reference shape.

        indices   : list of int
        traj_rows : list of (18,) arrays from make_traj(), one per index
        frozen    : bool or list of bool — freeze elastic shape (default True)

        Returns self.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import set_driven, DTYPE
        if isinstance(frozen, bool):
            frozen_flags = [frozen] * len(indices)
        else:
            frozen_flags = list(frozen)

        # Capture current deformed u for each particle that will be frozen
        if any(frozen_flags):
            P         = len(self._particles)
            N         = self._particles[0].N
            u_current = self._state['u'].numpy()   # (P, N, 2)
            if 'u_frozen' in self._params:
                u_frozen = self._params['u_frozen'].numpy().copy()
            else:
                u_frozen = np.zeros((P, N, 2), dtype=u_current.dtype)
            for idx, frz in zip(indices, frozen_flags):
                if frz:
                    u_frozen[idx] = u_current[idx]
            self._params['u_frozen'] = tf.constant(u_frozen, dtype=DTYPE)

        set_driven(self._params, list(indices), list(traj_rows),
                   frozen=frozen_flags)
        return self

    # ── checkpoint ────────────────────────────────────────────────────────────

    def save(self, path):
        """
        Save full system state to directory path.
        Returns pathlib.Path to checkpoint directory.
        """
        from src.epd.checkpoint import save_checkpoint
        return save_checkpoint(path, self)

    @classmethod
    def from_file(cls, path):
        """
        Load a System from a checkpoint directory and return it.
        Equivalent to System(...).load(path).
        """
        from src.epd.checkpoint import load_checkpoint
        return load_checkpoint(path, system=None)

    def load(self, path):
        """
        Reload state from checkpoint directory into this System.
        Returns self.
        """
        from src.epd.checkpoint import load_checkpoint
        loaded = load_checkpoint(path, system=None)
        # Copy all internals from loaded into self
        self.__dict__.update(loaded.__dict__)
        return self

    # ── utilities ─────────────────────────────────────────────────────────────

    def freeze(self, zero_velocity=True):
        """
        Zero velocities (v_cm, omega, u_dot) to prepare for production run.
        Returns self.
        """
        import tensorflow as tf
        from src.simulation.tf_sim import DTYPE
        st = self._state
        updates = {}
        if zero_velocity:
            updates['v_cm']  = tf.zeros_like(st['v_cm'])
            updates['omega'] = tf.zeros_like(st['omega'])
            updates['u_dot'] = tf.zeros_like(st['u_dot'])
        self._state = {**st, **updates}
        return self

    # ── rendering ─────────────────────────────────────────────────────────────

    @staticmethod
    def _outer_contour(xy, x_cm, r_c):
        """
        Outer capsule contour: each node offset radially outward from CM by r_c.
        Returns (N, 2) array of outer-contour points.
        """
        dr   = xy - x_cm[None, :]
        lens = np.linalg.norm(dr, axis=1, keepdims=True)
        lens = np.maximum(lens, 1e-12)
        return xy + r_c * (dr / lens)

    def render(self, ax=None, output_path=None, t_label=True, monitor_line=None):
        """
        Draw current particle configuration using the outer capsule contour
        (each node shifted outward by r_c from the particle CM).

        Parameters
        ----------
        monitor_line : str or None — health summary to print below title
        """
        created = ax is None
        if created:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.get_figure()

        ax.set_aspect('equal')
        ax.set_facecolor('#f0f0f0')

        # Axis limits: encompass the box with a small margin
        margin = max(self.Lx, self.Ly) * 0.03
        ax.set_xlim(-margin, self.Lx + margin)
        ax.set_ylim(-margin, self.Ly + margin)

        x_all   = self._state['x_all'].numpy()
        x_cm_np = self._state['x_cm'].numpy()
        P       = x_all.shape[0]
        r_c_arr = self._params['r_c_per_p'].numpy()

        # Distinguish elastic vs emulsion by particle type (drives edge styling)
        from src.simulation.emulsion_particle import EmulsionParticle
        is_emulsion = [isinstance(p, EmulsionParticle) for p in self._particles]

        # Fallback colormaps used only if _particle_colors isn't populated
        # (e.g. legacy callers that bypass initialize). Honoring _particle_colors
        # keeps render() consistent with make_movie() and respects user-set
        # palettes from set_color_palette().
        use_user_colors = (len(self._particle_colors) == P)
        colors_e  = plt.cm.tab20(np.linspace(0, 1, max(P, 1)))
        colors_em = plt.cm.cool(np.linspace(0.2, 0.9, max(P, 1)))
        ei = 0;  emi = 0

        # Periodic image offsets: draw copies that straddle box edges
        px = self._config.get('periodic_x', False)
        py = self._config.get('periodic_y', False)
        dx_offsets = [-self.Lx, 0.0, self.Lx] if px else [0.0]
        dy_offsets = [-self.Ly, 0.0, self.Ly] if py else [0.0]
        img_offsets = [(dx, dy) for dx in dx_offsets for dy in dy_offsets]

        for i in range(P):
            outer = self._outer_contour(x_all[i], x_cm_np[i], r_c_arr[i])
            if is_emulsion[i]:
                c = self._particle_colors[i] if use_user_colors else colors_em[emi % len(colors_em)]
                emi += 1
                fc, ec, lw, alpha = c, 'steelblue', 0.9, 0.70
            else:
                c = self._particle_colors[i] if use_user_colors else colors_e[ei % len(colors_e)]
                ei  += 1
                fc, ec, lw, alpha = c, 'k', 0.7, 0.72

            for odx, ody in img_offsets:
                shifted = outer + np.array([odx, ody])
                # Only draw copies that overlap the visible canvas
                if (shifted[:, 0].max() < -margin or shifted[:, 0].min() > self.Lx + margin or
                        shifted[:, 1].max() < -margin or shifted[:, 1].min() > self.Ly + margin):
                    continue
                ax.fill(shifted[:, 0], shifted[:, 1], color=fc, alpha=alpha)
                xs = np.append(shifted[:, 0], shifted[0, 0])
                ys = np.append(shifted[:, 1], shifted[0, 1])
                ax.plot(xs, ys, color=ec, lw=lw)

        # Draw simulation objects (walls, Couette cell, etc.)
        self._render_objects(ax)

        # Box outline (dashed)
        ax.plot([0, self.Lx, self.Lx, 0, 0],
                [0, 0, self.Ly, self.Ly, 0],
                'b--', lw=0.8, alpha=0.4)

        if t_label:
            t_val = float(self._state['t'])
            sc    = int(self._state['step'])
            title = f't={t_val:.3f}  step={sc}  φ={self.phi_outer:.3f}'
            if monitor_line:
                title += f'\n{monitor_line}'
            ax.set_title(title, fontsize=8)

        ax.tick_params(labelsize=7)

        if output_path is not None:
            fig.set_size_inches(6, 6)
            fig.savefig(output_path, dpi=110)
            plt.close(fig)

        return fig

    def _render_objects(self, ax, t=None):
        """Draw registered objects on ax at time t (default: current state t)."""
        import matplotlib.patches as mpatches
        from src.simulation.contact_primitives import LineSegment, Arc
        if t is None:
            t = float(self._state.get('t', 0.0)) if self._state else 0.0

        for obj in self._objects:
            rs    = getattr(obj, '_render_style', {})
            color = rs.get('color', 'k')
            lw    = rs.get('linewidth', 1.5)
            alpha = rs.get('alpha', 1.0)
            prims = obj.resolved(t=float(t))
            for pd in prims:
                prim = pd['prim']
                if isinstance(prim, LineSegment):
                    ax.plot([prim.p0[0], prim.p1[0]],
                            [prim.p0[1], prim.p1[1]],
                            '-', color=color, lw=lw, alpha=alpha)
                elif isinstance(prim, Arc):
                    if prim.angle_range is not None:
                        import math as _math
                        t0_deg = _math.degrees(prim.angle_range[0])
                        t1_deg = _math.degrees(prim.angle_range[1])
                        arc_patch = mpatches.Arc(
                            prim.center, 2 * prim.radius, 2 * prim.radius,
                            angle=0, theta1=t0_deg, theta2=t1_deg,
                            color=color, lw=lw, alpha=alpha)
                        ax.add_patch(arc_patch)
                    else:
                        circ = mpatches.Circle(
                            prim.center, prim.radius,
                            fill=False, ec=color, lw=lw, alpha=alpha,
                            linestyle='-' if not prim.convex else '--')
                        ax.add_patch(circ)

    @staticmethod
    def _capsule_outline_polygon(x, r_c, n_arc=6):
        """
        Rounded capsule outline used by make_squeeze_gif.
        Matches the rendering in the TF reference test scripts exactly.
        x     : (N, 2) node positions
        r_c   : float capsule contact radius
        Returns (M, 2) closed polygon ready for Matplotlib Polygon patch.
        """
        N  = len(x)
        xn = np.roll(x, -1, axis=0)
        e  = xn - x
        L  = np.linalg.norm(e, axis=1, keepdims=True)
        t  = e / np.where(L > 1e-15, L, 1.0)
        n_edge = np.column_stack([t[:, 1], -t[:, 0]])
        pts = []
        for i in range(N):
            i_prev = (i - 1) % N
            a1 = np.arctan2(n_edge[i_prev, 1], n_edge[i_prev, 0])
            a2 = np.arctan2(n_edge[i,      1], n_edge[i,      0])
            da = (a2 - a1) % (2 * np.pi)
            if da < 1e-10:
                da = 2 * np.pi
            if da <= np.pi:
                n_pts  = max(2, int(n_arc * da / np.pi) + 1)
                thetas = np.linspace(a1, a1 + da, n_pts)
                arc    = x[i] + r_c * np.column_stack([np.cos(thetas),
                                                        np.sin(thetas)])
            else:
                arc = x[i] + r_c * np.vstack([n_edge[i_prev], n_edge[i]])
            pts.append(arc)
        outline = np.vstack(pts)
        return np.vstack([outline, outline[:1]])

    def make_movie(self, output_path, fps=10, n_arc=6,
                   xlim=None, ylim=None, title=None,
                   bitrate=None, dpi=120,
                   frames=None, callback_data=None):
        """
        Render an animated movie from self.frames and self.callback_data.

        Output format is selected by file extension:
            .gif           → matplotlib pillow writer (lossless, big files)
            .mp4 / .webm   → matplotlib ffmpeg writer (compressed, much smaller)
                              ffmpeg binary is sourced from the imageio_ffmpeg
                              package if not already on PATH.

        Particle shapes come from self.frames[i]['x_all'] (snapshots).
        Particle fill colors from self._particle_colors (auto-assigned or set via
        p.color attribute).
        Wall / object styles from obj._render_style (set via obj.set_render()).
        Overlay text from self.callback_data[i].get('text', '').

        A time-series panel is added automatically when callback_data contains
        'wall_strain' and/or 'eps1'/'eps2' keys.

        Parameters
        ----------
        bitrate : int or None  – ffmpeg target bitrate in kbps (mp4/webm only).
                                 None → matplotlib default (~1024 kbps for h264).
        dpi     : int          – figure dpi for rasterisation (default 120).
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import Polygon as MplPolygon
        import matplotlib.animation as animation
        import pathlib
        from src.simulation.contact_primitives import LineSegment, Arc

        snap_frames = frames if frames is not None else self.frames
        cb          = callback_data if callback_data is not None else self.callback_data
        if not snap_frames:
            raise RuntimeError("No frames provided; call run() first or pass frames=.")

        out_path = pathlib.Path(output_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        n_frames    = len(snap_frames)
        P           = len(self._particles)
        r_c_arr     = self._params['r_c_per_p'].numpy()

        # Auto xlim / ylim: scan particles + object extents
        if xlim is None or ylim is None:
            all_x, all_y = [], []
            for fr in snap_frames:
                xy = fr['x_all'].reshape(-1, 2)
                all_x.extend(xy[:, 0]); all_y.extend(xy[:, 1])
            sample_ts = {snap_frames[0]['t'], snap_frames[-1]['t']}
            if len(snap_frames) > 2:
                sample_ts.add(snap_frames[len(snap_frames)//2]['t'])
            for t_s in sample_ts:
                for obj in self._objects:
                    for pd in obj.resolved(t=t_s):
                        prim = pd['prim']
                        if isinstance(prim, LineSegment):
                            all_x += [prim.p0[0], prim.p1[0]]
                            all_y += [prim.p0[1], prim.p1[1]]
            r_c_max = float(np.max(r_c_arr)) if len(r_c_arr) else 0.2
            xspan   = max(float(np.max(all_x)) - float(np.min(all_x)), 1e-3)
            yspan   = max(float(np.max(all_y)) - float(np.min(all_y)), 1e-3)
            margin  = r_c_max + 0.12 * max(xspan, yspan)
            if xlim is None:
                xlim = (float(np.min(all_x)) - margin, float(np.max(all_x)) + margin)
            if ylim is None:
                ylim = (float(np.min(all_y)) - margin, float(np.max(all_y)) + margin)

        # margin used for ghost-image visibility culling
        r_c_max = float(np.max(r_c_arr)) if len(r_c_arr) else 0.2
        margin  = r_c_max + 0.05 * max(xlim[1] - xlim[0], ylim[1] - ylim[0])

        # Detect time-series panel from callback keys
        has_ts = (cb and isinstance(cb[0], dict) and
                  any(k in cb[0] for k in ('wall_strain', 'eps1')))

        if has_ts:
            fig, (ax_sim, ax_ts) = plt.subplots(1, 2, figsize=(11, 5.5))
        else:
            fig, ax_sim = plt.subplots(figsize=(5.5, 5.5))

        ax_sim.set_aspect('equal')
        ax_sim.set_xlim(*xlim)
        ax_sim.set_ylim(*ylim)
        ax_sim.set_xlabel('x / R₀')
        ax_sim.set_ylabel('y / R₀')
        if title:
            ax_sim.set_title(title)

        # Particle patches — one primary + periodic ghost copies per particle.
        # img_offsets: list of (dx, dy) shifts; primary is always (0, 0).
        px = self._config.get('periodic_x', False)
        py = self._config.get('periodic_y', False)
        dx_list = [-self.Lx, 0.0, self.Lx] if px else [0.0]
        dy_list = [-self.Ly, 0.0, self.Ly] if py else [0.0]
        img_offsets = [(dx, dy) for dx in dx_list for dy in dy_list]

        # p_patches[pi][k] = MplPolygon for particle pi at offset img_offsets[k]
        p_patches = []
        for pi in range(P):
            c   = self._particle_colors[pi] if pi < len(self._particle_colors) else 'cornflowerblue'
            row = []
            for k, (dx, dy) in enumerate(img_offsets):
                pat = MplPolygon(np.zeros((4, 2)), closed=True,
                                 fc=c, ec='k', lw=0.8, alpha=0.65, zorder=3,
                                 visible=(dx == 0.0 and dy == 0.0))
                ax_sim.add_patch(pat)
                row.append(pat)
            p_patches.append(row)

        # Object artists (one Line2D or Circle per primitive; created at t0)
        obj_artists = []   # list of (obj, artist, kind, prim_idx_in_obj)
        t0 = snap_frames[0]['t']
        for obj in self._objects:
            rs    = getattr(obj, '_render_style', {})
            color = rs.get('color', '#555555')
            lw    = rs.get('linewidth', 2.0)
            alpha = rs.get('alpha', 0.85)
            for ki, pd in enumerate(obj.resolved(t=t0)):
                prim = pd['prim']
                if isinstance(prim, LineSegment):
                    line, = ax_sim.plot(
                        [prim.p0[0], prim.p1[0]], [prim.p0[1], prim.p1[1]],
                        '-', color=color, lw=lw, alpha=alpha, zorder=1)
                    obj_artists.append((obj, line, 'seg', ki))
                elif isinstance(prim, Arc):
                    if prim.angle_range is not None:
                        t0_deg = np.degrees(prim.angle_range[0])
                        t1_deg = np.degrees(prim.angle_range[1])
                    else:
                        t0_deg, t1_deg = 0.0, 360.0
                    arc_artist = mpatches.Arc(
                        prim.center, 2 * prim.radius, 2 * prim.radius,
                        angle=0, theta1=t0_deg, theta2=t1_deg,
                        ec=color, lw=lw, alpha=alpha, zorder=1,
                        linestyle='-' if not prim.convex else '--')
                    ax_sim.add_patch(arc_artist)
                    obj_artists.append((obj, arc_artist, 'arc', ki))

        txt = ax_sim.text(0.02, 0.98, '', transform=ax_sim.transAxes,
                          fontsize=8, va='top', family='monospace')

        # Time-series panel setup
        ts_lines = []
        vline    = None
        t_arr    = None
        if has_ts:
            t_arr       = np.array([fr['t'] for fr in snap_frames])
            strain_arr  = (np.array([c.get('wall_strain', 0.0) for c in cb])
                           if 'wall_strain' in cb[0] else None)
            eps_arr     = (np.array([0.5*(c.get('eps1', 0.0) + c.get('eps2', 0.0))
                                     for c in cb])
                           if 'eps1' in cb[0] else None)

            ax_ts.set_xlabel('t')
            ax_ts.set_ylabel('Value')
            ax_ts.set_title('Time series')
            ax_ts.set_xlim(0, float(t_arr[-1]))
            y_max = 0.02
            if strain_arr is not None: y_max = max(y_max, float(strain_arr.max()) * 1.3)
            if eps_arr    is not None: y_max = max(y_max, float(eps_arr.max())    * 1.3)
            ax_ts.set_ylim(-0.005, y_max)
            ax_ts.axhline(0, color='k', lw=0.5)
            vline = ax_ts.axvline(0, color='gray', lw=0.8, ls='--')
            if strain_arr is not None:
                l, = ax_ts.plot([], [], 'b-', lw=1.5, label='wall strain')
                ts_lines.append((l, strain_arr))
            if eps_arr is not None:
                l, = ax_ts.plot([], [], 'r-', lw=1.5, label='ε_p (avg)')
                ts_lines.append((l, eps_arr))
            ax_ts.legend(fontsize=7)

        plt.tight_layout()

        # Flat list of all particle patches (for blit return value)
        all_p_patches = [pat for row in p_patches for pat in row]

        # Per-object resolved cache to avoid rebuilding each resolved() per artist
        def _get_obj_prim(obj, t_now, idx):
            prims = obj.resolved(t=t_now)
            return prims[idx] if idx < len(prims) else None

        def _update(i):
            fr  = snap_frames[i]
            cbd = cb[i] if i < len(cb) else {}
            t_now = fr['t']

            # Particles — primary + periodic ghost copies
            for pi in range(P):
                outline = self._capsule_outline_polygon(
                    fr['x_all'][pi], float(r_c_arr[pi]), n_arc)
                for k, (odx, ody) in enumerate(img_offsets):
                    shifted = outline + np.array([odx, ody])
                    pat = p_patches[pi][k]
                    # Visibility: only draw copies that overlap the canvas
                    in_canvas = (shifted[:, 0].max() > xlim[0] - margin and
                                 shifted[:, 0].min() < xlim[1] + margin and
                                 shifted[:, 1].max() > ylim[0] - margin and
                                 shifted[:, 1].min() < ylim[1] + margin)
                    pat.set_xy(shifted)
                    pat.set_visible(in_canvas)

            # Objects (moving walls etc.)
            for (obj, artist, kind, ki) in obj_artists:
                pd = _get_obj_prim(obj, t_now, ki)
                if pd is None:
                    continue
                prim = pd['prim']
                if kind == 'seg' and isinstance(prim, LineSegment):
                    artist.set_data([prim.p0[0], prim.p1[0]],
                                    [prim.p0[1], prim.p1[1]])
                elif kind == 'arc' and isinstance(prim, Arc):
                    artist.set_center(tuple(prim.center))

            # Overlay text
            txt.set_text(cbd.get('text', ''))

            artists = all_p_patches + [a for (_, a, _, _) in obj_artists] + [txt]

            # Time series
            if has_ts:
                vline.set_xdata([t_now])
                for (l, arr) in ts_lines:
                    l.set_data(t_arr[:i+1], arr[:i+1])
                artists = artists + [vline] + [l for (l, _) in ts_lines]

            return artists

        ani = animation.FuncAnimation(fig, _update, frames=n_frames,
                                      interval=100, blit=True)

        # Pick writer from output extension; configure ffmpeg path on demand.
        ext = out_path.suffix.lower()
        if ext == '.gif':
            ani.save(str(out_path), writer='pillow', fps=fps, dpi=dpi)
        elif ext in ('.mp4', '.webm', '.mov', '.mkv'):
            import shutil
            ffmpeg_exe = shutil.which('ffmpeg')
            if ffmpeg_exe is None:
                try:
                    import imageio_ffmpeg
                    ffmpeg_exe = imageio_ffmpeg.get_ffmpeg_exe()
                except ImportError:
                    raise RuntimeError(
                        "Need ffmpeg for non-gif output. "
                        "Install with `pip install imageio-ffmpeg`.")
                matplotlib.rcParams['animation.ffmpeg_path'] = ffmpeg_exe
            writer = animation.FFMpegWriter(
                fps=fps,
                codec='libx264' if ext in ('.mp4', '.mov', '.mkv') else 'libvpx',
                bitrate=bitrate if bitrate is not None else -1,
                extra_args=['-pix_fmt', 'yuv420p'] if ext != '.webm' else None,
            )
            ani.save(str(out_path), writer=writer, dpi=dpi)
        else:
            raise ValueError(f"unsupported output extension {ext!r}; "
                             f"use .gif, .mp4, .webm, .mov, or .mkv")
        plt.close(fig)
        return out_path

    def make_gif(self, output_path, frame_paths, fps=10):
        """
        Assemble PNG frames into a GIF.
        Returns pathlib.Path to output GIF.
        """
        try:
            import imageio
        except ImportError:
            raise ImportError("imageio required for make_gif; pip install imageio")

        output_path = pathlib.Path(output_path)
        frames = [imageio.imread(str(fp)) for fp in frame_paths]
        imageio.mimsave(str(output_path), frames, fps=fps)
        return output_path

    # ── properties ────────────────────────────────────────────────────────────

    @property
    def t(self):
        """Current absolute simulation time."""
        return float(self._state['t'])

    @property
    def step_count(self):
        """Total number of steps taken."""
        return int(self._state['step'])

    def _accessible_area(self):
        """
        Accessible area for packing-fraction calculations.
        For unconstrained / periodic systems: Lx*Ly - obstacle areas.
        For systems with an 'exterior' container polygon (e.g. HopperRegion):
          container polygon area - obstacle areas.
        """
        base = self.Lx * self.Ly
        for obj in self._objects:
            if getattr(obj, '_exclusion', None) == 'exterior':
                poly = obj.region_polygon(0.0)
                if poly is not None:
                    v = np.asarray(poly['vertices'])
                    x, y = v[:, 0], v[:, 1]
                    area = 0.5 * abs(
                        np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))
                    base = min(base, area)
        excl = sum(obj.exclusion_area() for obj in self._objects)
        return base - excl

    @property
    def phi_outer(self):
        """Outer-perimeter packing fraction (obstacle areas excluded from denominator)."""
        from src.epd.initializer import compute_phi_outer
        r_c_arr = self._params['r_c_per_p'].numpy()
        return compute_phi_outer(self._state, self.Lx, self.Ly, r_c_arr,
                                 accessible_area=self._accessible_area())

    @property
    def phi_box(self):
        """Box-area packing fraction (π*R0_mean² * P / (Lx*Ly))."""
        P       = len(self._particles)
        R0_mean = np.mean([p.R0 for p in self._particles])
        return P * np.pi * R0_mean**2 / (self.Lx * self.Ly)

    @property
    def particles(self):
        """Read-only list of CapsuleParticle objects."""
        return list(self._particles)

    @property
    def state(self):
        """Current full state dict (read-only reference)."""
        return self._state

    # ── post-init particle control ────────────────────────────────────────────

    def select_particles(self, criterion, fraction=0.1, **kwargs):
        """
        Select particle indices by a spatial criterion.

        criterion : 'top_fraction'     — top fraction by y_cm
                    'bottom_fraction'  — bottom fraction by y_cm
                    'left_fraction'    — left fraction by x_cm
                    'right_fraction'   — right fraction by x_cm
                    'all'              — all particle indices
                    list of int        — explicit index list (passed through)
        fraction  : float (0,1) — fraction of total particles to select

        Returns list of int particle indices.
        """
        if isinstance(criterion, (list, tuple, np.ndarray)):
            return list(criterion)

        x_cm = self._state['x_cm'].numpy()   # (P, 2)
        P    = x_cm.shape[0]
        n    = max(1, int(round(fraction * P)))

        if criterion == 'all':
            return list(range(P))
        elif criterion == 'top_fraction':
            order = np.argsort(x_cm[:, 1])[::-1]
            return list(order[:n].tolist())
        elif criterion == 'bottom_fraction':
            order = np.argsort(x_cm[:, 1])
            return list(order[:n].tolist())
        elif criterion == 'left_fraction':
            order = np.argsort(x_cm[:, 0])
            return list(order[:n].tolist())
        elif criterion == 'right_fraction':
            order = np.argsort(x_cm[:, 0])[::-1]
            return list(order[:n].tolist())
        else:
            raise ValueError(f"Unknown selection criterion: '{criterion}'")

    def particles_set_motion(self, indices, motion_spec, frozen_shape=True):
        """
        Set driven motion for specific particle indices (called after initialize()).

        indices      : list of int
        motion_spec  : MotionSpec
        frozen_shape : bool — freeze elastic DOFs for driven particles (default True)

        Modifies params in-place. Returns self.
        """
        from src.simulation.tf_sim import set_driven
        traj_rows = [motion_spec.to_traj_row() for _ in indices]
        set_driven(self._params, list(indices), traj_rows,
                   frozen=[frozen_shape] * len(indices))
        return self
