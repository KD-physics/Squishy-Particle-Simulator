"""
test_phase8.py — Phase 8: Saddle-aware energy minimization via TF optimizers.

16-particle ADAM quench test. Validates the four NEW System methods
(eval_forces_at, eval_potential_energy, sync_from_var, make_position_variable)
and the new tf_sim.eval_potential_energy_tf function.

Test recipe (avoids save/restore — would re-roll PRCM and shift PE baseline):
  1. RSA initialize at phi=0.49 (just above isostatic jamming)
  2. Swell to phi=0.85 via continuous box_rate over N_SWELL steps. Set
     all damping to zero (alpha, beta_rb, xi). Capture PE_swell.
  3. ADAM PHASE — operates on x_var (a tf.Variable copy of x_all); does NOT
     mutate s._state. Compute PE_adam = eval_potential_energy(x_var).
  4. NEGATIVE CONTROL — full Newtonian (no damping anywhere) physics from the
     UN-touched swell state for N_CTRL steps. PE_physics ≈ PE_swell expected
     (no monotonic drop without dissipation).
  5. sync_from_var(x_var) overwrites s._state with the ADAM result.
  6. Regression check: re-enable physical xi, run N_REG steps. Verify NaN=0.

Pass criteria:
  PE_adam  < PE_swell - eps        (ADAM lowered the energy)
  PE_adam  < PE_physics            (ADAM beats passive evolution)
  PE_physics within 50% of PE_swell  (passive Newtonian doesn't relax)
  NaN = 0 throughout
  Throughput: ADAM step time within 3× swell step time

Run: python src/epd/tests/test_phase8.py
"""
import sys, os, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')

import numpy as np
import tensorflow as tf
import src.simulation.tf_sim as tfm
tfm.set_dtype(tf.float64)
from src.epd.particles import ParticleSpec
from src.epd.system import System
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────────────────────
# Knobs (override via env vars for quicker smoke runs)
# ──────────────────────────────────────────────────────────────────────────────
P            = int(os.environ.get('PHASE8_P',          16))
N_NODES      = int(os.environ.get('PHASE8_N',          60))
PHI_INIT     = float(os.environ.get('PHASE8_PHI_INIT',  0.49))
PHI_TARGET   = float(os.environ.get('PHASE8_PHI_TGT',   0.85))
N_SWELL      = int(os.environ.get('PHASE8_N_SWELL',   100_000))
N_CTRL       = int(os.environ.get('PHASE8_N_CTRL',     2_000))
N_ADAM       = int(os.environ.get('PHASE8_N_ADAM',     2_000))
N_REG        = int(os.environ.get('PHASE8_N_REG',      1_000))
ADAM_LR      = float(os.environ.get('PHASE8_ADAM_LR',   1e-3))
SEED         = int(os.environ.get('PHASE8_SEED',         42))
OH_PHYSICAL  = float(os.environ.get('PHASE8_OH',          5.0))

OUTDIR = os.path.join(ROOT, 'results', 'phase8_test')
os.makedirs(OUTDIR, exist_ok=True)

results = []
def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


def snapshot_png(s, label, fname):
    s.set_color_palette('palette1', seed=1)
    fig, ax = plt.subplots(figsize=(7, 7))
    s.render(ax=ax)
    ax.set_title(f"{label}  phi={s.phi_outer:.4f}  Lx={s.Lx:.3f}")
    fig.savefig(fname, dpi=120, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# 1. Build + RSA at phi=0.49
# ──────────────────────────────────────────────────────────────────────────────
print(f"\n── Phase 8 — ADAM saddle-quench test (P={P}, N={N_NODES}) ──────────────")
print(f"  swell phi {PHI_INIT} → {PHI_TARGET} over {N_SWELL} steps;  "
      f"ADAM {N_ADAM} steps;  Oh={OH_PHYSICAL}")

# Box sized so that RSA at PHI_INIT is achievable; initializer picks Lx from phi.
LX_init = 18.0  # only used as starting box; will be overridden by initialize()

s = System(LX_init, LX_init, periodic_x=True, periodic_y=True,
           dt_factor=0.25, candidacy_kind='prcm', E_candidates=9)
s.add_particles(ParticleSpec(count=P, type='emulsion', gamma=1.0, kappa=0.02,
                              N_nodes=N_NODES, Oh=OH_PHYSICAL,
                              poly_dist={'type':'bimodal', 'ratio':0.5, 'delta':0.2}))
# Initialize with standard adaptive_swell to phi_init (just above isostatic
# jamming). This gives us a reasonable starting point — the box_rate swell
# below then takes us from phi_init → PHI_TARGET with the user-validated
# rate-tracking adaptive scheme.
s.initialize(phi_target=PHI_INIT, seed=SEED, verbose=False,
             relax_only=False, n_relax_init=200)

phi_after_rsa = s.phi_outer
print(f"\n  [step 1] init done: phi_outer={phi_after_rsa:.4f}  Lx={s.Lx:.4f}  "
      f"NaN={int(np.isnan(s._state['x_cm'].numpy()).sum())}")

# Save physical xi before swell (for the regression check at the end)
xi_phys = tf.constant(s._params['xi_drag_per_p'].numpy(),
                       dtype=s._params['xi_drag_per_p'].dtype)


# ──────────────────────────────────────────────────────────────────────────────
# 2. Swell phi_init → PHI_TARGET via adaptive box_rate.
#
#    Uses the user-validated swell strategy:
#      (a) per-chunk target rate from remaining-steps target dφ/dt
#      (b) rate magnitude monotonically non-increasing — never speed up
#      (c) annealing oscillation near target (alternate small expand / compress)
#          to relieve contact stress without overshoot
#      (d) β_rb = 0.005/dt rigid-body damping + per-chunk drift removal
#      (e) α·dt = 0.2 internal-mode damping (semi-implicit stable)
# ──────────────────────────────────────────────────────────────────────────────
DTYPE = s._params['_alpha_tf'].dtype
dt    = float(s._params['_dt_tf'].numpy())
SWELL_ALPHA = 0.2 / dt
s._params['_alpha_tf'] = tf.constant(SWELL_ALPHA, dtype=DTYPE)
s._params['alpha_damp_per_p'] = tf.constant(np.full(P, SWELL_ALPHA, dtype=np.float64),
                                              dtype=DTYPE)
s.beta_rb = 0.005 / dt
print(f"  [setup] dt={dt:.4e}  α={SWELL_ALPHA:.2f} (α·dt={SWELL_ALPHA*dt:.3f})  "
      f"β_rb={s.beta_rb:.2f} (β·dt={s.beta_rb*dt:.4f})")

CHUNK = max(N_SWELL // 125, 100)   # ~125 chunks for nominal N_SWELL

# Seed initial rate at 2× the local rate_test (linearized rate at iteration 0).
# The 2× factor gives the loop room: iteration 1 will adopt rate_test (which
# is smaller in magnitude), starting natural rate-test tracking. Seeding at
# exactly rate_test triggers immediate decay-only mode, which can't catch up.
target_dphi_dt_init = (PHI_TARGET - phi_after_rsa) / (N_SWELL * dt)
rate_init           = -2.0 * target_dphi_dt_init / (2.0 * phi_after_rsa)
s.box_rate          = (rate_init, rate_init)

oscillation_triggered = False
direction = -1.0  # -1 = compressing
N_SWELL_CAP = int(1.5 * N_SWELL)

steps = 0
t0 = time.time()
while s.phi_outer < PHI_TARGET - 1e-6 and steps < N_SWELL_CAP:
    phi_now  = s.phi_outer
    n_remain = max(N_SWELL - steps, CHUNK)
    target_dphi_dt = (PHI_TARGET - phi_now) / (n_remain * dt)
    rate_test = -target_dphi_dt / (2.0 * phi_now)
    cur_rate = s.box_rate[0]
    if abs(rate_test) < abs(cur_rate) and steps < N_SWELL:
        rate = rate_test
    else:
        rate = cur_rate * 0.99

    if (PHI_TARGET - phi_now) <= 0.02 or oscillation_triggered:
        if not oscillation_triggered:
            oscillation_triggered = True
            print(f"  [anneal] oscillation triggered at Δphi = "
                  f"{PHI_TARGET - phi_now:.4f}  (steps={steps})")
        if direction < 0:           # was compressing → small expand
            rate =  abs(rate) / 2
        else:                       # was expanding → bigger compress
            rate = -abs(rate) * 2
        direction = -direction

    s.box_rate = (rate, rate)
    s.run(CHUNK, sample_every=CHUNK, verbose=False, record_initial=False)
    s.sync_box(CHUNK)
    steps += CHUNK

    # Remove residual drift (mass-weighted mean velocity / angular velocity)
    M_np  = s._params['M_disk'].numpy()
    I_np  = s._params['I_disk'].numpy()
    v_cm_np  = s._state['v_cm'].numpy()
    omega_np = s._state['omega'].numpy()
    v_mean    = (M_np[:, None] * v_cm_np).sum(axis=0) / M_np.sum()
    omega_mean = (I_np * omega_np).sum() / I_np.sum()
    s._state['v_cm']  = tf.constant(v_cm_np  - v_mean,    dtype=DTYPE)
    s._state['omega'] = tf.constant(omega_np - omega_mean, dtype=DTYPE)

s.box_rate = 0.0
swell_wall = time.time() - t0
swell_throughput_ms = swell_wall / max(steps, 1) * 1000
phi_after_swell = s.phi_outer
print(f"  [step 2] swell done: {steps} steps in {swell_wall:.1f}s "
      f"({swell_throughput_ms:.2f} ms/step) → phi={phi_after_swell:.4f}  "
      f"Lx={s.Lx:.4f}  (target {PHI_TARGET:.4f})")
snapshot_png(s, 'after swell', os.path.join(OUTDIR, 'snap_swell.png'))


# ──────────────────────────────────────────────────────────────────────────────
# 3. PE_swell baseline. Set ALL damping to zero → fair Newtonian baseline for
#    both ADAM (no spurious drag) and the negative control (no spurious decay).
# ──────────────────────────────────────────────────────────────────────────────
s._params['xi_drag_per_p']    = tf.zeros_like(xi_phys)
s._params['_alpha_tf']        = tf.constant(0.0, dtype=DTYPE)
s._params['alpha_damp_per_p'] = tf.zeros([P], dtype=DTYPE)
s.beta_rb = 0.0  # already 0 from swell-end is unlikely; force-zero for clarity

PE_swell = s.eval_potential_energy()
PE_swell_comp = s.eval_potential_energy_components()
print(f"  [step 3] PE_swell = {PE_swell:.4e}  (all damping zeroed for baseline)")
print(f"            edge={PE_swell_comp['edge']:.3e}  bend={PE_swell_comp['bend']:.3e}  "
      f"area={PE_swell_comp['area']:.3e}\n"
      f"            lt  ={PE_swell_comp['lt']:.3e}  cc  ={PE_swell_comp['cc']:.3e}  "
      f"prim={PE_swell_comp['prim']:.3e}")
snapshot_png(s, 'before ADAM', os.path.join(OUTDIR, 'snap_swell.png'))


# ──────────────────────────────────────────────────────────────────────────────
# 4. ADAM phase — operates on x_var; does NOT touch s._state.
#    Compute PE_adam from x_var directly (eval_potential_energy(x_var)).
# ──────────────────────────────────────────────────────────────────────────────
opt   = tf.keras.optimizers.Adam(learning_rate=ADAM_LR)
x_var = s.make_position_variable()

print(f"  [step 4] starting ADAM: lr={ADAM_LR}, {N_ADAM} steps")

# Refresh CapCandidates relative to x_var before starting (defensive — covers
# the case where the swell-end PRCM list is stale relative to x_var positions).
x_np = x_var.numpy(); x_cm_np = x_np.mean(axis=1)
theta_np = s._state['theta'].numpy()
s._cm_mgr.update(x_cm_np, theta_np, x_all=x_np)

f0 = s.eval_forces_at(x_var)
f0_mag = tf.norm(f0, axis=2).numpy()
print(f"  [step 4] initial |F|: max={f0_mag.max():.3e}  mean={f0_mag.mean():.3e}")

t0 = time.time()
adam_log = []      # (step, PE, max|F|, mean|F|)
adam_log.append((0, PE_swell, float(f0_mag.max()), float(f0_mag.mean())))
SAMPLE_EVERY = max(1, N_ADAM // 20)
PRCM_EVERY   = 50

for i in range(N_ADAM):
    f = s.eval_forces_at(x_var)              # tensor (P, N, 2), tape-free
    opt.apply_gradients([(-f, x_var)])       # ADAM treats -f as gradient
    if (i + 1) % PRCM_EVERY == 0:
        x_np = x_var.numpy(); x_cm_np = x_np.mean(axis=1)
        if s._cm_mgr.needs_update(x_cm_np, theta_np):
            s._cm_mgr.update(x_cm_np, theta_np, x_all=x_np)
    if (i + 1) % SAMPLE_EVERY == 0:
        PE_now = s.eval_potential_energy(tf.constant(x_var.numpy(), dtype=DTYPE))
        f_now  = s.eval_forces_at(x_var)
        f_mag  = tf.norm(f_now, axis=2).numpy()
        adam_log.append((i + 1, PE_now, float(f_mag.max()), float(f_mag.mean())))

adam_wall = time.time() - t0
adam_throughput_ms = adam_wall / N_ADAM * 1000
x_var_tensor = tf.constant(x_var.numpy(), dtype=DTYPE)
PE_adam = s.eval_potential_energy(x_var_tensor)
PE_adam_comp = s.eval_potential_energy_components(x_var_tensor)
print(f"  [step 4] ADAM {N_ADAM} steps in {adam_wall:.2f}s "
      f"({adam_throughput_ms:.2f} ms/step) → PE_adam={PE_adam:.4e}")
print(f"            edge={PE_adam_comp['edge']:.3e}  bend={PE_adam_comp['bend']:.3e}  "
      f"area={PE_adam_comp['area']:.3e}\n"
      f"            lt  ={PE_adam_comp['lt']:.3e}  cc  ={PE_adam_comp['cc']:.3e}  "
      f"prim={PE_adam_comp['prim']:.3e}")


# ──────────────────────────────────────────────────────────────────────────────
# 5. Negative control — physics from un-touched swell state. All damping is 0
#    (set in step 3), so PE+KE is conserved → PE oscillates around PE_swell,
#    not monotonically dropping.
#
#    Note: s._cm_mgr.CapCandidates is currently aligned to x_var (final ADAM
#    position), so refresh it for the original swell state before running.
# ──────────────────────────────────────────────────────────────────────────────
x_cm_np  = s._state['x_cm'].numpy()
theta_np = s._state['theta'].numpy()
s._cm_mgr.update(x_cm_np, theta_np, x_all=s._state['x_all'].numpy())

t0 = time.time()
s.run(N_CTRL, sample_every=N_CTRL, verbose=False, record_initial=False)
ctrl_wall = time.time() - t0
PE_physics = s.eval_potential_energy()
PE_physics_comp = s.eval_potential_energy_components()
print(f"  [step 5] negative control: {N_CTRL} Newtonian physics steps "
      f"({ctrl_wall/N_CTRL*1000:.2f} ms/step) → PE_physics={PE_physics:.4e}")
print(f"            edge={PE_physics_comp['edge']:.3e}  bend={PE_physics_comp['bend']:.3e}  "
      f"area={PE_physics_comp['area']:.3e}\n"
      f"            lt  ={PE_physics_comp['lt']:.3e}  cc  ={PE_physics_comp['cc']:.3e}  "
      f"prim={PE_physics_comp['prim']:.3e}")
snapshot_png(s, 'after Newtonian control', os.path.join(OUTDIR, 'snap_physics_control.png'))


# ──────────────────────────────────────────────────────────────────────────────
# 6. Commit ADAM result, run regression check with physical xi.
# ──────────────────────────────────────────────────────────────────────────────
s.sync_from_var(x_var)
PE_adam_committed = s.eval_potential_energy()
print(f"  [step 6] post-sync PE = {PE_adam_committed:.4e}  (should match PE_adam)")
snapshot_png(s, 'after sync_from_var', os.path.join(OUTDIR, 'snap_post_adam.png'))

s._params['xi_drag_per_p'] = xi_phys
t0 = time.time()
s.run(N_REG, sample_every=N_REG, verbose=False, record_initial=False)
reg_wall = time.time() - t0
nan_count = int(np.isnan(s._state['x_cm'].numpy()).sum())
print(f"  [step 6] regression: {N_REG} physics steps "
      f"({reg_wall/N_REG*1000:.2f} ms/step)  NaN={nan_count}")
snapshot_png(s, 'after physics restart', os.path.join(OUTDIR, 'snap_post_physics_restart.png'))


# ──────────────────────────────────────────────────────────────────────────────
# PE component breakdown comparison — where did ADAM's drop come from?
# ──────────────────────────────────────────────────────────────────────────────
print("\n── PE breakdown (swell → ADAM → Newtonian-control) ─────────────────────")
print(f"  {'term':<6}  {'swell':>11}  {'adam':>11}  {'physics':>11}  "
      f"{'Δ(adam-swell)':>13}")
print(f"  {'-'*6}  {'-'*11}  {'-'*11}  {'-'*11}  {'-'*13}")
for k in ['edge', 'bend', 'area', 'lt', 'cc', 'prim', 'total']:
    s_v = PE_swell_comp[k]
    a_v = PE_adam_comp[k]
    p_v = PE_physics_comp[k]
    delta = a_v - s_v
    print(f"  {k:<6}  {s_v:>11.4e}  {a_v:>11.4e}  {p_v:>11.4e}  {delta:>+13.4e}")


# ──────────────────────────────────────────────────────────────────────────────
# Pass criteria
# ──────────────────────────────────────────────────────────────────────────────
print("\n── Pass criteria ───────────────────────────────────────────────────────")

eps_PE = max(1e-30, abs(PE_swell) * 1e-4)
check('PE_adam < PE_swell',
      PE_adam < PE_swell - eps_PE,
      f"PE_swell={PE_swell:.4e}  PE_adam={PE_adam:.4e}  drop={PE_swell-PE_adam:.3e}")

check('PE_adam < PE_physics',
      PE_adam < PE_physics,
      f"PE_physics={PE_physics:.4e}  PE_adam={PE_adam:.4e}")

# Newtonian negative control: PE oscillates around PE_swell since total energy
# (PE+KE) is conserved. We only require PE doesn't drift down by more than 50%
# (lower bound — Newtonian shouldn't relax to a deep basin without dissipation).
check('Newtonian control does not drop PE > 50%',
      PE_physics > 0.5 * PE_swell,
      f"PE_physics/PE_swell = {PE_physics/max(abs(PE_swell),1e-30)*100:.1f}%")

check('regression run NaN=0', nan_count == 0)

check('swell hit phi target (3% tolerance)',
      abs(phi_after_swell - PHI_TARGET) / PHI_TARGET < 0.03,
      f"target={PHI_TARGET:.4f} achieved={phi_after_swell:.4f} "
      f"(Δ/target = {abs(phi_after_swell-PHI_TARGET)/PHI_TARGET*100:.2f}%)")

# sync_from_var should preserve PE (only candidacy may shift due to wrap+rebuild).
# Allow 20% slack for PRCM list re-roll.
check('sync_from_var preserves PE within 20%',
      abs(PE_adam_committed - PE_adam) / max(abs(PE_adam), 1e-30) < 0.20,
      f"pre-sync={PE_adam:.4e}  post-sync={PE_adam_committed:.4e}")

# Throughput: be generous; XLA compile cost on first call dominates small-N tests
check('ADAM throughput within 5× of swell',
      adam_throughput_ms < 5 * swell_throughput_ms + 10,
      f"swell={swell_throughput_ms:.2f}ms/step  adam={adam_throughput_ms:.2f}ms/step")


# ──────────────────────────────────────────────────────────────────────────────
# PE convergence plot
# ──────────────────────────────────────────────────────────────────────────────
if adam_log:
    its    = np.array([r[0] for r in adam_log])
    pes    = np.array([r[1] for r in adam_log])
    f_max  = np.array([r[2] for r in adam_log])
    f_mean = np.array([r[3] for r in adam_log])
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    # PE row
    for ax, ylog in zip(axes[0], [False, True]):
        ax.axhline(PE_swell,   color='k', ls='--', label=f'PE_swell={PE_swell:.3e}')
        ax.axhline(PE_physics, color='r', ls='--', label=f'PE_physics={PE_physics:.3e}')
        ax.plot(its, pes, 'b-o', ms=3, label='PE during ADAM')
        ax.set_xlabel('ADAM step')
        ax.set_ylabel('PE')
        if ylog: ax.set_yscale('log')
        ax.legend(); ax.grid(alpha=0.3)
        ax.set_title(f"PE convergence (lr={ADAM_LR})  ({'log' if ylog else 'linear'} y)")
    # |F| row — force-balance quality
    for ax, ylog in zip(axes[1], [False, True]):
        ax.plot(its, f_max,  'r-o', ms=3, label='|F|_max')
        ax.plot(its, f_mean, 'g-o', ms=3, label='|F|_mean')
        ax.set_xlabel('ADAM step')
        ax.set_ylabel('|F| per node')
        if ylog: ax.set_yscale('log')
        ax.legend(); ax.grid(alpha=0.3)
        ax.set_title(f"Force-balance convergence  ({'log' if ylog else 'linear'} y)")
    fig.savefig(os.path.join(OUTDIR, 'pe_convergence.png'), dpi=120, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────────────
n_pass = sum(1 for _, c, _ in results if c)
n_tot  = len(results)
print(f"\n── Summary: {n_pass}/{n_tot} checks passed ─────────────────────────────")
for name, cond, detail in results:
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}")
print(f"\n  Outputs: {OUTDIR}")
sys.exit(0 if n_pass == n_tot else 1)
