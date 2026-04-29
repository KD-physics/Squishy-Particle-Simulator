# Rigid-Body Decomposition: Aim, Benchmarks, and Metrics

**Status:** Implementation in progress — see `PLAN.md` §Phase 5  
**Code target:** `src/simulation/capsule_shell.py` — `CapsuleParticle` + `CapsuleSim.step()`  
**Benchmark scripts:** `src/validation/rb_benchmark_*.py`

---

## 1. Aim

The current integrator treats every perimeter node as an independent degree of
freedom:

```python
p.f -= alpha_damp * m_node * p.v   # BUG: damps ALL velocity, including rigid translation
a = p.f / m_node
p.v += a * dt
p.x += p.v * dt
```

This is sufficient for the symmetric flat-wall calibration (no net translation
or rotation by construction), but it introduces two physical errors that must be
fixed before any multi-particle simulation is valid:

| Error | Description | Consequence |
|-------|-------------|-------------|
| E1 | Shell inertia $I_\text{shell} = mR_0^2$ instead of disk inertia $I_\text{disk} = \tfrac{1}{2}mR_0^2$ | Particle rotates half as fast as it should |
| E2 | Nodal damping applies to rigid-body velocity, not only elastic velocity | Particle decelerates spuriously in free space; rigid-body momentum is not conserved |

**The fix:** decompose each particle's motion into a rigid part
$(\mathbf{x}_\text{cm}, \theta)$ and an elastic part $\mathbf{u}_i$ at every
timestep.  Integrate the rigid part with undamped Verlet using the correct
solid-disk mass $M = \rho_f \pi R_0^2$ and moment of inertia
$I = \tfrac{1}{2}MR_0^2$.  Apply damping only to $\dot{\mathbf{u}}_i$
(elastic oscillations).

The decomposition at each step:
```
x_cm  = (1/N) Σ x_i                          # centroid
θ     = argmin_φ Σ |x_i − x_cm − R(φ) X_i|² # best-fit rotation (SVD/polar)
u_i   = R(−θ)(x_i − x_cm) − X_i             # elastic displacement in body frame
```

Rigid-body resultants from nodal contact forces:
```
F     = Σ f_i^contact                         # net force on CM
T     = Σ (x_i − x_cm) × f_i^contact         # net torque about CM
```

Integration:
```
# Rigid body (undamped, correct mass/inertia)
v_cm  += (F / M) * dt
x_cm  += v_cm * dt
ω     += (T / I) * dt
θ     += ω * dt

# Elastic (damped, in co-moving frame)
u̇_i  -= alpha_damp * u̇_i * dt              # damp only deformation velocity
u̇_i  += (f_i^elastic / m_node) * dt
u_i   += u̇_i * dt

# Reconstruct perimeter
x_i   = x_cm + R(θ)(X_i + u_i)
```

---

## 2. Benchmarks and Pass Criteria

### Benchmark 0 — Free particle (precursor)

**Setup:** Single particle, no walls, no contacts.  Give it an initial CM
velocity $v_0 = (1, 0)$ and an initial spin $\omega_0$.  Run for 100 particle
transit times ($t = 100 R_0 / v_0$).

**Expected physics:** No forces → constant momentum, constant spin.

| Metric | Pass criterion | Failure mode in buggy code |
|--------|---------------|---------------------------|
| $\|\mathbf{p}_\text{cm}(t) - \mathbf{p}_\text{cm}(0)\| / \|\mathbf{p}_\text{cm}(0)\|$ | $< 10^{-10}$ (machine precision) | Exponential decay from spurious drag |
| $\|\omega(t) - \omega_0\| / \omega_0$ | $< 10^{-10}$ | Decay from spurious rotational drag |
| Shape: $\max_i \|u_i(t)\| / R_0$ | $< 10^{-12}$ | Should be exactly zero — particle started circular |

This benchmark **must pass before any collision benchmark is run**.

---

### Benchmark 1 — Single wall impact

**Setup:** Single particle at rest.  A horizontal wall moves upward at constant
velocity $v_w = 0.1 \sqrt{K_\text{fluid}/(\rho_f R_0)}$ (one-tenth the wave
speed), contacts the particle, compresses it by $\varepsilon_\text{max} = 0.06$
(per-particle strain), then **stops** (wall velocity set to zero).  Run until
the particle has fully separated and the shape has relaxed.

**Three phases to track:**

| Phase | Condition | Expected behavior |
|-------|-----------|-------------------|
| Contact | $F_\text{contact} > 0$ | CM accelerates; elastic energy builds; shape deforms |
| Post-stop, in contact | Wall stopped, particle still touching | Elastic restoring force continues to accelerate CM; particle separates with $v_\text{cm} > v_w$ for elastic case |
| Free flight | $F_\text{contact} = 0$ | CM moves at constant velocity; shape oscillates then damps |

**Metrics:**

| Metric | Pass criterion |
|--------|---------------|
| Momentum in free flight: $\|d\mathbf{p}_\text{cm}/dt\|$ | $< 10^{-10} \cdot M \cdot v_w / R_0$ (machine zero) |
| Coefficient of restitution $e = v_\text{cm,out} / (2 v_w)$ | $0 < e \leq 1$; $e \to 1$ as $\alpha_d \to 0$ |
| Shape relaxation: $\max_i \|u_i\| \to 0$ | Monotonically decreasing after separation |
| Energy closure: $\Delta \text{KE} + \Delta \text{PE} + W_\text{diss} = 0$ | Residual $< 10^{-6}$ (relative) |
| $W_\text{diss}$ attribution: damping only on $\dot{u}_i$, not on $v_\text{cm}$ | $W_\text{diss}^\text{rigid} = 0$ in free flight |

The COR sweep $e(\alpha_d)$ at fixed $q$ and $v_w$ should be a smooth
monotonically decreasing function from $e(0) = 1$ to $e \to 0$ as
$\alpha_d \to \infty$.  This is a physical prediction of the model and
its shape should be reported.

---

### Benchmark 2 — Two-particle head-on collision

**Setup:** Two identical particles approaching head-on with equal and opposite
velocities $\pm v_0$ in the CM frame.  Run through the full compression–rebound
cycle until both particles are separated and shape-relaxed.  Then repeat with
mass ratio 1:4 (different $R_0$ or different $\rho_f$).

**Metrics:**

| Metric | Pass criterion |
|--------|---------------|
| Total momentum $\mathbf{p}_1 + \mathbf{p}_2 = \text{const}$ | Deviation $< 10^{-10} \cdot M v_0$ at all times |
| COR $e = \|\mathbf{v}_\text{rel,out}\| / \|\mathbf{v}_\text{rel,in}\|$ | $e \in (0, 1]$; matches Benchmark 1 for same $\alpha_d$ |
| Equal-mass elastic ($\alpha_d = 0$): particles exchange velocities | $|v_{1,\text{out}} - v_{2,\text{in}}| / v_0 < 10^{-6}$ |
| Energy closure: $\Delta \text{KE} = -W_\text{diss}$ | Residual $< 10^{-6}$ (relative) |
| $W_\text{diss} = \int \alpha_d m_\text{node} \|\dot{\mathbf{u}}_i\|^2 \, dt$ | Equals $\Delta \text{KE}$ to $< 10^{-4}$ (relative) |
| Off-centre impact: net rotation $\omega$ at separation | $|\omega| R_0 / v_0 \ll 1$ (frictionless, nearly circular) |

The equal-mass elastic velocity exchange is the clearest possible signal:
if Newton's 3rd law and momentum conservation are correctly implemented,
this is exact (not approximate).

---

## 3. Dissipation Attribution (energy-budget diagnostic)

The key invariant that distinguishes the correct from the buggy integrator is:

**Correct integrator:**
$$W_\text{diss} = \int_0^T \alpha_d \sum_i m_\text{node}\, \|\dot{\mathbf{u}}_i\|^2\, dt$$

**Buggy integrator:**
$$W_\text{diss}^\text{buggy} = \int_0^T \alpha_d \sum_i m_\text{node}\, \|\mathbf{v}_i\|^2\, dt
= W_\text{diss} + \underbrace{\alpha_d M \|v_\text{cm}\|^2 \cdot T_\text{contact}}_{\text{spurious rigid-body loss}}$$

Track both quantities in the benchmark scripts.  In the buggy code,
$W_\text{diss}^\text{buggy} > W_\text{diss}$ because rigid-body kinetic energy
is incorrectly consumed.  In the correct code, the two are identical.

---

## 4. Notes on Rotation

For frictionless EPD particles, contact forces are normal-only.  For a nearly
circular disk, the moment arm $\mathbf{r}_i \approx R_0\hat{n}_i$ and the force
$\mathbf{f}_i \parallel \hat{n}_i$, so $T_i = \mathbf{r}_i \times \mathbf{f}_i \approx 0$.
Residual torques scale as $\Delta A / A_0 \leq 6\%$ in the calibration regime.

**Expected outcome of off-centre Benchmark 2:** $\omega R_0 / v_0 \lesssim 0.06$.
This is small but non-zero — quantifying it is part of the benchmark.

---

*Document authored: 2026-04-19. Update when implementation is complete and metrics are measured.*
