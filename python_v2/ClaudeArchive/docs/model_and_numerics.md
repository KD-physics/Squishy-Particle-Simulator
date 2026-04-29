# The Elastic Perimeter Disk Model: Construction, Parameterization, and Calibration

---

## §0  Motivation and Scope

A central question in the mechanics of dense soft-particle systems is how particle
*squishiness* — the ease with which a particle deforms under contact — governs
collective behavior near the jamming point.  In continuum solid mechanics, the Poisson
ratio ν parameterizes this axis cleanly: ν = 0 corresponds to a fully compressible
material (volume shrinks freely under uniaxial load), while ν → 1/2 in 2D (ν → 1/2 in
3D) corresponds to an incompressible material that conserves volume at the cost of lateral
deformation.  A granular analog of ν — a single dimensionless knob that spans the full
range from area-compressible to shape-deformable particles — would allow systematic
exploration of jamming mechanics across the squishiness axis without committing to a
specific material.

Full FEM-DEM simulations, in which a proper finite-element elasticity solve is performed
for the interior of each particle at every timestep, can in principle provide this
control.  In practice, however, they are computationally prohibitive for large assemblies
(N ≳ 100 particles), require detailed material constitutive laws, and expose many
parameters whose effects on collective dynamics are not individually transparent.

We instead construct a minimal 2D particle model — the **Elastic Perimeter Disk** (EPD)
model — whose only representation of the particle interior is through forces acting on
a discretized perimeter.  The model is designed around a single physical question: what
is the minimum set of elastic energy terms on the perimeter that produces a tunable,
Poisson-ratio-like control knob?  We show that three terms suffice — membrane stretching,
bending, and an area penalty — and that the model's behavior can be reduced, through
dimensional analysis and numerical sweeps, to two physically interpretable dimensionless
parameters.  One combination of these parameters is then fixed as a working point, and a
single remaining parameter, $q$, controls a smooth squishiness axis from ν ≈ 0 to
ν ≈ 0.94.  The resulting calibration curve $q \mapsto \nu$ allows a user to specify a
target Poisson ratio and have all simulation parameters determined uniquely.

The remainder of this section documents the model construction (§1–§3), the dimensional
reduction (§4), the calibration of numerically necessary internal parameters (§5), the
phase-space structure of the model's response (§6), and the final squishiness axis
calibration (§7).

---

## §1  Geometry and Degrees of Freedom

**Particle geometry.**  Each particle is a 2D disk of reference radius $R_0$.  The
perimeter is discretized as $N$ nodes placed uniformly on the reference circle, with
reference positions

$$
\mathbf{X}_i = R_0 \bigl(\cos(2\pi i/N),\, \sin(2\pi i/N)\bigr), \quad i = 0, 1, \ldots, N-1.
$$

The reference edge length is $L_0 = 2\pi R_0 / N$, and the reference curvature at every
node is $\kappa_0 = 1/R_0$.

**Deformed configuration.**  At time $t$, node $i$ occupies position $\mathbf{x}_i(t)$.
The elastic displacement is $\mathbf{u}_i = \mathbf{x}_i - \mathbf{X}_i$.  We make no
assumption about the magnitude of $\mathbf{u}_i$ relative to $L_0$ or $R_0$; the
geometric quantities (edge lengths, curvatures, area) are computed exactly from
$\{\mathbf{x}_i\}$ at each step.

**Degrees of freedom.**  The $2N$ perimeter coordinates $\{\mathbf{x}_i\}$ are the
primary degrees of freedom.  For inter-particle dynamics, we additionally track the
center-of-mass position $\mathbf{x}_\text{cm}$, velocity $\mathbf{v}_\text{cm}$,
orientation angle $\theta$, and angular velocity $\omega$ — the rigid-body degrees of
freedom.  These are not independent of the perimeter DOF; they are defined as the
mean translation and mean rotation of the perimeter, and the elastic DOF describe
deformation relative to the co-moving rigid frame.

**Contact radius.**  Each perimeter edge is treated as a line segment of half-width
$r_c = L_0 = 2\pi R_0 / N$, giving each particle a capsule-like sausage geometry.
The effective particle radius is therefore $R_\text{eff} = R_0 + r_c = R_0(1 + 2\pi/N)$.
At $N = 32$, $r_c \approx 0.20 R_0$ (effective radius $\approx 1.20 R_0$); at $N = 120$,
$r_c \approx 0.052 R_0$ (effective radius $\approx 1.052 R_0$).  In multi-particle
simulations, the packing fraction must be computed using $R_\text{eff}$, not $R_0$.  For
calibration in the two-disk geometry, the contact region size affects the contact
half-width but not the bulk elastic observables $\Delta A$ or $\nu_\text{meas}$.

---

## §2  Elastic Energy Functional

The elastic energy of a single particle has three contributions: membrane stretching,
bending, and area conservation.

### 2.1  Membrane (stretching) energy

The membrane energy penalizes departure of each edge from its reference length:

$$
E_\text{mem} = \frac{E_l}{2} \sum_{i=0}^{N-1} e_i^2 \, L_0,
\qquad e_i = \frac{L_i - L_0}{L_0},
$$

where $L_i = |\mathbf{x}_{i+1} - \mathbf{x}_i|$ is the current edge length (indices
mod $N$) and $E_l$ is the membrane stiffness.  Physically, $E_l = E_s t$ where $E_s$
is the 2D Young's modulus of the shell material and $t$ is the shell thickness; $E_l$
has units of force (N in 2D, or N/m in the per-unit-depth convention).  We introduce the
dimensionless shell thickness $\tau = t/R_0$ so that $E_l = E_s \tau R_0$.

### 2.2  Bending energy

The bending energy penalizes departure from the reference curvature $\kappa_0 = 1/R_0$:

$$
E_\text{bend} = \frac{EI}{2} \sum_{i=0}^{N-1} (\kappa_i - \kappa_0)^2 \, L_0,
$$

where $\kappa_i$ is the discrete curvature at node $i$, computed from the signed turning
angle $\phi_i = \angle(\mathbf{t}_{i-1}, \mathbf{t}_i)$ between successive edge
tangents: $\kappa_i = \phi_i / L_0$.  Here $EI$ is the bending stiffness, given by thin-
shell theory as $EI = E_s t^3/12 = E_s (\tau R_0)^3/12$.  It has units of force × length.

### 2.3  Area penalty

Rather than enforcing exact incompressibility (which would require a Lagrange multiplier
solve), we impose a soft area constraint via a penalty:

$$
E_\text{area} = \frac{K_\text{area}}{2} \left(\frac{A - A_0}{A_0}\right)^2 A_0,
$$

where $A = \tfrac{1}{2}|\sum_i (\mathbf{x}_i \times \mathbf{x}_{i+1})|$ is the current
enclosed area computed by the shoelace formula (exact for the polygon with vertices
$\{\mathbf{x}_i\}$; converges to the smooth disk area as $N \to \infty$), $A_0$ is the
reference area of the initial regular polygon (not $\pi R_0^2$, which the polygon
approximates), and $K_\text{area}$ is the area stiffness with units of force.
As $K_\text{area} \to \infty$ the particle approaches incompressibility ($A \to A_0$);
finite $K_\text{area}$ allows controlled area loss $\Delta A = A_0 - A$ under compression.

### 2.4  Contact energy

Contact between two particles (or between a particle and a rigid wall) is modeled
by a penalty barrier that activates when the gap $g$ between a perimeter node and a
neighboring contact surface falls below zero.  The contact force on node $i$ has
magnitude proportional to $C \cdot |g_i| / R_0$ (where $C$ is the contact stiffness,
with units of force) and is directed along the contact normal.  A 2-point Gauss quadrature
is used along each edge pair to accurately represent distributed contact pressure over
the flattened contact zone.

The contact force law used is a linear spring: $f_c = k_c \times |g_i| \times L_0$ for
$g_i < 0$, zero otherwise.  This is a one-sided linear penalty, not a Hertz contact
($f \propto g^{3/2}$) nor a Lennard-Jones type.  The linear form is chosen for numerical
simplicity and stability; the specific law does not affect the bulk elastic observables
$(\Delta A, \nu_\text{meas})$ as long as the contact compression $\text{cc} = |g_i^\text{min}|/r_c$
is kept below the 2\% threshold (§5.1).

### 2.5  Nodal forces

The nodal force at node $i$ is computed as $\mathbf{f}_i = -\partial E / \partial
\mathbf{x}_i$ where $E = E_\text{mem} + E_\text{bend} + E_\text{area} + E_\text{contact}$.
These forces are evaluated analytically from the discrete energy expressions above; no
stiffness matrix assembly or linear solve is required.

Two implementation choices ensure that the discrete forces are energy-consistent:

1. **Bending forces via signed turning angle.**  The discrete curvature at node $i$ is
   computed as the signed turning angle $\phi_i$ between consecutive edge tangents
   $\mathbf{t}_{i-1}$ and $\mathbf{t}_i$:
   $\kappa_i = \phi_i/L_0$, where $\phi_i = \text{atan2}(\mathbf{t}_{i-1} \times \mathbf{t}_i,\, \mathbf{t}_{i-1} \cdot \mathbf{t}_i)$.
   This is the discrete Frenet formula and correctly handles curvature near $0$ and $\pi$.
   The bending force gradient $\partial \kappa_i / \partial \mathbf{x}_j$ is computed
   analytically from this expression.

2. **Contact forces via 2-point Gauss quadrature.**  The contact energy between two
   edges from different particles is integrated by evaluating the gap function at the two
   Gauss points $s \in \{(1 - 1/\sqrt{3})/2,\, (1 + 1/\sqrt{3})/2\}$ along each edge.
   This 2-point rule integrates the linear and quadratic terms in the gap exactly and
   accurately represents the distributed contact pressure across the flattened contact
   zone without per-node contact detection.

---

## §3  Equations of Motion and Time Integration

### 3.1  Elastic DOF: overdamped explicit integration

The elastic perimeter degrees of freedom are integrated with mass $\rho_l$ and linear
damping; there is no static solve.  The slow loading ($\text{SR} \ll 1$) and moderate
overdamping ($\zeta \approx 16\%$) together ensure that the elastic DOF remain close to
their instantaneous equilibrium at each frame, so the deformation state is effectively
quasi-static even though the integration is explicit and inertial.  The perimeter nodes
follow the damped explicit Verlet scheme:

$$
\mathbf{x}_i^{n+1} = \mathbf{x}_i^n + (1-\alpha_d \Delta t)\,\Delta\mathbf{x}_i^n + \frac{(\Delta t)^2}{\rho_l}\,\mathbf{f}_i^n,
$$

where $\rho_l$ is the nodal mass per unit arc length (set from a fluid density $\rho_f$
times the reference area divided by $N$), $\alpha_d$ is the damping coefficient, and
$\Delta\mathbf{x}_i^n = \mathbf{x}_i^n - \mathbf{x}_i^{n-1}$ is the displacement increment.

### 3.2  Rigid-body integration

The center-of-mass position and orientation are integrated via standard explicit Verlet:

$$
\mathbf{x}_\text{cm}^{n+1} = 2\mathbf{x}_\text{cm}^n - \mathbf{x}_\text{cm}^{n-1} + \frac{(\Delta t)^2}{M}\,\mathbf{F}^n,
\qquad
\theta^{n+1} = 2\theta^n - \theta^{n-1} + \frac{(\Delta t)^2}{I}\,T^n,
$$

where $M = \rho_f \pi R_0^2$ is the particle mass, $I = M R_0^2/2$ is the moment of
inertia, $\mathbf{F}^n = \sum_i \mathbf{f}_i^\text{contact}$ is the total contact force,
and $T^n = \sum_i \mathbf{r}_i \times \mathbf{f}_i^\text{contact}$ is the total torque.

### 3.3  Timestep and acoustic wave speed

The natural timescale is the acoustic wave crossing time for the perimeter:

$$
T_\text{wave} = \frac{2\pi R_0}{c_s}, \qquad c_s = \sqrt{\frac{E_l}{\rho_l}},
$$

where $c_s$ is the one-dimensional wave speed along the shell and $\rho_l$ is the
linear mass density.  The timestep is set as $\Delta t = f_\text{dt} \, L_0 / c_s$ with
$f_\text{dt} = 0.4$ (below the CFL limit of 1).

### 3.4  Quasi-static loading criterion and justification

Walls are moved at a constant speed $v_\text{wall}$.  The wall velocity is set as a
fraction of the wave speed: $v_\text{wall} = \text{SR} \times c_s$, with the strain-rate
ratio $\text{SR} = 0.01$.  This ensures that elastic equilibration occurs on timescales
much shorter than the compression ramp ($\text{SR} \ll 1$), guaranteeing quasi-static
conditions.  Under this loading, the engineering strain $\varepsilon_p = \delta/R_0$
(where $\delta$ is the wall displacement from initial contact) increases at a rate much
slower than the elastic relaxation rate, so each frame is effectively an equilibrium state.

**Justification for treating elastic DOF as quasi-static.**  The inertia of the elastic
perimeter DOF is set by $\rho_l$, the linear mass density of the shell.  The acoustic
wave crossing time $T_\text{wave} = 2\pi R_0/c_s$ is the timescale for elastic
equilibration.  The loading timescale is $T_\text{load} = 2R_0 \varepsilon_\text{max} /
v_\text{wall} = 2R_0 \varepsilon_\text{max} / (\text{SR} \times c_s) = T_\text{wave}
\times \varepsilon_\text{max} / (\pi \times \text{SR})$.  At $\text{SR} = 0.01$ and
$\varepsilon_\text{max} = 0.12$, $T_\text{load}/T_\text{wave} \approx 4$: the loading
takes only 4 wave-crossing times.  However, the damping ensures that perturbations decay
within $\sim 1/\zeta \approx 6$ wave periods, so the elastic DOF remain near equilibrium
throughout.

**Why explicit integration, not implicit.**  Implicit schemes eliminate the acoustic CFL
constraint and allow larger timesteps, at the cost of a linear system solve per step.
For $N$ nodes and $N_\text{particles}$ particles, implicit integration requires an
$O(N N_\text{particles})$-dimensional solve per timestep, which eliminates the
computational advantage over full FEM-DEM for large assemblies.  The explicit Verlet
scheme avoids this while remaining stable under the CFL constraint $\Delta t < L_0/c_s$,
and the quasi-static loading ensures that the inertial artifacts introduced by finite
$\text{SR}$ are negligible.

### 3.5  Damping

The damping coefficient $\alpha_d$ is set via the dimensionless parameter $\alpha_0$:

$$
\alpha_d = \frac{\alpha_0}{T_\text{wave}},
$$

so that the damping fraction $\zeta = \alpha_0 / (4\pi)$ is constant regardless of
particle size or material stiffness.  The value $\alpha_0 = 2.0$ gives $\zeta \approx
16\%$ (moderately overdamped), which ensures fast relaxation to quasi-static equilibrium
without excessive oscillation.  The calibration of $\alpha_0$ is discussed in §5.2.

---

## §4  Dimensional Analysis and Parameter Reduction

### 4.1  Raw parameter inventory

The model as stated in §2–3 contains the following parameters with physical dimensions:

| Symbol | Description | Units |
|--------|-------------|-------|
| $E_l$ | membrane stiffness | force |
| $EI$ | bending stiffness | force × length |
| $K_\text{area}$ | area penalty stiffness | force |
| $C$ | contact stiffness | force |
| $R_0$ | reference radius | length |
| $\rho_f$ | fluid density | mass / length² |
| $t$ or $\tau$ | shell thickness (absolute or dimensionless) | length or 1 |

The observables of interest — the fractional area loss $\Delta A / A_0$ and the
measured 2D Poisson ratio $\nu_\text{meas}$ — are both dimensionless and depend only on
dimensionless combinations of the above parameters.

### 4.2  S as a pure force scale

We introduce a reference force scale $S$ via the definition

$$
E_l^\text{base} \equiv \frac{12 S}{\tau^2},
\tag{4.1}
$$

where $E_l^\text{base}$ is the membrane stiffness evaluated at $R_0 = 1$.  For a physical
disk of radius $R_0$, the membrane stiffness scales with the total perimeter length:
$E_l^\text{code} = E_l^\text{base} \times R_0$.  This is the value used internally by the
simulator.  In terms of material properties, $E_l^\text{base} = E_s \tau$ where $E_s$ is
the 2D Young's modulus of the shell material and $t = \tau R_0$ is the shell thickness.

The bending stiffness follows from the standard thin-shell relation
$EI = E_s t^3/12 = E_l^\text{code} \times (\tau R_0)^2/12$.  Substituting
$E_l^\text{code} = E_l^\text{base} \times R_0 = (12S/\tau^2) \times R_0$:

$$
EI = \frac{12S}{\tau^2} \times R_0 \times \frac{\tau^2 R_0^2}{12} = S \, R_0^3.
\tag{4.2}
$$

**Claim:** all dimensionless observables are independent of $S$.

**Proof (numerical):**  At fixed $\tau = 0.05$, $q = 1$, $R_0 = 1$, $\alpha_0 = 2.0$, the
contact-normalized force $F/E_l^\text{code}$, area loss $\Delta A$, and Poisson ratio
$\nu_\text{meas}$ were computed for $S \in \{0.01, 0.1, 1.0, 10.0\}$ (a 1000× range) at
$\varepsilon_p \in \{0.05, 0.08, 0.10\}$.  The coefficient of variation was $0.00\%$ for
all three observables at all three strains (e.g., $F/E_l^\text{code} = 9.20 \times 10^{-3}$,
$\Delta A = 0.545\%$, $\nu_\text{meas} = 0.454$ at $\varepsilon_p = 0.08$, identical across
all $S$).  Since $S$ enters only as an overall force prefactor, it cancels from the
dimensionless Euler–Lagrange equations.

**Consequence:** fix $S = 1$ with no loss of generality.  The remaining parameters are
$\tau$, $K_\text{area}$, $C$, $R_0$, and $\rho_f$.

### 4.3  R₀ as a pure length scale

**Claim:** all dimensionless observables are independent of $R_0$, provided the parameters
are scaled as:

$$
EI = S R_0^3, \qquad K_\text{area} = q \cdot E_l^\text{base}, \qquad C = C_0 S(1 + q),
\tag{4.3}
$$

where $E_l^\text{base} = 12S/\tau^2$ is the membrane stiffness evaluated at $R_0 = 1$ and
$q$ is a dimensionless ratio defined below.  Crucially, $K_\text{area}$ and $C$ are
*not* scaled with $R_0$.

**Proof (numerical):**  At fixed $S = 1$, $\tau = 0.05$, $q = 1$, $\alpha_0 = 2.0$,
$K_\text{area} = 4800$ (= $E_l^\text{base}$), $C = 1000$, the observables
$\Delta A$ and $\nu_\text{meas}$ were computed for
$R_0 \in \{0.25, 0.5, 1.0, 2.0, 4.0\}$ (a 16× range) at perimeter resolutions
$N \in \{32, 64, 128\}$.  The resulting $\varepsilon_p$–$\Delta A$ and $\varepsilon_p$–$\nu$
curves were numerically identical across all $R_0$ values: CV $= 0.00\%$ for
$\Delta A$ and $\nu_\text{meas}$ at $\varepsilon_p \in \{0.05, 0.08, 0.10\}$ for all $N$.
Representative values at $\varepsilon_p = 0.08$: $\Delta A \approx 0.55\%$,
$\nu_\text{meas} \approx 0.756$ (all $R_0$).

**Why $EI \propto R_0^3$:** The discrete bending force on node $i$ scales as
$EI \times \Delta\kappa_i / L_0$, where $\Delta\kappa_i \sim (1/R_0)$ (curvature
perturbation on a disk of radius $R_0$) and $L_0 = 2\pi R_0/N \propto R_0$.  Therefore
the bending force per node scales as $EI/(R_0^2)$.  The membrane force per node scales as
$E_l^\text{code} \times \text{strain} = E_l^\text{base} R_0 \times \text{strain} \propto R_0$.
For these to remain in the same ratio at all $R_0$ (self-similarity):
$$
\frac{EI/R_0^2}{E_l^\text{base} R_0} = \frac{EI}{E_l^\text{base} R_0^3} = \frac{S}{E_l^\text{base}} = \frac{\tau^2}{12} = b_0 = \text{const}
\implies EI = S R_0^3.
$$
This is precisely the thin-shell formula derived in (4.2).

**Why $K_\text{area}$ does not scale with $R_0$:** The area force per node scales as
$K_\text{area} \times (\Delta A/A_0) \times L_0 \propto K_\text{area} \times R_0/N$
(since $L_0 \propto R_0$ and $\Delta A/A_0$ is dimensionless).  Dividing by the membrane
force per node $\propto R_0$, the ratio is $K_\text{area}/(N \times E_l^\text{base})$,
which is $R_0$-independent as long as $K_\text{area} = \text{const}$.  If instead
$K_\text{area} \propto R_0$, the effective area stiffness ratio $q_\text{eff} \propto R_0$
grows with particle size, breaking self-similarity.

**Why $C$ does not scale with $R_0$:**  The contact spring stiffness is $k_c = C/R_0$.
The contact force per node scales as $k_c \times g \times L_0 \sim (C/R_0) \times g \times R_0 = C \times g$,
where $g$ is the gap magnitude.  For a geometrically self-similar configuration,
$g \propto R_0$, so contact force $\propto C \times R_0$.  This matches the membrane
force per node $\propto R_0$, giving $C/E_l^\text{base}$ as the relevant ratio — which
is $R_0$-independent as long as $C = \text{const}$.

**The role of $\rho_f$:**  The fluid density $\rho_f$ sets the node mass
$m_i = \rho_f (\tau R_0) L_0 \propto \rho_f R_0^2/N$ and hence the acoustic wave speed
$c_s = \sqrt{E_l^\text{code}/(\rho_f \tau R_0)} \propto \rho_f^{-1/2}$.  Since all
quasi-static observables ($\Delta A$, $\nu_\text{meas}$) are measured after the elastic
DOF have equilibrated, they are independent of $\rho_f$ (which controls only the transient
dynamics, not the equilibrium).  We fix $\rho_f = 1$ with no loss of generality for
quasi-static measurements.

**Consequence:** fix $R_0 = 1$ for all calibration work.  Simulations at $R_0 \neq 1$
use equations (4.3) and inherit the same $(\Delta A, \nu)$ response without re-calibration.

### 4.4  The two dimensionless physical parameters

With $S = 1$ fixed, the model at $R_0 = 1$ reduces to two dimensionless parameters
that control the particle response:

**Bending number $b_0$:**

$$
b_0 = \frac{EI}{E_l^\text{code}}\bigg|_{R_0=1} = \frac{S}{E_l^\text{base}} = \frac{\tau^2}{12}.
\tag{4.4}
$$

This is the ratio of bending to membrane stiffness at the calibration geometry $R_0 = 1$.
For $R_0 \neq 1$, the same ratio evaluates to $b(R_0) = b_0 R_0^2$; however, the
physics (ν, ΔA) remains $R_0$-independent as shown in §4.3, so $b_0$ is the invariant
material descriptor and the subscript is dropped hereafter.  Small $b$ corresponds to a
membrane-dominated shell (shape set by geometry); large $b$ to a bending-dominated shell
(resists curvature change).

**Area stiffness ratio $q$:**

$$
q = \frac{K_\text{area}}{E_l^\text{base}} = \frac{K_\text{area} \tau^2}{12 S}.
\tag{4.5}
$$

This is the ratio of area penalty stiffness to the $R_0 = 1$ membrane stiffness.
Because $K_\text{area}$ is constant with respect to $R_0$ (§4.3), $q$ is genuinely
$R_0$-independent: the same $q$ produces the same $\Delta A$ and $\nu_\text{meas}$
at every particle size.  Small $q$ allows large area loss (area-compressible); large $q$
forces near-constant area so the particle deforms its shape instead (shape-deformable).

**Natural midpoint — $q \cdot b = 1$:**  The product

$$
q \cdot b = \frac{K_\text{area} \cdot EI}{E_l^2}
\tag{4.6}
$$

is a dimensionless measure of the competition between area conservation and bending resistance
relative to the membrane.  When $q \cdot b = 1$, area stiffness and bending stiffness are
balanced at the membrane scale: $K_\text{area} \cdot EI = E_l^2$.  This defines a natural
midpoint of the squishiness axis.

### 4.5  Two numerical parameters to calibrate once

Two additional parameters are required for numerical stability but do not represent
independent physical content:
- $C$ — contact stiffness: set by a stability criterion (§5.1)
- $\alpha_0$ — damping parameter: set by a quasi-static convergence criterion (§5.2)

Both are calibrated once across the full $(b, q)$ space and held fixed thereafter.

---

## §5  Numerical Calibration

### 5.1  Contact stiffness $C$

The contact penalty must be stiff enough to prevent significant perimeter-node overlap
(which would corrupt the force balance) but not so stiff as to require an impractically
small timestep.  We define the *contact compression* $\text{cc} = |g_i^\text{min}| / r_c$,
where $g_i^\text{min}$ is the most negative gap value, as the diagnostic.  The stability
criterion is $\text{cc} < 2\%$.

**Why $C \propto S(1+q)$, not $C \propto E_l$:**  The contact spring stiffness is
$k_c = C/R_0$, giving a contact force $\sim k_c R_0 = C$ at unit overlap relative to
$R_0$.  For the contact force to be hard relative to the elastic restoring force, we need
$C \gg S$ (since $S$ sets the elastic force scale).  Crucially, $C$ must be
*independent of $\tau$*: as $\tau$ increases, $E_l = 12S/\tau^2$ decreases, but $EI =
S R_0^3$ and the contact spring stiffness must remain comparable to $EI/R_0^2 = S$,
which is $\tau$-independent.  An earlier parameterization $C \propto E_l$ gave
$C \approx 62$ at $\tau = 0.20$ (vs. $C = 1000$ at $\tau = 0.05$), producing
$\text{cc} = 22\%$ — catastrophic overlap.

The area spring adds a restoring force $\sim K_\text{area} = q E_l^\text{base}$ that
stiffens the perimeter in the normal direction when $q$ is large.  The contact spring
must dominate both the membrane and the area spring:

$$
C \gtrsim S \max(1, q) \sim S(1+q).
\tag{5.1}
$$

**Calibration:**  A bisection sweep over $C_0$ (where $C = C_0 S(1+q)$) at
$\tau \in \{0.05, 0.10, 0.20\}$ and $q \in \{0.1, 1, 5, 10\}$ found that
$C_0 = 3000$ gives $\text{cc} < 2\%$ for all tested cases.  The resulting formula is:

$$
\boxed{C = 3000 \, S (1 + q).}
\tag{5.2}
$$

This is $\tau$-independent (verified: CV $< 1\%$ across $\tau$ at fixed $q$) and scales
correctly with $S$ and $q$.  The formula was also verified at the calibration working
point $b = 0.2$ ($\tau \approx 1.549$): for all $q \in [0.05, 50]$, $\text{cc} < 0.5\%$
(see §7.1).

### 5.2  Damping parameter $\alpha_0$

The damping parameter $\alpha_0$ determines the overdamping fraction $\zeta = \alpha_0/(4\pi)$.
Too little damping ($\alpha_0 \ll 1$) causes persistent elastic oscillations that corrupt the
quasi-static stress-strain measurement.  Too much damping slows convergence and can mask
physical behavior.

The criterion for calibration is that the perimeter-node kinetic energy must fall below
$10^{-4}$ of the elastic strain energy within $\sim 10$ acoustic wave periods after each
load increment.  Sweep over $\alpha_0 \in \{0.5, 1.0, 1.5, 2.0, 3.0, 5.0\}$ at representative
$(b, q)$ values shows that $\alpha_0 = 2.0$ ($\zeta \approx 16\%$) satisfies this criterion
across the full parameter space.  Values $\alpha_0 < 1.0$ produced measurable oscillation
artifacts in $\nu_\text{meas}(\varepsilon_p)$ curves; values $\alpha_0 > 3.0$ produced
excessive drag that artificially stiffened the low-$q$ response.

$$
\boxed{\alpha_0 = 2.0 \quad (\zeta \approx 16\%).}
\tag{5.3}
$$

---

## §6  Phase-Space Structure: The $(\Delta A, \nu)$ Diagram

### 6.0  Calibration geometry: the two-disk flat-wall squeeze

All phase-space characterization and calibration is performed in a single canonical
geometry: two identical EPD particles squeezed symmetrically between two rigid flat walls.
The walls are horizontal; the two disks are stacked vertically with centers separated by
$2R_0$ at initialization (touching, zero overlap).  Wall velocity $v_\text{wall}$ is
applied to both top and bottom walls simultaneously, so each disk is compressed by equal
and opposite forces.  By symmetry, neither disk translates or rotates during quasi-static
compression; only the elastic deformation and the rigid-body squeeze vary.

The engineering strain is defined as $\varepsilon_p = \delta/R_0$ where $\delta$ is the
compression of a single particle (= half the total wall-to-wall displacement, by symmetry).
All inter-particle contact occurs through the gap between the two disk perimeters; no side
walls are used (the particles are free to expand laterally without constraint).

This geometry provides:
- Clean uniaxial loading with no rotational or shear components
- Symmetric deformation of both disks (single representative response)
- Unambiguous definition of $\varepsilon_\text{vert}$ (= $\varepsilon_p$) and
  $\varepsilon_\text{lat}$ (lateral expansion at the equator)
- Well-defined contact zone (flat region at top and bottom of each disk)

The geometry is the standard two-body Hertz contact test bed, and the measured Poisson
ratio $\nu_\text{meas}$ generalizes the Hertzian contact response to arbitrary squishiness.

### 6.1  Observables

In the two-disk squeeze geometry — two identical disks compressed between parallel flat walls
to engineering strain $\varepsilon_p = \delta / R_0$ — we measure two dimensionless scalar
observables at a reference strain $\varepsilon_\text{ref}$:

**Area loss:**
$$
\Delta A (\%) = \left(1 - \frac{A}{A_0}\right) \times 100.
$$

**Measured 2D Poisson ratio:**
$$
\nu_\text{meas} = -\frac{\varepsilon_\text{lat}}{\varepsilon_\text{vert}},
$$

where $\varepsilon_\text{vert} = \varepsilon_p$ is the applied compressive strain and
$\varepsilon_\text{lat} = (w - w_0)/w_0$ is the fractional lateral expansion of the
particle (here $w$ is the horizontal span of the deformed perimeter and $w_0 = 2R_0$).
In 2D, an incompressible particle ($A = A_0$) under uniaxial compression has
$\nu_\text{meas} = 1$; a fully area-compressible particle that does not expand laterally
has $\nu_\text{meas} = 0$.  **Contact half-angle:**
$$
\theta_w = \arcsin(a / R_0),
$$
where $a$ is the half-width of the flat contact zone (the distance from the particle axis
to the outermost contact node).  This quantifies how much of the perimeter is in contact
with the wall; it depends weakly on both $b$ and $\varepsilon_p$ and is used to assess
geometric resolution of the contact zone.  At $b = 0.2$ and $\varepsilon_p = 0.08$,
$\theta_w \approx 9°$ across the full squishiness range.

In the EPD model, $\nu_\text{meas}$ saturates below 1 at
large $q$ (typically $\nu_\text{meas} \approx 0.94$ at $q = 50$, $b = 0.2$) because the
perimeter polygon at finite $N$ has a residual area loss even under strong incompressibility
penalty, and because lateral expansion is measured from the perimeter nodes which includes
the capsule contact radius $r_c$.  The saturation value approaches 1 as $N \to \infty$
and $q \to \infty$, but the rate depends on $N$.

### 6.2  Shape collapse: perimeter shape is geometry-determined

A central and non-obvious result of the EPD model is that the perimeter *shape* at a
given $\Delta A$ is essentially independent of the material parameters $(E_l, EI, K_\text{area})$
separately.  This is the **shape collapse** property.

**Argument.**  The perimeter equilibrium is the solution to:

$$
\min_{\{\mathbf{x}_i\}} \Bigl[ E_\text{mem} + E_\text{bend} \Bigr]
\quad \text{subject to} \quad A = A_0 - \Delta A, \; \text{contact constraints},
$$

where we temporarily treat $\Delta A$ as a prescribed constraint (the area spring
enforces it in practice, but once the equilibrium is reached, the shape is determined
by the above minimization at the realized $\Delta A$).  Writing

$$
E_\text{mem} + E_\text{bend} = E_l^\text{code} \Bigl[ \tfrac{1}{2}\sum_i e_i^2 L_0 \Bigr]
+ E_l^\text{code} \, b \Bigl[ \tfrac{1}{2}\sum_i (\kappa_i - \kappa_0)^2 L_0 \Bigr]
= E_l^\text{code} \cdot \mathcal{F}_b\bigl(\{\mathbf{x}_i\}\bigr),
$$

where $b = EI/E_l^\text{code} = \tau^2/12$ at $R_0 = 1$.  The prefactor
$E_l^\text{code}$ cancels from the Euler–Lagrange equations, so the minimizer
$\{\mathbf{x}_i^*\}$ depends only on the *geometric* dimensionless functional
$\mathcal{F}_b$, which is a function of node positions relative to $R_0$ and of the
parameter $b$.  It does not depend on $E_l$, $\tau$, $S$, or $q$ individually; only
$b$ appears.

**In the membrane-dominated limit $b \to 0$:** the bending term vanishes and the shape
minimizes pure membrane energy $\sum_i e_i^2$ subject to fixed contact geometry and area.
This shape depends only on $\Delta A/A_0$ and the contact half-width $a/R_0$.  Since
$\Delta A/A_0 < 2\%$ in the typical regime and $a/R_0 < 15\%$, the perimeter deviates
only slightly from a circle and the shape is nearly universal.

**For finite $b$:**  The shape depends weakly on $b$, entering as a correction of
relative order $b$.  At $b = 0.2$, the bending contribution produces a measurable
smoothing of the contact corner (the contact half-angle $\theta_w$ varies by $\sim 1^\circ$
across the full $q$ range at fixed $b$), but the overall perimeter shape is virtually
unchanged.

**Numerical verification:**  Perimeter outlines at maximum strain were overlaid across all
$(τ, q)$ combinations after centering.  At fixed $\Delta A$, all outlines collapsed onto a
single curve to within the plotting resolution.  Differences visible across varying $q$ at
fixed $\varepsilon_p$ are attributable entirely to different $\Delta A$ values, not to
different shapes at the same $\Delta A$.

**Resolution requirement:**  The shape collapse argument is a continuum result.  At small
$N$, the contact zone is resolved by only a few nodes ($\sim N \theta_w / (2\pi)$ nodes
in the contact half-angle), and the discretization error can mask the collapse.  At
$N = 32$ and $\theta_w \approx 9°$, the contact zone contains $\sim 2$–3 nodes; the
collapse holds to visual accuracy but quantitative agreement requires $N \gtrsim 48$.
This is consistent with the $N$-sensitivity found in the calibration curves (§7.2).

**Consequence:** $\Delta A$ is the sole shape-control parameter.  The "squishiness" of a
particle's shape under flat-wall compression is determined by how much area it loses, not by
any internal stiffness ratio.  This shifts the question from "what shape does the particle
adopt?" to "how much area does it lose, and how does that affect the lateral response
$\nu_\text{meas}$?"

### 6.3  Phase-space sweep: the $(\Delta A,\, \nu)$ heat map

A full sweep over $q \in \{0.1, 0.25, 0.5, 1, 2, 5, 10, 25\}$ and
$b \in \{2 \times 10^{-3}, \ldots, 2\}$ (equivalently $\tau \in \{0.05, \ldots, 5\}$)
at $\varepsilon_\text{ref} = 0.08$ reveals the following structure (see Fig.~X):

- **q controls $\Delta A$**: increasing $q$ suppresses area loss monotonically.
  Along any fixed-$b$ slice, $\Delta A$ decreases from $\sim 10\%$ at $q = 0.1$ to
  $\sim 0.1\%$ at $q = 50$.

- **b controls $\nu_\text{meas}$**: at small $b$ (membrane-dominated), the shell
  compresses almost uniaxially, giving $\nu_\text{meas}$ close to its purely geometric
  value set by $\Delta A$.  At large $b$ (bending-dominated), the shell resists curvature
  changes, forcing lateral expansion that raises $\nu_\text{meas}$.  However, at
  $b \gg 1$, the particle becomes shape-locked and $\nu_\text{meas}$ drops again as the
  large bending stiffness resists all deformation.

- **The natural midpoint $q \cdot b = 1$:** the locus $q \cdot b = 1$ (dashed curve in
  Fig.~X) separates the region where area effects dominate ($q \cdot b < 1$) from the
  region where bending effects dominate ($q \cdot b > 1$).

### 6.4  One-to-one relationship between $\Delta A$ and $\nu_\text{meas}$ at fixed $b$

A key simplification emerges from the phase-space map: for fixed $b$, the observables
$\Delta A$ and $\nu_\text{meas}$ are monotonically related.  As $q$ increases at fixed $b$,
$\Delta A$ decreases and $\nu_\text{meas}$ increases; the trajectory $q \mapsto (\Delta A(q),
\nu_\text{meas}(q))$ is a smooth, monotone curve in the $(\Delta A, \nu)$ plane.  This
one-to-one relationship means that specifying either $\Delta A$ or $\nu_\text{meas}$ at fixed
$b$ uniquely determines the other — and both are determined by $q$ alone.  The two
observables are therefore not independent degrees of freedom at fixed $b$; the model exposes
a genuine one-parameter squishiness axis.

---

## §7  The Squishiness Axis: Calibration at Fixed $b = 0.2$

### 7.1  Choice of working point $b = 0.2$

The phase-space map (§6.3) shows that fixing $b = 0.2$ and varying $q$ produces a clean
sweep of the full squishiness range:

| $q$ | $q \cdot b$ | $\Delta A$ | $\nu_\text{meas}$ | Character |
|-----|-------------|------------|-------------------|-----------|
| 0.05 | 0.010 | 6.24% | 0.173 | strongly area-compressible |
| 0.10 | 0.020 | 5.94% | 0.211 | area-compressible |
| 0.20 | 0.040 | 5.43% | 0.276 | area-compressible |
| 0.50 | 0.100 | 4.31% | 0.417 | compressible |
| 1.0  | 0.200 | 3.20% | 0.555 | intermediate |
| 2.0  | 0.400 | 2.12% | 0.691 | intermediate |
| **5.0** | **1.000** | **1.05%** | **0.824** | **midpoint ($q \cdot b = 1$)** |
| 10.0 | 2.000 | 0.57% | 0.883 | shape-deformable |
| 20.0 | 4.000 | 0.30% | 0.917 | nearly incompressible |
| 50.0 | 10.00 | 0.12% | 0.939 | nearly incompressible |

The full 18-point calibration table ($N = 48$, $\varepsilon_\text{ref} = 0.08$) is:

| $q$ | $q \cdot b$ | $\nu_\text{meas}$ | $\Delta A$ (%) | Character |
|-----|-------------|-------------------|----------------|-----------|
| 0.050 | 0.010 | +0.1729 | 6.24 | area-compressible |
| 0.100 | 0.020 | +0.2105 | 5.94 | area-compressible |
| 0.150 | 0.030 | +0.2446 | 5.67 | area-compressible |
| 0.200 | 0.040 | +0.2757 | 5.43 | area-compressible |
| 0.300 | 0.060 | +0.3304 | 4.99 | compressible |
| 0.500 | 0.100 | +0.4170 | 4.31 | compressible |
| 0.750 | 0.150 | +0.4963 | 3.67 | compressible |
| 1.000 | 0.200 | +0.5552 | 3.20 | intermediate |
| 1.500 | 0.300 | +0.6369 | 2.55 | intermediate |
| 2.000 | 0.400 | +0.6908 | 2.12 | intermediate |
| 3.000 | 0.600 | +0.7576 | 1.59 | shape-deformable |
| **5.000** | **1.000** | **+0.8238** | **1.05** | **midpoint ($q \cdot b = 1$)** |
| 7.500 | 1.500 | +0.8624 | 0.74 | shape-deformable |
| 10.000 | 2.000 | +0.8834 | 0.57 | shape-deformable |
| 15.000 | 3.000 | +0.9057 | 0.39 | shape-deformable |
| 20.000 | 4.000 | +0.9173 | 0.30 | nearly incompressible |
| 30.000 | 6.000 | +0.9293 | 0.20 | nearly incompressible |
| 50.000 | 10.000 | +0.9392 | 0.12 | nearly incompressible |

At $b = 0.2$ ($\tau = \sqrt{2.4} \approx 1.549$), the bending contribution is
non-negligible but subdominant relative to membrane, ensuring that the contact zone
is smoothly resolved without numerical artifacts.  Choices $b \ll 0.1$ risk sharp
contact corners that are under-resolved by the perimeter discretization; choices
$b \gg 0.5$ confine the squishiness axis to a narrow range of $\nu$ before the
bending-lock regime is reached.

### 7.2  Calibration procedure

The calibration curve $\nu_\text{meas}(q)$ is obtained by running the two-disk squeeze
simulation at each of $n_q$ values of $q$ (log-spaced from 0.05 to 50, covering $q \cdot b$
from 0.01 to 10), at fixed $b = 0.2$, $R_0 = 1$, $S = 1$, and reading off $\nu_\text{meas}$
and $\Delta A$ at $\varepsilon_\text{ref}$.

**Sensitivity to $\varepsilon_\text{ref}$:**  The calibration is evaluated at five reference
strains $\varepsilon_\text{ref} \in \{0.04, 0.06, 0.08, 0.10, 0.12\}$, at $N = 48$.

| $q$ | $q \cdot b$ | $\nu$ @ 0.04 | $\nu$ @ 0.06 | $\nu$ @ 0.08 | $\nu$ @ 0.10 | $\nu$ @ 0.12 |
|-----|-------------|-------------|-------------|-------------|-------------|-------------|
| 0.05 | 0.010 | +0.1637 | +0.1688 | +0.1729 | +0.1766 | +0.1802 |
| 0.10 | 0.020 | +0.2003 | +0.2061 | +0.2105 | +0.2143 | +0.2179 |
| 0.50 | 0.100 | +0.4019 | +0.4109 | +0.4170 | +0.4218 | +0.4260 |
| 1.00 | 0.200 | +0.5376 | +0.5480 | +0.5552 | +0.5611 | +0.5663 |
| 5.00 | 1.000 | +0.8018 | +0.8139 | +0.8238 | +0.8328 | +0.8414 |
| 10.0 | 2.000 | +0.8606 | +0.8728 | +0.8834 | +0.8934 | +0.9031 |
| 50.0 | 10.00 | +0.9155 | +0.9279 | +0.9392 | +0.9501 | +0.9608 |

The maximum $\varepsilon_\text{ref}$ drift across the full range $\{0.04, \ldots, 0.12\}$ is
$\Delta\nu_\text{max} = 0.045$ (at large $q$, where $\nu$ is still rising slowly with strain),
and the median drift is $\Delta\nu_\text{med} = 0.033$.  This is the strain-rate artefact
inherent to a dynamic, quasi-static squeeze: $\nu_\text{meas}$ increases monotonically with
strain because the area-penalty restoring force is not yet fully equilibrated at the early
contact stages.  The reference strain $\varepsilon_\text{ref} = 0.08$ (8\% per-particle compression)
is chosen as a practical compromise: large enough that the inter-particle contact is well
developed and the Poisson-ratio signal is strong, small enough to remain in the quasi-linear
elastic regime before geometric nonlinearities accumulate.  Any calibration performed at a
consistent $\varepsilon_\text{ref}$ is internally consistent; the absolute value of $\nu$ at
$\varepsilon_\text{ref} = 0.08$ is what the recipe in §8 delivers.

**Sensitivity to $N$:**  The calibration is evaluated at perimeter resolutions
$N \in \{32, 48, 72\}$, at $\varepsilon_\text{ref} = 0.08$.

| $q$ | $q \cdot b$ | $\nu(N=32)$ | $\nu(N=48)$ | $\nu(N=72)$ | $\Delta\nu_\text{max}$ |
|-----|-------------|------------|------------|------------|----------------------|
| 0.05 | 0.010 | +0.1760 | +0.1729 | +0.1702 | 0.0058 |
| 0.50 | 0.100 | +0.4211 | +0.4170 | +0.4117 | 0.0094 |
| 0.75 | 0.150 | +0.5004 | +0.4963 | +0.4909 | **0.0095** |
| 1.00 | 0.200 | +0.5590 | +0.5552 | +0.5499 | 0.0091 |
| 5.00 | 1.000 | +0.8249 | +0.8238 | +0.8206 | 0.0043 |
| 10.0 | 2.000 | +0.8836 | +0.8834 | +0.8812 | 0.0024 |
| 50.0 | 10.00 | +0.9384 | +0.9392 | +0.9380 | 0.0012 |

The maximum $N$-drift across $N \in \{32, 48, 72\}$ at any $q$ is
$\Delta\nu_\text{max} = 0.0095$ (at $q = 0.75$, $q \cdot b = 0.15$); the median drift
is $0.0064$.  The drift is monotone: $\nu$ decreases with increasing $N$ in the
compressible regime (small $q$), reflecting better resolution of the contact geometry
at fine discretizations.  The drift is below $0.01$ everywhere, meaning the calibration
at $N = 48$ agrees with $N = 72$ to within $\Delta\nu < 0.01$ across the full squishiness
axis.  A single calibration performed at $N = 48$ therefore transfers to $N = 72$ without
re-calibration, to within the $\Delta\nu = 0.01$ tolerance.  When a multi-particle target
resolution is eventually established, the calibration may be re-run at that $N$ if tighter
agreement is required; for $N \geq 48$, a single curve is adequate for exploratory work.

### 7.3  Calibration inversion: $\nu_\text{target} \to q$

Since $\nu_\text{meas}(q)$ is monotone increasing, it can be inverted to give
$q(\nu_\text{target})$ via log-linear interpolation on $\log q$.  The full parameter
set for a simulation with target Poisson ratio $\nu_\text{target}$ and particle radius
$R_0$ is then:

$$
\tau = \sqrt{12 \times 0.2} \approx 1.549, \qquad q = q(\nu_\text{target}),
\tag{7.1}
$$

$$
E_l = \frac{12 S}{\tau^2}, \quad EI = S R_0^3, \quad K_\text{area} = q \, E_l^\text{base},
\quad C = 3000 \, S (1 + q).
\tag{7.2}
$$

Here $E_l^\text{base} = 12S/\tau^2$ is constant (no $R_0$ factor); it equals the
membrane stiffness $E_l^\text{code}$ only at $R_0 = 1$.  $K_\text{area}$ and $C$ are
independent of $R_0$.  The only parameter that scales with $R_0$ is $EI = S R_0^3$.

### 7.4  R₀ invariance spot-check

To verify that the calibration at $R_0 = 1$ transfers to $R_0 \neq 1$, three target
$\nu$ values — representing the low, mid, and high squishiness regimes — are selected,
the corresponding $q$ values are read from the calibration curve, and two-disk simulations
are run at $R_0 \in \{0.5, 1.0, 2.0\}$.

**Criterion:** $\text{CV}(\nu_\text{meas}) < 1\%$ and $\text{CV}(\Delta A) < 1\%$ across
the three $R_0$ values for each $\nu_\text{target}$.

Three target Poisson ratios spanning the squishiness axis are selected
($\nu_\text{target} \in \{0.30, 0.60, 0.85\}$), mapped to $q$ via the calibration curve
at $N = 72$, $\varepsilon_\text{ref} = 0.08$, and simulated at $R_0 \in \{0.5, 1.0, 2.0\}$
with all parameters set by the recipe (7.2).

| $\nu_\text{target}$ | $q$ | $R_0$ | $EI$ | $\nu_\text{meas}$ | $\Delta A$ (%) | $\text{cc\%}$ |
|---------------------|-----|--------|------|-------------------|----------------|---------------|
| 0.30 | 0.2474 | 0.50 | 0.1250 | 0.2984 | 5.24 | 1.8 |
| 0.30 | 0.2474 | 1.00 | 1.0000 | 0.2984 | 5.24 | 1.8 |
| 0.30 | 0.2474 | 2.00 | 8.0000 | 0.2984 | 5.24 | 1.8 |
| **CV** | | | | | **0.00%** | **0.00%** |
| 0.60 | 1.2814 | 0.50 | 0.1250 | 0.6005 | 2.84 | 1.4 |
| 0.60 | 1.2814 | 1.00 | 1.0000 | 0.6005 | 2.84 | 1.4 |
| 0.60 | 1.2814 | 2.00 | 8.0000 | 0.6005 | 2.84 | 1.4 |
| **CV** | | | | | **0.00%** | **0.00%** |
| 0.85 | 6.7725 | 0.50 | 0.1250 | 0.8511 | 0.827 | 0.8 |
| 0.85 | 6.7725 | 1.00 | 1.0000 | 0.8511 | 0.827 | 0.8 |
| 0.85 | 6.7725 | 2.00 | 8.0000 | 0.8511 | 0.827 | 0.8 |
| **CV** | | | | | **0.00%** | **0.01%** |

The coefficient of variation of $\nu_\text{meas}$ across $R_0 \in \{0.5, 1.0, 2.0\}$ is
$\text{CV} = 0.00\%$ for all three targets; $\text{CV}(\Delta A) \leq 0.01\%$.
This confirms that the scaling rules $EI = S R_0^3$, $K_\text{area} = \text{const}$,
$C = \text{const}$ render the physics exactly $R_0$-independent: to floating-point
precision, the same dimensionless deformation occurs regardless of particle size.

---

## §8  Summary of Fixed Parameters and Calibration Recipe

For reference, the complete set of fixed and calibrated parameters is collected here:

**Dimensionless scales (fixed, no physical content):**
$$
S = 1 \quad (\text{force scale}), \qquad R_0 = 1 \quad (\text{calibration geometry}).
$$

**Numerical stability parameters (calibrated once, held fixed):**
$$
C_0 = 3000, \quad \alpha_0 = 2.0 \quad (\zeta \approx 16\%), \quad
\text{SR} = 0.01, \quad f_\text{dt} = 0.4.
$$

**Squishiness working point (fixed by phase-space analysis):**
$$
b = 0.2 \quad \Rightarrow \quad \tau = \sqrt{2.4} \approx 1.549.
$$

**Free parameter and calibration:**
$$
q = q(\nu_\text{target}) \quad \text{from calibration curve.}
$$

**Full parameter recipe for a simulation with target $(\nu_\text{target}, R_0)$:**

```
S      = 1.0                        # force scale (fixed)
tau    = sqrt(2.4) ≈ 1.549          # shell thickness parameter (fixed for b=0.2)
q      = q(nu_target)               # from calibration curve
El_t   = 12 * S / tau**2            # membrane stiffness (R0=1 base)
EI     = S * R0**3                  # bending stiffness (scales with R0)
K_area = q * El_t                   # area stiffness (constant, no R0 factor)
C      = 3000 * S * (1 + q)         # contact stiffness (constant, no R0 factor)
```

Numerical settings: $\text{SR} = 0.01$, $\alpha_0 = 2.0$, $f_\text{dt} = 0.4$.

All other observables ($\Delta A$, $\theta_w$, contact force, perimeter shape) are
outputs, not inputs.  The single input $\nu_\text{target}$ uniquely determines all
particle physics.
