# EPD Simulator — Python API Reference

2D Elastic Particle DEM (EPD) and Emulsion Droplet simulator.
Particles are closed elastic membranes (capsule shells) with a tunable squishiness knob ν.
The simulator is built on TensorFlow and runs efficiently on GPU.

**Design target:** GPU (Colab T4/A100, RTX 3080+). CPU runs are correct but slow.

---

## Contents

1. [Installation](#1-installation)
2. [Quick Start](#2-quick-start)
3. [ParticleSpec](#3-particlespec)
4. [System](#4-system)
5. [Boundary Objects](#5-boundary-objects)
6. [MotionSpec](#6-motionspec)
7. [State Tensors and Snapshots](#7-state-tensors-and-snapshots)
8. [Callbacks](#8-callbacks)
9. [TF Loop Compliance — Critical](#9-tf-loop-compliance--critical)
10. [Recording and Movies](#10-recording-and-movies)
11. [Live Parameter Updates](#11-live-parameter-updates)
12. [Checkpointing](#12-checkpointing)
13. [Background Flow](#13-background-flow)
14. [Known Limitations and Areas for Development](#14-known-limitations-and-areas-for-development)

---

## 1. Installation

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib imageio tensorflow
pip install -e .
```

Verify:

```python
import tensorflow as tf
from src.epd.particles import ParticleSpec
from src.epd.objects   import Wall
from src.epd.system    import System
from src.epd.motion    import MotionSpec
print('OK', tf.__version__)
```

**Do not install:** polyfempy, torch, meshio — these are for an axed FEM/NN path.

---

## 2. Quick Start

```python
import tensorflow as tf
import src.simulation.tf_sim as tf_sim_mod
tf_sim_mod.set_dtype(tf.float64)   # must be called once before any System

from src.epd.particles import ParticleSpec
from src.epd.objects   import Wall, Box
from src.epd.system    import System

# Build system
sys = System(Lx=20.0, Ly=20.0, periodic_x=True, periodic_y=True)

# Add particles
sys.add_particles(ParticleSpec(count=40, type='elastic', nu=0.5, N_nodes=60))

# Initialize: RSA seed → adaptive swell to phi_target
sys.initialize(phi_target=0.80, seed=42, verbose=True)

# Run
def my_callback(s):
    snap = s.snapshot()
    return {'t': float(snap['t'])}

sys.run(10000, sample_every=200, callback=my_callback, verbose=True)

# Render
sys.make_movie('output/run.gif', fps=10)
```

---

## 3. ParticleSpec

```python
from src.epd.particles import ParticleSpec
```

Specifies a population of particles. Multiple specs can be added to one System
(different types, sizes, or parameters).

### 3.1 Constructor

```python
ParticleSpec(
    count,                # int — number of particles
    type   = 'elastic',  # str — 'elastic' | 'emulsion' | 'rigid'

    # --- elastic high-level (preferred) ---
    nu     = None,        # float — effective Poisson ratio 0.18–0.94; derives q
    # --- emulsion high-level ---
    gamma  = None,        # float — line tension (surface energy/length); default 1.0
    kappa  = None,        # float — area compressibility ratio γ/(R0·K_area); default 0.02
    Oh     = None,        # float — Ohnesorge number (sets viscous drag; None = no drag)
    # --- low-level overrides ---
    q      = None,        # float — K_area/El_t; overrides nu if given
    tau_b  = 0.2,         # float — bending working point b (default 0.2)
    alpha_damp = None,    # float — damping coefficient; auto if None
    C      = None,        # float — contact hardness; auto if None

    # --- size distribution ---
    N_nodes   = 60,       # int   — perimeter nodes per particle; increase for smoother shapes
    R0_mean   = 1.0,      # float — mean radius; normalised to 1.0 exactly
    poly_dist = None,     # None | float | dict — see §3.3

    # --- motion / shape ---
    motion       = None,  # MotionSpec — drives CM; None = free particle
    frozen_shape = False, # bool — freeze elastic DOFs (rigid body only)
)
```

### 3.2 Parameter derivation

The user sets `nu` (elastic) or `kappa`/`gamma` (emulsion). Everything else is derived:

| Particle type | Key input | Derives |
|---------------|-----------|---------|
| `elastic` | `nu` (0.18–0.94) | `q` via calibration table, `El_t`, `K_area`, `C`, `alpha_damp` |
| `emulsion` | `gamma`, `kappa` | `K_area = gamma/(kappa*R0)`, `C = 500*gamma/R0`, `alpha_damp = 5.0` |
| `rigid`    | — | large `C`, same edge/area springs, shape frozen |

**Calibration table:** `results/calibration_sweep/calibration_data.json` (N=32, eps_ref=0.08).
Must be present for `type='elastic'` with `nu` set.

Derived values are available as:

```python
spec = ParticleSpec(count=10, type='elastic', nu=0.5)
print(spec.derived)
# {'TAU': ..., 'El_t': ..., 'K_area': ..., 'C': ..., 'alpha': ..., 'q': ..., ...}
```

### 3.3 Size distributions (`poly_dist`)

```python
poly_dist = None                              # monodisperse (default)
poly_dist = 0.05                              # Gaussian sigma=0.05 (shorthand)
poly_dist = {'type': 'gaussian', 'sigma': 0.05}
poly_dist = {'type': 'bimodal',  'ratio': 0.5, 'delta': 0.1}  # ratio=large fraction, delta=half-gap
poly_dist = {'type': 'explicit', 'values': [...]}              # len must == count
```

`sample_R0()` returns `(count,)` float64, mean=1.0 exactly.

### 3.4 Ohnesorge number / drag

```python
spec = ParticleSpec(count=10, type='emulsion', gamma=1.0, kappa=0.02, Oh=0.25)
# Terminal velocity under gravity g:
v_t = spec.terminal_velocity(g=0.05)
# Reverse: set Oh to achieve desired terminal velocity:
spec.set_terminal_velocity(v_t=0.1, g=0.05)
```

Drag is off by default (`Oh=None`). Setting `Oh` enables viscous drag per node.

---

## 4. System

```python
from src.epd.system import System
```

### 4.1 Constructor

```python
System(
    Lx, Ly,              # float — simulation box dimensions
    periodic_x = True,   # bool  — periodic boundary in x
    periodic_y = True,   # bool  — periodic boundary in y
    dt_factor  = 0.4,    # float — multiplier on auto-stable dt
    alpha_damp = None,   # float — global damping override; None = auto per-particle
    E_candidates = 128,  # int   — candidacy table width (columns); see §14
    skin       = 1.0,    # float — candidacy skin factor (× mean R0)
    g          = 0.0,    # float — gravitational acceleration (downward)
)
```

### 4.2 Building a system (registration)

```python
sys = System(Lx=20, Ly=20)
sys.add_particles(spec1)          # returns self — chainable
sys.add_particles(spec2)
sys.add_object(wall)              # returns self — chainable
```

Multiple `add_particles` calls are allowed. Particle indices are assigned in
registration order: spec1 first (0..count1-1), then spec2 (count1..count1+count2-1).

### 4.3 Initialization

```python
sys.initialize(
    phi_target   = 0.80,   # float — target packing fraction after swell
    seed         = 42,     # int   — RNG seed for RSA placement and size sampling
    verbose      = True,   # bool
    relax_only   = False,  # bool  — skip swell; settle in-place (use for fixed containers)
    n_relax_init = 200,    # int   — settle steps when relax_only=True
)
```

`initialize()` performs:
1. RSA placement (random sequential addition) at low phi
2. Adaptive box compression to `phi_target`
3. Builds TF state and params dicts
4. Builds CandidacyManager
5. Calls `clear_recording()`

**After `initialize()`, `sys.Lx` and `sys.Ly` may differ from the constructor values**
if auto-expansion was needed for RSA. Always read `sys.Lx`, `sys.Ly` after init.

`initialize_from_particles(particles, ...)` is also available for custom placement.

### 4.4 Running

```python
sys.run(
    N,                       # int  — total steps to run
    sample_every  = None,    # int  — output cadence AND chunk size; default = N
    callback      = None,    # callable(sys) → dict or None
    record_initial = True,   # bool — record snapshot before first step
    cand_check_interval = 10,# int  — candidacy rebuild poll interval
    verbose       = True,    # bool
)
```

`run()` returns `self` (fluent). Frames and callback data **accumulate** across multiple
`run()` calls. Call `clear_recording()` to reset without touching simulation state.

**Performance note:** `sample_every = max(1, N // 2000)` gives ~2000 output points
with minimal Python-side overhead. Per-step output (`sample_every=1`) is slow on GPU.

### 4.5 Key public attributes (post-init)

| Attribute | Type | Description |
|-----------|------|-------------|
| `sys.Lx`, `sys.Ly` | float | Box dimensions (may change after swell) |
| `sys.t` | float | Current simulation time |
| `sys.frames` | list[dict] | Recorded snapshots (appends on each `run()`) |
| `sys.callback_data` | list[dict] | Callback return dicts |
| `sys.diag` | list[dict] | Per-chunk plumbing diagnostics |
| `sys._dt` | float | Integration timestep (auto-stable) |
| `sys._state` | dict[str, tf.Tensor] | Live TF state (see §7) |
| `sys._params` | dict[str, tf.Tensor] | Live TF force parameters |
| `sys._particles` | list | CapsuleParticle / EmulsionParticle objects |

### 4.6 Utility methods

```python
sys.freeze(zero_velocity=True)     # zero all velocities (pre-production run)
sys.clear_recording()              # wipe frames/callback_data/diag, keep state
sys.snapshot()                     # returns current state as numpy dict (see §7)
sys.render(ax=None, output_path=None)  # draw current frame
sys.make_movie(path, fps=10, ...)  # render GIF from sys.frames (non-destructive)
```

---

## 5. Boundary Objects

```python
from src.epd.objects import (
    Wall, OscillatingWall, ArcWall,
    Box, Channel, CouetteCell, CircleObstacle,
    RegularPolygon, SquareObstacle, HopperRegion,
    CustomObject,
)
```

All objects are registered via `sys.add_object(obj)` before `initialize()`.

### 5.1 Wall (single line segment)

```python
Wall(
    p0,             # (2,) — start endpoint
    p1,             # (2,) — end endpoint
    normal = None,  # (2,) — inward normal (toward particles); auto if None
)
wall.set_motion(motion_spec)    # attach motion
wall.set_k_pen(k)               # penalty multiplier (default 1.0)
wall.set_r_c_wall(r)            # wall capsule radius (default 0)
wall.set_render(color, linewidth, alpha, fill, zorder)
```

### 5.2 OscillatingWall

```python
OscillatingWall(
    y0,            # float — rest y position
    half_w,        # float — half-width (spans [-half_w, half_w] in x)
    A,             # float — oscillation amplitude
    omega,         # float — angular frequency (rad/s)
    is_top = True, # bool  — True: normal=[0,-1]; False: normal=[0,+1]
    r_c_wall = 0.0,
)
# y(t) = y0 + sign * A * sin(omega * t)
# OscillatingWall is TF-native: prim_data is static, motion encoded in parameters.
# No prim_data rebuild needed at runtime — this is the most GPU-efficient moving wall.
```

### 5.3 ArcWall (circular boundary)

```python
ArcWall(
    center,          # (2,) — arc center
    radius,          # float
    convex = True,   # True = convex obstacle (pushes outward); False = concave container (pushes inward)
    angle_range = None,  # (theta_min, theta_max) radians or None (full circle)
)
```

### 5.4 Pre-built composite objects

All accept `set_motion()`, `set_k_pen()`, `set_render()`.

```python
Box(width, height, x0=0, y0=0, theta=0, exclusion='exterior')
    # 4 walls; particles inside; typical closed container

Channel(width, height, x0=0, y0=0, exclusion='exterior')
    # top + bottom walls only; pair with periodic_x=True

CouetteCell(inner_radius, outer_radius, x0=0, y0=0)
    # annular cell; particles between the two arcs

CircleObstacle(radius, x0=0, y0=0)
    # solid disk obstacle; particles excluded from interior

RegularPolygon(sides, radius, x0=0, y0=0, theta=0, exclusion='interior')
    # N-sided polygon; exclusion='interior' = obstacle, 'exterior' = container

SquareObstacle(side, x0=0, y0=0, theta=0)
    # square obstacle with Polygon primitive (shared group ID, avoids double-counting)

HopperRegion(x_c, W_out, W_res, theta_deg=30, h_res=42, y_bot=0, corner_radius=0.4)
    # funnel + reservoir; open top and bottom; particles seed inside polygon

CustomObject(x0=0, y0=0, theta=0, exclusion=None)
    # assemble walls manually: obj.add_primitive(Wall(...))
```

### 5.5 Exclusion semantics

`exclusion='exterior'` — accessible region is **inside** the polygon (container).
`exclusion='interior'` — accessible region is **outside** the polygon (obstacle).

The seeder uses `obj.contains(pt)` to reject placements, and `obj.exclusion_area()`
to correct the phi denominator.

---

## 6. MotionSpec

```python
from src.epd.motion import MotionSpec
```

Three modes:

```python
# Mode 1 — parametric DC + AC (most common; TF-native for OscillatingWall)
MotionSpec(vy=-0.1)                              # constant downward velocity
MotionSpec(vx_dc=0.5, vx_ac=1.0, freq_x=2.0)   # DC + oscillating x
MotionSpec(omega_dc=0.3, r_ref=(5.0, 5.0))      # spin about r_ref

# Mode 2 — callable (TF-native when f(t) uses TF ops)
MotionSpec(vx=lambda t: tf.sin(omega * t), vy=0.0)

# Mode 3 — pre-sampled (Python fallback; slower)
MotionSpec.from_samples(vx_fn=lambda t: np.cos(t), dt=0.001, duration=10.0)
```

**Attach to object or particle:**
```python
wall.set_motion(MotionSpec(vy=-0.1))
spec.set_motion(MotionSpec(vx=1.0))   # drives CM of all particles in spec
```

**Displacement integral** (analytic for Mode 1, numerical otherwise):
```python
dx, dy, dtheta = ms.displacement(t)   # used by Wall.resolved(t)
```

---

## 7. State Tensors and Snapshots

### 7.1 State dict (`sys._state`)

All tensors are `tf.float64` unless noted.

| Key | Shape | Description |
|-----|-------|-------------|
| `x_cm`  | `(P, 2)` | CM positions (x, y) |
| `v_cm`  | `(P, 2)` | CM velocities |
| `theta` | `(P,)` | Rotation angles (rad) |
| `omega` | `(P,)` | Angular velocities (rad/s) |
| `u`     | `(P, N, 2)` | Node displacements in body frame |
| `u_dot` | `(P, N, 2)` | Node velocities in body frame |
| `t`     | scalar | Simulation time |
| `step`  | scalar `int64` | Step counter |

P = total number of particles, N = nodes per particle (same for all; set by `N_nodes`).

**Caution:** `x_all = x_cm[:, None, :] + R(theta) @ u` gives world-frame node positions.
`snapshot()` computes and returns `x_all` as a convenience.

### 7.2 Snapshot dict

`sys.snapshot()` returns a Python dict with numpy arrays:

| Key | Shape | Notes |
|-----|-------|-------|
| `x_cm`  | `(P, 2)` float64 | CM positions |
| `v_cm`  | `(P, 2)` float64 | CM velocities |
| `theta` | `(P,)` float64 | Orientations |
| `omega` | `(P,)` float64 | Angular velocities |
| `x_all` | `(P, N, 2)` float64 | World-frame node positions |
| `t`     | scalar float64 | Current time |
| `step`  | scalar int64 | Current step |

Snapshots are stored in `sys.frames` (list of dicts, one per recorded step).

### 7.3 Params dict (`sys._params`)

Key per-particle physics tensors (all `(P,)` float64 unless noted):

| Key | Shape | Description |
|-----|-------|-------------|
| `El_t_per_p`      | `(P,)` | Edge spring stiffness |
| `K_area_per_p`    | `(P,)` | Area spring stiffness |
| `C_per_p`         | `(P,)` | Contact penalty |
| `alpha_damp_per_p`| `(P,)` | Damping coefficient |
| `xi_drag_per_p`   | `(P,)` | Drag per unit arc length |
| `r_c_per_p`       | `(P,)` | Node capsule radius |
| `shape_frozen`    | `(P,)` | 1.0 = frozen shape, 0.0 = free |
| `candidates`      | `(P, E)` int32 | Candidacy table (E = E_candidates) |

**Do not modify `_state` or `_params` directly during a `run()` call.** Safe to patch
between `run()` calls via `tf.constant` assignments (see §9).

---

## 8. Callbacks

```python
def my_callback(sys) -> dict:
    snap = sys.snapshot()          # always safe to call
    x_cm = snap['x_cm']           # (P, 2) numpy
    return {
        't'          : float(snap['t']),
        'mean_y'     : float(x_cm[:, 1].mean()),
        'text'       : f"t={snap['t']:.2f}",   # shown in make_movie overlay
        'wall_strain': ...,        # if present, shown in make_movie time-series panel
    }

sys.run(N, sample_every=200, callback=my_callback)
```

**The callback is called once per chunk** (every `sample_every` steps), not every step.
Callback return dicts are stored in `sys.callback_data`.

**Callbacks are read-only with respect to the TF graph.** You can read state via
`sys.snapshot()` but must not call `sys.run()` or modify `sys._state` inside the callback.
For state modification between chunks, run in a loop (see §9.2).

**Special callback return keys:**
- `'text'` — string overlay on movie frames
- `'wall_strain'` — scalar; plotted in movie time-series panel
- `'eps1'`, `'eps2'` — scalars; averaged and plotted in movie panel

---

## 9. TF Loop Compliance — Critical

The core integrator runs as a compiled TF `while_loop`. This means:

**You cannot call Python-side logic inside the TF graph.** The step function
(`run_fast`) is traced once and reused. Any Python code you want to run must go
in the callback (which executes between chunks, outside the graph).

### 9.1 Use the built-in `run()` method

The recommended pattern for almost all use cases:

```python
sys.run(N, sample_every=S, callback=cb)
```

This is correct, TF-loop-compliant, and efficient. The graph is traced once on
the first `run_fast` call and reused thereafter.

**Avoid** writing a manual Python step loop:

```python
# WRONG — defeats TF compilation, extremely slow on GPU
for step in range(N):
    sys.run_fast(1)
    do_something(sys._state)
```

### 9.2 State patching between `run()` calls

To modify state between batches (e.g., re-pin particle positions in a T1-event
benchmark), run in chunks and patch between calls:

```python
import tensorflow as tf
import numpy as np

for chunk in range(N_chunks):
    sys.run(CHUNK_SIZE, sample_every=SAMPLE, record_initial=(chunk == 0))

    # Safe to patch here — between run() calls, outside TF graph
    x_cm_np = sys._state['x_cm'].numpy().copy()
    x_all_np = sys._state['x_all'].numpy().copy() if 'x_all' in sys._state else None
    v_cm_np  = sys._state['v_cm'].numpy().copy()

    x_cm_np[0]  = target_cm          # re-pin particle 0
    v_cm_np[0]  = np.zeros(2)

    sys._state['x_cm'] = tf.constant(x_cm_np, dtype=tf.float64)
    sys._state['v_cm'] = tf.constant(v_cm_np, dtype=tf.float64)
```

**State patching invalidates the compiled runner** (TF retraces on next `run_fast`
call). Patch as rarely as possible to minimise retrace cost.

### 9.3 Moving walls

For moving walls in TF-loop-compliant code:
- `OscillatingWall` — fully TF-native; no prim_data rebuild needed; preferred.
- `Wall.set_motion(MotionSpec(vy=-c))` — parametric DC motion; TF-native.
- Callable `MotionSpec` with TF ops (e.g., `lambda t: tf.sin(omega*t)`) — TF-native.
- Python-callable `MotionSpec` (numpy functions) — evaluated outside TF loop; requires
  prim_data rebuild every chunk; avoid for performance-critical work.

---

## 10. Recording and Movies

### 10.1 Frame accumulation

`sys.frames` is a Python list of snapshot dicts. It **accumulates** across `run()` calls.
`make_movie` reads from `sys.frames` and does **not** clear it.

```python
sys.run(N1, ...)   # frames: [0..N1//S]
sys.run(N2, ...)   # frames: [0..N1//S, N1//S+1..(N1+N2)//S]

sys.make_movie('full.gif', fps=10)   # renders all accumulated frames

sys.clear_recording()   # wipe frames/callback_data/diag; keep state
sys.run(N3, ...)        # fresh recording from current state
```

### 10.2 make_movie

```python
sys.make_movie(
    output_path,         # str or Path — .gif output
    fps       = 10,      # int — frames per second
    n_arc     = 6,       # int — arc interpolation points per node gap
    xlim      = None,    # (xmin, xmax) or None (auto)
    ylim      = None,    # (ymin, ymax) or None (auto)
    title     = None,    # str — title overlay or None
)
```

Output: animated GIF. Requires `imageio` (`pip install imageio`).
If `callback_data` contains `'wall_strain'` or `'eps1'`/`'eps2'` keys, a time-series
panel is added below the particle view.

---

## 11. Live Parameter Updates

Parameters can be updated after `initialize()`. Changes take effect on the next step.

### 11.1 Via ParticleSpec setters (preferred)

```python
spec.set_Oh(0.5)             # update Ohnesorge; propagates xi_drag to TF params
spec.set_kappa(0.05)         # emulsion: update kappa; propagates K_area
spec.set_nu(0.7)             # elastic: update nu; propagates q, El_t, K_area
spec.set_xi(0.3)             # direct drag coefficient
spec.set_terminal_velocity(v_t=0.1, g=0.05)
```

These push changes to the live TF `_params` immediately via `_push_to_system()`.

### 11.2 Via `sys.set_param()`

```python
sys.set_param('Oh',    0.5,  particles='all')
sys.set_param('kappa', 0.05, particles=[0, 1, 2])
sys.set_param('nu',    0.7,  particles=range(10))
sys.set_param('U_max', 2.0)  # update background flow amplitude
```

Supported names: `'Oh'`, `'kappa'`, `'nu'`, `'xi'`, `'K_area'`, `'El_t'`, `'U_max'`.

### 11.3 Background flow

```python
sys.U_background = None                                    # quiescent (default)
sys.U_background = ('shear',       {'rate': 0.1})          # simple shear dv_x/dy = rate
sys.U_background = ('constant',    {'U': (1.0, 0.0)})      # uniform flow
sys.U_background = ('parabolic',   {'U_max': 2.0, 'H': 10.0})
sys.U_background = ('extensional', {'rate': 0.05})
sys.U_background = ('poiseuille_v',{'U_max': 2.0, 'H': 10.0, 'x_c': 0.0})  # vertical Poiseuille
```

Can be set before or after `initialize()`; takes effect immediately post-init.

---

## 12. Checkpointing

### 12.1 Full checkpoint (recommended)

```python
sys.save('results/run1/checkpoint/')     # saves to directory
sys2 = System.from_file('results/run1/checkpoint/')  # restore
```

### 12.2 Lightweight state save/restore

```python
path = sys.save_state('results/run1/state.npz')   # saves state tensors + box dims
sys.restore_state('results/run1/state.npz')        # overwrites state; must be init'd
```

`restore_state` rescales registered objects to match the checkpoint box and rebuilds
the CandidacyManager. Use to skip expensive swell on repeated runs.

---

## 13. Background Flow

See §11.3 for the setter API. Flow presets and their parameter keys:

| Preset | Parameters | Description |
|--------|-----------|-------------|
| `'zero'` | — | Quiescent (default) |
| `'constant'` | `U: (Ux, Uy)` | Uniform flow |
| `'shear'` | `rate: γ̇` | Simple shear: u_x = γ̇ · y |
| `'parabolic'` | `U_max, H` | Parabolic channel flow |
| `'extensional'` | `rate: ε̇` | Extensional: u_x = ε̇·x, u_y = -ε̇·y |
| `'poiseuille_v'` | `U_max, H, x_c` | Vertical Poiseuille (flow in y) |

Background flow adds a drag-like coupling: nodes are pushed toward the local flow
velocity. Requires `Oh > 0` (drag enabled) to have physical effect.

---

## 14. Known Limitations and Areas for Development

**Candidacy manager (C++ extension, CPU-side):**
- Currently builds an `(P, E_candidates)` table where `E_candidates=128` by default.
  This is intentionally oversize for correctness. Shrinking it (to ~20) would reduce
  GPU memory bandwidth cost by ~6x but requires tuning for the system geometry.
- Update frequency is set by `cand_check_interval` (default 10 steps). A skin-distance
  threshold (rebuild only when a particle moves more than `skin/2`) is not yet
  implemented. Adding it would remove CPU-GPU stall for fast GPU steps.
- Single-threaded. Parallel builds across particles are straightforward but not
  yet implemented.

**All particles must have the same N_nodes:**
The TF state tensor `u` has shape `(P, N, 2)` with a fixed N. Mixed `N_nodes`
populations are not supported in a single System.

**Polydisperse systems with very different R0:**
The timestep `dt` is set from the smallest R0 (most constraining). Wide size
distributions force very small dt; consider using `dt_factor < 0.4` cautiously.

**N_nodes and accuracy tradeoff:**
The ν calibration table was built at N=32. Using higher `N_nodes` (60, 120, 240)
gives smoother shapes and higher accuracy, but the calibration interpolation is
based on N=32 measurements. For quantitative work with `N_nodes != 32`, a
re-calibration sweep is recommended.

**`make_movie` for large frame counts:**
All frames are held in memory (`sys.frames`). For very long runs, memory can become
an issue. `clear_recording()` between batches and calling `make_movie` before clearing
is the current workaround.

**No multi-experiment batching:**
Each `System` instance holds one experiment. Running M experiments in parallel requires
M separate processes. A vectorized `(E, P, N, 2)` architecture (Isaac Gym style) is
a future direction that would give near-linear scaling up to GPU memory limits.
See `FUTURE_DIRECTIONS.md` for discussion.

**Driven particles (`set_driven_particles`):**
Low-level API that requires `make_traj()` tuples from `tf_sim.py`. Not yet wrapped
in a clean high-level interface. Use `ParticleSpec(motion=MotionSpec(...))` instead
where possible.
