# EPD Simulator — Python API Reference

2D Elastic Particle DEM (EPD) and Emulsion Droplet simulator.
Particles are closed elastic membranes (capsule shells) with a tunable squishiness knob ν.
The simulator is built on TensorFlow and runs efficiently on GPU.

**Design target:** GPU (Colab T4/A100, RTX 3080+). CPU runs are correct but slow.

**Methods writeup:** see `papers/summary_of_methods/main.pdf` for the full derivation
of the elastic and emulsion force models, calibration of the ν knob, the closing-rate
stability watchdog, and the time-step sweep results.

**Worked examples:** `notebooks/getting_started.ipynb` walks through the canonical
patterns (squeeze tests, hopper flow, Couette shear, parameter sweeps).
`notebooks/reproduce_paper.ipynb` regenerates every paper figure end-to-end.

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
14. [Closing-rate Stability Watchdog](#14-closing-rate-stability-watchdog)
15. [Parallel Sweeps — Running N Experiments Concurrently](#15-parallel-sweeps--running-n-experiments-concurrently)
16. [Known Limitations and Areas for Development](#16-known-limitations-and-areas-for-development)

---

## 1. Installation

### 1a. Google Colab (recommended for first-time users)

The setup cell at the top of either notebook clones this repo, installs Python
deps, and builds the C++ candidacy-manager extension automatically. Open
`notebooks/getting_started.ipynb` in Colab and run the first cell:

```python
!git clone https://github.com/<user>/<repo>.git
%cd <repo>/python_v2
!pip install -q numpy scipy matplotlib imageio tensorflow
!pip install -q .
```

Choose a GPU runtime (Runtime → Change runtime type → T4 / A100 / L4) for
real performance.

### 1b. Local install (Linux / WSL2)

The local environment is a conda env so the CUDA toolkit can be pinned to a
TensorFlow-compatible version.

```bash
# create the env (Python 3.10 + TF 2.15 + bundled CUDA 12.x)
conda create -n DPM python=3.10 -y
conda activate DPM
pip install 'tensorflow[and-cuda]==2.15.*'
pip install scipy matplotlib pandas imageio scikit-build-core pybind11 imageio-ffmpeg

# build + install the project (also builds the C++ extension)
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
print('GPUs:', tf.config.list_physical_devices('GPU'))
```

The C++ candidacy-manager extension is built via `scikit-build-core` (configured
in `pyproject.toml` and `CMakeLists.txt`); `pip install -e .` builds it
automatically. If you see `ModuleNotFoundError: src.simulation._candidacy_cpp`,
re-run `pip install -e .` to rebuild.

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
    Lx, Ly,                  # float — simulation box dimensions
    periodic_x   = True,     # bool  — periodic boundary in x
    periodic_y   = True,     # bool  — periodic boundary in y
    dt_factor    = 1.5,      # float — multiplier on the CFL-stable dt_max
    alpha_damp   = None,     # float — global damping override; None = auto per-particle
    E_candidates = 32,       # int   — candidacy table width; 32 emulsion / 64 dense+elastic
    skin         = 1.0,      # float — candidacy skin factor (× mean R0)
    g            = 0.0,      # float — gravitational acceleration (downward)
    # Closing-rate stability watchdog (see §14)
    disp_advisory = 0.05,    # float — one-shot warning when ρ_contact crosses
    disp_strong   = 0.10,    # float — per-chunk warning when ρ_contact crosses
    disp_critical = 0.20,    # float — divergence imminent; ρ above this means reduce dt
)
```

**Defaults updated for v2:**
- `dt_factor=1.5` for emulsion (validated up to Bo=0.10 with ~2× safety margin
  on ρ_contact). For elastic capsules use `dt_factor=1.0`. For mixed populations
  take the minimum of the per-population recommendations.
- `E_candidates=32` is sized for typical hopper / shear flows at φ ≈ 0.40.
  Dense packings or stiff/elastic systems may need `E=64`. The watchdog reports
  `max_used` after each chunk — bump if it nears the cap.

The CFL-bound `dt_max` is computed from the smallest particle's stiffest mode
(edge spring for elastic; capillary wave for emulsion). See §4.4 for the
`sys.dt` / `sys.dt_factor` sync API used to switch between absolute and
CFL-relative time-step control.

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

### 4.4a Time-step control: `sys.dt` / `sys.dt_factor` sync API

Two equivalent entry points control the integration step. Setting either
property updates the other and propagates through every TF tensor that carries
the value (`self._dt_tf`, `params['_dt_tf']`).

```python
sys.dt_factor = 1.5    # multiplicative scale on the CFL bound dt_max
sys.dt        = 0.002  # absolute dt; dt_factor is recomputed as dt / dt_max
```

Use the factor form when reasoning relative to stability (the watchdog operates
on factors); use the absolute form when matching simulated time across heterogeneous
runs. Both are safe to call between `run()` calls. The sync uses the cached
`sys._dt_max` (the CFL bound, set during `initialize()` and preserved by
`restore_state()`).

### 4.4b Closing-rate stability watchdog

Every step, the integrator returns `metrics['max_closing_ratio']` — the maximum
per-contact closing rate normalised by the contact pair's combined `r_c`. The
watchdog tracks the running max per chunk and fires advisory/strong/critical
warnings (thresholds set on the `System` constructor; see §14 for the full
description and recommended `dt_factor` defaults).

```python
print(sys._max_disp_ratio_run)   # ρ_contact peak across all run() calls so far
```

Warnings print to stdout and never halt the run — they only report. The watchdog
covers both particle-particle and particle-wall contacts.

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
| `x_all` | `(P, N, 2)` | World-frame node positions (kept in sync with x_cm/theta/u) |
| `x_cm`  | `(P, 2)` | CM positions (x, y) |
| `v_cm`  | `(P, 2)` | CM velocities |
| `theta` | `(P,)` | Rotation angles (rad) |
| `omega` | `(P,)` | Angular velocities (rad/s) |
| `u`     | `(P, N, 2)` | Node displacements in body frame |
| `u_dot` | `(P, N, 2)` | Node velocities in body frame |
| `X_ref` | `(P, N, 2)` | Body-frame rest reference geometry (immutable per particle) |
| `t`     | scalar | Simulation time |
| `step`  | scalar `int64` | Step counter |

P = total number of particles, N = nodes per particle (same for all; set by `N_nodes`).

`x_all = x_cm[:, None, :] + R(theta) @ (X_ref + u)` gives the world-frame node
positions and is maintained live by the integrator (no need to recompute
post-hoc). `X_ref` is the rest geometry — set at `initialize()` and never
modified by integration.

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
    output_path,         # str or Path — output file (.gif, .mp4, .webm, .mov, .mkv)
    fps      = 10,       # int — frames per second
    n_arc    = 6,        # int — arc interpolation points per node gap
    xlim     = None,     # (xmin, xmax) or None (auto-fit)
    ylim     = None,     # (ymin, ymax) or None (auto-fit)
    title    = None,     # str — title overlay or None
    bitrate  = None,     # int — ffmpeg bitrate in kbps for mp4/webm; None = matplotlib default
    dpi      = 120,      # int — figure dpi for rasterisation
    frames   = None,     # list — override sys.frames (e.g., render only the latest chunk)
    callback_data = None,# list — override sys.callback_data
)
```

**Output format is selected by file extension:**

- `.gif` → matplotlib `pillow` writer (lossless, large files)
- `.mp4`, `.webm`, `.mov`, `.mkv` → matplotlib `ffmpeg` writer (h264 / vp9; ~5–10×
  smaller than the equivalent GIF). The ffmpeg binary is sourced from the
  `imageio_ffmpeg` package automatically — no system `ffmpeg` install needed.

For long runs, prefer `.mp4` — a 500-frame N=1200 hopper render is ~3 MB as
mp4 vs ~20 MB as gif.

The `frames=` kwarg lets you render an arbitrary slice (e.g., live previews
during a long run, where you render only the latest chunk's snapshots) without
disturbing `sys.frames`. See `notebooks/getting_started.ipynb` for the pattern.

Walls are drawn from each registered object's primitives (`LineSegment`, `Arc`
with `angle_range`); particles are drawn as the true capsule outer contour
(Minkowski sum of the N-node polygon with the disk of radius `r_c`).

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

## 14. Closing-rate Stability Watchdog

The integrator returns, alongside each step's new state, the maximum
*closing rate* of any active contact during that step normalised by the
pair's combined contact radius:

$$
\rho_{\rm contact} = \max\!\left\{
  \max_{(a,b)\,\in\,\text{CapCand}}\,
    \frac{\| \Delta \mathbf{x}_a - \Delta \mathbf{x}_b \|}{r_{c,a} + r_{c,b}},\;
  \max_{(a,s)}\,
    \frac{\| \Delta \mathbf{x}_a \| + |\mathbf{v}_{w,s}|\,\Delta t}{r_{c,a} + r_{c,w,s}}
\right\}
$$

where the first max is over particle-particle candidate pairs and the second
over particle-wall primitive pairs. The metric is **symmetric** — bulk
translation that moves contacting particles together (e.g., free fall) leaves
$\Delta \mathbf{x}_a - \Delta \mathbf{x}_b \approx 0$ and does not trigger the
watchdog.

### 14.1 Three thresholds

| Threshold | Default | Behaviour |
|---|---:|---|
| advisory | 0.05 | One-shot warning per run. "Approaching the closing-rate limit." |
| strong   | 0.10 | Per-chunk warning. "Contacts likely to unwind in a few steps." |
| critical | 0.20 | Per-chunk warning. Numerically unstable; integration will diverge. |

Each warning prints the offending value and a recommended `dt_factor`
reduction `f_new = f_dt × (0.05 / ρ_contact)`. **The watchdog never halts
the run — it only reports.** You're responsible for restarting with a smaller
`dt_factor` if a critical warning fires.

### 14.2 Recommended `dt_factor` defaults

From the Bo×physics sweep at the hopper test bed (P=300, T=2 τ₀ matched-time
comparison):

| Working point | `dt_factor` headroom at f=1.5 | Recommended default |
|---|---|---|
| Emulsion, Bo ≤ 0.05, κ = 0.02 | ρ_contact ≈ 0.009, 5× below advisory | **1.5** |
| Emulsion, Bo = 0.10 | ρ_contact ≈ 0.025, 2× below advisory | 1.5 (mid-flow) / 1.0 (transient) |
| Elastic, ν = 0.5 | ρ_contact ≈ 0.05, at advisory | **1.0** |
| Mixed populations | take the minimum over per-population recommendations | |
| Polydisperse | take the minimum over the smallest particles | |

Transient phases (first few thousand steps after release from rest) are more
demanding than steady-state mid-flow because of larger initial relative
velocities; lower `dt_factor` for the first ~1k steps then ramp up.

See `papers/summary_of_methods/main.pdf` §11 for the full sweep tables and
discussion.

---

## 15. Parallel Sweeps — Running N Experiments Concurrently

When you want a parameter sweep (e.g. Bo at five values, or a κ × Oh matrix),
the obvious approach is `for sys in systems: sys.run(...)` — that serializes
through the GPU. `parallel_run(systems, n_steps, sample_every)` runs N
independent `System`s concurrently in one fused TF graph. After the call, each
system's `frames`, `_state`, and watchdog history are populated normally —
`make_movie`, `save_state`, etc. all work the same as a single-experiment run.

```python
from src.simulation.parallel import parallel_run

# Build N Systems however you want — could share code, could be totally different
systems = [build_test_bed(P=300, Bo=Bo)[0] for Bo in [0.05, 0.07, 0.10, 0.12, 0.15]]

# One call runs them concurrently on the GPU
parallel_run(systems, n_steps=80000, sample_every=150)

# After: each System is fully populated (frames, state, watchdog history)
for sys, Bo in zip(systems, [0.05, 0.07, 0.10, 0.12, 0.15]):
    sys.save_state(f'results/sweep_Bo{Bo}.npz')
    sys.make_movie(f'results/sweep_Bo{Bo}.mp4', xlim=(2, 28), ylim=(-2, 22))
```

### 15.1 Constraints

**Hard** (call fails otherwise):
- All systems share the same `dtype` (fp64 throughout).
- Same `n_steps` and `sample_every` for the call (kernel runs in lockstep).

**Soft** (efficiency only — XLA SIMD-fuses best when shapes match):
- Same `P` and `N_nodes` across systems.

**Free to vary per system:** `Bo`, `γ`, `κ`, `Oh`, `ν`, `dt_factor`, walls,
box, periodic config, RSA seed, particle type (emulsion / elastic / mixed).

### 15.2 When this actually pays off

The speedup comes from filling spare GPU compute capacity. On hardware where
our single-system kernels already saturate fp64 throughput (e.g. consumer
RTX cards in fp64 at our typical N), running N experiments in parallel just
costs N× the time — the GPU has no headroom to amortize. On A100/H100,
fp64 throughput is 20–65× higher than RTX 3080, so the same sweep runs with
significant headroom and shows real speedup (~3–7× depending on hardware
and problem size).

| Target | Expected 10-experiment speedup at N=300 |
|---|---:|
| RTX 3080 (this machine) | ~1× (fp64 saturated) |
| A100 (Colab) | ~5–7× |
| H100 | ~7–9× |

**Rule of thumb:** worth using on A100/H100 for sweep workloads; on consumer
GPUs the framework still works but speedup is marginal.

The N count is fixed at trace time: each unique N pays one ~10 s XLA compile
cost (cached for the rest of the session). The framework adds zero risk to
single-experiment code — `parallel_run` is purely additive.

See `notebooks/getting_started.ipynb` (last section) for a runnable 4-experiment
benchmark that prints sequential vs parallel timing and verifies trajectory
correctness.

---

## 16. Known Limitations and Areas for Development

**Candidacy manager (C++ extension, CPU-side):**
- Builds an `(K, E_candidates)` table where `E_candidates=32` by default
  (was `128` in v1). The skin-distance threshold (`skin/2`) is now wired up:
  the candidacy is only rebuilt when a particle has moved more than half the
  skin since the last rebuild — typical mid-flow workloads update every
  20–50 steps instead of every 10.
- Single-threaded. Per-experiment OpenMP parallelisation across the `parallel_run`
  loop is straightforward but not yet implemented; would help only at very high
  N_exp on data-center GPUs.

**All particles must have the same N_nodes:**
The TF state tensor `u` has shape `(P, N, 2)` with a fixed N. Mixed `N_nodes`
populations are not supported in a single System.

**Polydisperse systems with very different R0:**
The timestep `dt` is set from the smallest R0 (most constraining). Wide size
distributions force very small `dt`; the dt/dt_factor sync API (§4.4a) lets you
tune this without rebuilding the System.

**N_nodes and accuracy tradeoff:**
The ν calibration table was built at N=32. Using higher `N_nodes` (60, 120, 240)
gives smoother shapes and higher accuracy, but the calibration interpolation is
based on N=32 measurements. For quantitative work with `N_nodes != 32`, a
re-calibration sweep is recommended.

**`make_movie` for large frame counts:**
All frames are held in memory (`sys.frames`). For very long runs, memory can
become an issue. The pattern in `notebooks/getting_started.ipynb` saves the
frames to a `.npz` cache after each chunk and uses the `frames=` kwarg to render
just the latest segment for live previews — avoiding the memory build-up.

**Multi-experiment parallel runs are supported (see §15):**
Use `from src.simulation.parallel import parallel_run` to run N independent
Systems concurrently on one GPU. Speedup is hardware-dependent (~1× on RTX 3080
in fp64; ~3–7× on A100/H100). A vectorized `(N_exp, P, N, 2)` architecture
(Isaac-Gym-style) was considered and rejected — `parallel_run` achieves the same
fused-graph result with zero kernel-layer changes.

**Driven particles (`set_driven_particles`):**
Low-level API that requires `make_traj()` tuples from `tf_sim.py`. Not yet wrapped
in a clean high-level interface. Use `ParticleSpec(motion=MotionSpec(...))` instead
where possible.
