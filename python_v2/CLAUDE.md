# Deformable Particle Model (DPM) — Claude Code Operating Instructions

## Session Startup Checklist
1. Read `ClaudeArchive/CLAUDE.md` and `ClaudeArchive/HANDOFF.md` for codebase background and architecture
2. Read `HANDOFF.md` for current session state and next target
3. Read `PLAN.md` for current phase breakdown and task status
4. Activate conda env: `conda activate DPM`
5. Check in with Ken to confirm scope before beginning any work

---

## Project Overview
EPD (Elastic Perimeter Disk) — a 2D discrete element simulator for deformable soft particles,
parameterized by a single dimensionless squishiness parameter ν (effective Poisson ratio).
Built on TensorFlow (GPU-accelerated physics loop via tf.while_loop) with a C++ pybind11
extension (CandidacyManager) for neighbor search.

See `ClaudeArchive/` for full architectural context. That archive is historical reference only —
do not treat it as an instruction set for this session.

---

## Environment
- **Conda env:** `DPM` (Python 3.10, TensorFlow 2.15, CUDA 12.x)
- **Target platforms:** Local RTX 3080 (development), Google Colab A100 (scale testing)
- **Colab compatibility is a hard requirement** — no dependencies or syntax that doesn't run on Colab Python 3.10 / TF 2.15
- All work done within the DPM env; no other envs

---

## Hard Constraints — Never Violate

### Physics are Sacred
Force calculation logic cannot change. This includes:
- Edge spring forces (perimeter elasticity)
- Bending forces
- Area spring (volume conservation)
- Contact penalty forces (inter-capsule)
- Surface tension / capillary pressure (emulsion model)

### Stack is Fixed
TensorFlow + C++ (pybind11) only. Never introduce:
- PyTorch or any other ML framework
- FEM libraries (polyfempy, etc.)
- meshio, torch, or any previously-axed dependency
- Any library not already in requirements.txt without explicit approval

### No Workarounds — Fix the Class
If a method or class is broken, **patch the method/class itself** in a backward-compatible way.
The fix must live in the class permanently so all callers benefit.
Never bypass a class with a one-off workaround in a notebook cell or script.
If a plotting method doesn't work, fix the renderer. If a save method fails, fix the checkpoint.
The notebook cells define the intended API surface — they are the spec, not the thing to rewrite.

### Backward Compatibility
All patches must not break existing callers. When adding functionality to a class,
use default arguments so existing call sites continue to work unchanged.

---

## What Can Change (Optimization Targets)
- Time integration scheme (e.g., upgrading to velocity Verlet for larger dt)
- Candidacy list management: rebuild frequency, skin distance threshold, matrix shape
- TF graph structure and compilation strategy
- C++ neighbor search: threading, BVH improvements, multi-experiment indexing
- Batch experiment dimension (N_exp) — Phase 4

---

## Working Style

This is a **collaborative workflow**, not fully autonomous execution:

1. **Plan together:** Ken and Claude scope the next target through conversation, then Claude updates `PLAN.md` with the agreed breakdown
2. **Execute autonomously:** Ken says "work phase X" or "work items Y–Z" and Claude works to that stopping point without interruption
3. **Patch autonomously:** Minor patches and backward-compatible fixes → just do it, document in `HANDOFF.md`
4. **Discuss before restructuring:** If a fix requires non-trivial restructuring (not a minor patch), surface it to Ken before touching anything
5. **Update `HANDOFF.md`** at the end of every work session and at natural stopping points

### Decision Threshold
| Action | Behavior |
|--------|----------|
| Bug fix, backward-compatible patch | Autonomous — just do it |
| Adding a method or argument to a class | Autonomous — just do it |
| Profiling, running tests, benchmarking | Autonomous |
| Any change to time integration or candidacy logic | Autonomous (these are valid targets) |
| Restructuring a class or module | Discuss with Ken first |
| Changing a function signature in a breaking way | Discuss with Ken first |
| Changes touching multiple modules simultaneously | Discuss with Ken first |

---

## Key Files
| File | Purpose |
|------|---------|
| `VISION.md` | Phase roadmap — stable reference, don't modify |
| `PLAN.md` | Living task tracker — update as work progresses |
| `HANDOFF.md` | Session state — update at end of every session |
| `ClaudeArchive/CLAUDE.md` | Architectural reference — read-only, historical |
| `ClaudeArchive/HANDOFF.md` | Old session state — historical reference only |
| `notebooks/reproduce_paper.ipynb` | Primary Phase 1 target |
| `src/simulation/tf_sim.py` | Core TF physics loop |
| `src/simulation/candidacy_manager.py` | C++ neighbor search wrapper |
| `src/cpp/candidacy_manager.hpp` | C++ CandidacyManager implementation |

---

## Running Tests
```bash
conda activate DPM
cd /path/to/python
python src/epd/tests/test_4_1.py   # objects
python src/epd/tests/test_4_2.py   # motion specs
python src/epd/tests/test_4_3.py   # ParticleSpec / calibration
python src/epd/tests/test_4_4.py   # System initialization
python src/epd/tests/test_4_5.py   # recording + callbacks
python src/epd/tests/test_4_6.py   # checkpoint
```
Always run the relevant test after patching a class to confirm no regressions.
