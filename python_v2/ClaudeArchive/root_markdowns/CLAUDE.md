# CLAUDE.md — Master Instructions for Claude Code

> **Science question and design philosophy:** See `VISION.md`.
> **Execution roadmap:** See `PLAN.md`.
> This file contains technical instructions for the AI agent.

## ⛔ SCOPE — READ BEFORE ANYTHING ELSE

This project is the **EPD (Elastic Particle DEM) capsule shell model** and the
**emulsion droplet model**. It is a pure DEM (Discrete Element Method) approach.

**The FEM solver and NN-DtN surrogate were explicitly axed.** Do not:
- Run PolyFEM or any FEM solver
- Generate FEM training data (src/data_gen/ is dead weight)
- Train or use any neural network (src/nn/ is dead weight)
- Modify disk_mesh.py, fem_elastic.py, two_disk_contact.py, three_disk_contact.py,
  four_disk_contact.py, single_disk_contact.py, or polyfem_shear_benchmark.py

The live, active codebase is:
- `src/simulation/capsule_shell.py` — EPD model core (CapsuleParticle, CapsuleSim)
- `src/simulation/emulsion_particle.py` — Emulsion droplet model
- `src/simulation/ellipse_particle.py` — Ellipse/capsule shape extensions
- `src/simulation/square_particle.py` — Square particle extension
- `src/validation/` — EPD and emulsion benchmarks (all the rb_*, emulsion_*, calibration_* scripts)
- `papers/summary_of_methods/` — Methods paper for the EPD model (LaTeX source + PDF)

**Phase 2 (emulsion droplet model, all benchmarks pass) is the current completion point.**
Do not autonomously start new phases beyond what is listed in PLAN.md.

---

## 🤖 AUTONOMY RULES

- **Ask before starting new phases** not listed as pending in HANDOFF.md
- **Never run FEM solvers, polyfempy, or NN training** — these are axed
- **If a command fails**, try the documented fallback, then try one reasonable alternative.
  Only stop if all fallbacks are exhausted and the blocker is truly unresolvable.
- **When a task is complete**, update HANDOFF.md and LOG.md, then wait for instructions.
- **The only reason to stop** is a hard blocker. Write it clearly in HANDOFF.md and LOG.md.

---

## ⚡ FIRST ACTIONS ON ANY SESSION START

1. Read `HANDOFF.md` — identifies current phase, last command, immediate next step
2. Read the relevant phase section of `PLAN.md` — understand the waypoints you're working toward
3. Scan the last 5–10 entries in `LOG.md` — understand recent progress and any open issues
4. Then execute the **Immediate Next Step** listed in `HANDOFF.md`

Do not start fresh. Do not re-run completed steps. Resume exactly where the last session ended.

---

## PROJECT OVERVIEW

**Goal:** Build and validate a 2D **EPD (Elastic Particle DEM)** model in which each disk
is a capsule shell — a closed elastic membrane — with a tunable squishiness knob ν.
The model spans the full behavioral range from glass-like (stiff, compressible) to
rubber-like (soft, area-conserving).

The EPD model is a **purely DEM approach**: no FEM interior solves, no NN surrogates.
Each disk's shape is represented by N perimeter nodes connected by spring-like edge,
bending, and area forces. Contact between particles uses a penalty-based repulsion.
Rigid-body dynamics (translation + rotation) are integrated explicitly.

A related emulsion droplet model uses the same membrane mechanics with different
physics (surface tension replaces edge springs; no interior pressure).

**The science:** How does particle deformability (squishiness) affect collective dynamics,
rearrangements, and stress relaxation in dense 2D packings near the jamming transition?
The EPD model provides a computationally cheap, physically motivated tool to explore this.

---

## EPD MODEL ARCHITECTURE

### Particle: CapsuleParticle

Each particle is a 2D closed elastic membrane (capsule shell) described by N perimeter
nodes. The forces on each node are:

| Force term | Parameter | Role |
|-----------|-----------|------|
| Edge spring | El_t (membrane stiffness) | Resists stretching of perimeter |
| Bending spring | τ | Controls contact flatness (turning angle at contact) |
| Area spring | K_area = q × El_t | Resists area change (incompressibility) |
| Contact penalty | C | Prevents node overlap between particles / wall |
| Damping | α_damp | Dissipates energy (capsule ringing) |

### The Squishiness Knob

The single user-level parameter is **ν** (effective Poisson ratio, 0.18–0.94).
Everything else is derived:

```python
q = np.exp(np.interp(nu_target, nu_arr, np.log(q_arr)))  # from calibration table
TAU    = np.sqrt(12.0 * 0.2)    # ≈ 1.5492, fixed at b=0.2 working point
K_area = q * (12.0 * S / TAU**2)
C      = 3000.0 * S * (1.0 + q)
alpha_damp = 2.0 / T_wave
```

Calibration table: `results/calibration_sweep/calibration_data.json` (18 q values, N=32).
See `memory/project_nu_knob.md` for the full recipe.

### Rigid Body Integration

Each particle has 3 rigid-body DOFs: x_cm, y_cm, θ (rotation).
The step_rb() integrator (in capsule_shell.py) decomposes contact impulses into
rigid-body impulses exactly, with machine-precision momentum conservation.
Elastic internal forces are computed from node displacements relative to the
rigid-body frame.

---

## CANONICAL PARAMETER RECIPE

For any (τ, q) with S=1, R0=1:

| Parameter | Formula | Notes |
|-----------|---------|-------|
| El_t_base | 12 / τ² | Membrane stiffness at S=R0=1 |
| K_area    | q × El_t_base | Area spring |
| C         | 3000 × S × (1 + q) | Contact hardness — τ-INDEPENDENT |
| α₀        | 2.0 | Damping; ζ ≈ 16% constant for all (τ,q) |
| dt_factor | 0.4 | Auto-stable for all (τ,q) |

The b=0.2 working point: TAU = √(12×0.2) ≈ 1.549, which gives flat contacts
and physically motivated contact flatness vs. squishiness behavior.

---

## EMULSION MODEL

The emulsion droplet model (src/simulation/emulsion_particle.py) uses the same
perimeter node framework but with:
- Surface tension γ instead of edge springs (constant tension along arc)
- No bending stiffness (τ=0 in emulsion limit)
- Fluid pressure inside droplet (area spring → Laplace pressure)

Validated against: capillary-wave damping (ω₂ = √6), T1-event three-droplet
squeeze, and falling droplet with physical Bond number.

---

## GROUND TRUTH: 2D HERTZIAN CONTACT

Used for Phase 0 validation. These are 2D plane strain equations.

**Effective radius:** `R_eff = R/2` for equal disks

**Contact half-width:** `a = sqrt(4 * F * R_eff * (1 - nu²) / (π * E))`

**Approach displacement:**
`δ ≈ (2F(1-ν²)/πE) × (ln(4R_eff/a) - 0.5)`

---

## HANDOFF PROTOCOL — MANDATORY AFTER EVERY MILESTONE

After completing any waypoint or significant step:

1. **Update HANDOFF.md:**
   - Set `Current Phase` and `Status`
   - Record `Last successful command` (exact command or script name)
   - Write `Immediate next step`
   - Note any blockers or issues

2. **Append to LOG.md** (append to top, never delete):
   ```
   ### [YYYY-MM-DD HH:MM] — Phase X.Y: Brief description
   **Action:** What was done
   **Command:** Key command or script
   **Result:** Outcome
   **Metrics:** Any numbers (errors, counts, timings)
   **Next:** Immediate next action
   ```

3. **Update PLAN.md:** Check off completed waypoints with ✅

---

## ENVIRONMENT SETUP

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy matplotlib pandas h5py tqdm
pip install scikit-fem  # FEM library for Hertz benchmark only
pip install imageio     # for movie rendering
```

Verify: `python -c "import numpy, scipy, skfem, h5py; print('OK')"`

**Do NOT install:** polyfempy, triangle, meshio, torch — these are for the axed NN/FEM path.

---

## FILE ORGANIZATION

```
project/
├── CLAUDE.md                    ← YOU ARE HERE
├── HANDOFF.md                   ← Read first on session start
├── PLAN.md                      ← Waypoints and gate metrics
├── LOG.md                       ← Append-only milestone log
├── VISION.md                    ← Science question and philosophy
│
├── src/
│   ├── simulation/
│   │   ├── capsule_shell.py     ← EPD model core (ACTIVE)
│   │   ├── emulsion_particle.py ← Emulsion model (ACTIVE)
│   │   ├── ellipse_particle.py  ← Ellipse/capsule shapes (ACTIVE)
│   │   ├── square_particle.py   ← Square particle (ACTIVE)
│   │   ├── contact_primitives.py← Contact geometry (ACTIVE)
│   │   └── two_disk.py          ← Two-disk squeeze (ACTIVE)
│   │
│   │   [DEAD — FEM/NN approach, do not use]
│   │   ├── disk_mesh.py
│   │   ├── fem_elastic.py
│   │   ├── two_disk_contact.py
│   │   ├── three_disk_contact.py
│   │   ├── four_disk_contact.py
│   │   ├── single_disk_contact.py
│   │   └── polyfem_shear_benchmark.py
│   │
│   ├── validation/              ← EPD and emulsion benchmarks (ACTIVE)
│   │   └── [dtn_*.py files are DEAD — NN/FEM path]
│   │
│   ├── nn/                      ← DEAD — NN-DtN surrogate (axed)
│   └── data_gen/                ← DEAD — FEM training data generation (axed)
│
├── papers/summary_of_methods/   ← EPD model methods paper (ACTIVE)
│   ├── main.tex / main.pdf      ← Main paper
│   └── supplement.tex / .pdf   ← Supplementary material
│
├── results/                     ← EPD simulation results, plots, movies
└── data/
    └── processed/               ← DEAD — FEM/NN training data (ignore)
```

---

## PERMISSIONS

Claude Code has autonomy to:
- Read, write, delete any file in this workspace
- Install packages via pip (EPD-relevant packages only)
- Run EPD/emulsion simulations and validation scripts
- Generate plots and movies

Claude Code must NOT autonomously:
- Start phases beyond what HANDOFF.md lists as the immediate next step
- Run FEM solvers, NN training, or data generation for NN
- Modify the paper (papers/) without being asked
