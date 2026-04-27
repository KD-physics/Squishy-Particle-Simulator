# Future Directions — Performance and Scalability

This document captures architectural improvement ideas discussed after Phase 6 completion.
These are not immediate tasks — they are a roadmap of options ordered roughly by
implementation cost vs. expected payoff.

---

## Context: Current Performance Baseline

- Single simulation, unoptimized TF graph, no candidate tuning
- ~4 hrs for a 1000-particle emulsion sim on Colab A100
- Roughly matches a carefully hand-optimized MATLAB CPU-serial DEM at ~400 particles / 2 days
- GPU is absorbing significant inefficiency "for free" through parallelism
- The code was built for correctness and static-graph compliance, not performance

**Key insight:** On CPU, the bottleneck is compute (FLOPs). On GPU, the bottleneck is
memory bandwidth — the cost question is always "what am I moving through GPU memory,
and can I move less of it?"

---

## Tier 1 — Candidacy Manager Optimization (moderate effort, ~2-3x gain)

### 1a. Shrink the candidate matrix

Current: `(N_particles, N_nodes)` — every node of every particle is a candidate entry.
Target: `(N_particles, ~20)` — only the ~20 closest neighbor nodes per particle.

The win is not skipping computation (GPU runs null pairs in parallel for free) — it is
**memory bandwidth**. A smaller candidate tensor means less data moving through GPU
memory on every kernel launch. Estimated gain: 2-3x on the pair-pair kernel.

### 1b. Reduce rebuild frequency (skin-distance threshold)

Current: candidate list rebuilt every N steps (possibly too often).
Classic DEM fix: rebuild only when the fastest particle has moved more than `skin/2`
since the last rebuild. On GPU where each step is fast (~50 µs), a CPU rebuild every
step can stall the GPU. Adaptive rebuild removes that stall.

### 1c. Multithread the C++ manager

Each rebuild is independent and single-threaded. Could parallelize across particles
with `std::thread` or `openmp`. However, since the rebuild is only ~10-20% of total
wall time, the total gain is modest (~5-10%). Lower priority than 1a and 1b.

---

## Tier 2 — Subprocess Parallelism (zero code changes, ~3-6x gain)

Run M independent `System` instances in separate Python subprocesses, each with its
own TF session and candidacy manager. Coordinate from a parent process via
`multiprocessing.Pool`.

```python
from multiprocessing import Pool

def run_one(params):
    sys = System(...)
    sys.add_particles(ParticleSpec(..., nu=params['nu']))
    sys.initialize(...)
    sys.run(N, ...)
    return collect_results(sys)

with Pool(processes=8) as pool:
    results = pool.map(run_one, param_list)
```

**Why it works:** Each subprocess is fully isolated — no shared state, no contention.
The candidacy managers never see each other. The OS time-slices CPU cores naturally,
and since rebuilds are short and infrequent, contention is negligible.

**GPU sharing:** Requires `tf.config.experimental.set_memory_growth(gpu, True)` in
each subprocess so sessions allocate VRAM lazily rather than grabbing everything at init.

**Ceiling:**
- RTX 3080 (10 GB): ~6-8 parallel sims before VRAM is exhausted
- Colab A100 (40-80 GB): ~15-20 parallel sims

**Why not linear scaling like Isaac Gym:** Each subprocess launches its own independent
TF kernel streams. The GPU context-switches between them rather than seeing one batched
operation. Gain is real but bounded by GPU utilization per sim — if a single sim uses
30% of the GPU, 3 parallel sims use 90% and you get ~3x. Beyond that, diminishing returns.

---

## Tier 3 — Vectorized Multi-Experiment Architecture (significant rewrite, ~20-30x gain)

The Isaac Gym / Brax architecture applied to EPD: bake an experiment dimension `E`
directly into the TF graph so that E experiments run as one batched tensor operation.

### The idea

State tensors gain a leading experiment dimension:

| Current | Vectorized |
|---------|-----------|
| `x_all: (P, N, 2)` | `x_all: (E, P, N, 2)` |
| `v_cm: (P, 2)` | `v_cm: (E, P, 2)` |
| `candidates: (P, K)` | `candidates: (E, P, K)` |

One TF graph, one kernel launch per step, covers all E experiments simultaneously.
The GPU sees one big batched operation — no context switching, no competing sessions.

### User-facing API concept

```python
sys = System(LX, LY, n_experiments=24)
sys.add_particles(ParticleSpec(
    count=100, type='emulsion',
    nu=np.linspace(0.2, 0.9, 24),   # one value per experiment
))
sys.run(N, sample_every=100, callback=cb)
# results indexed as sys.frames[t]['x_cm']  →  shape (24, P, 2)
```

Experiments share identical graph structure (same P, N, same force types) but differ
in parameters. The experiment dimension is embarrassingly parallel by construction —
experiments never interact.

### Candidacy manager in this architecture

Each experiment has different particle positions so needs its own candidate rebuild.
Options:
- **Sequential:** E rebuilds on CPU before each GPU step (still fast if rebuild is rare)
- **Multithreaded:** E independent rebuilds on E threads (no GIL issue if C++ releases it)
- **Batched:** Rewrite manager to accept a batch dimension and output `(E, P, K)` directly

The batched option is the cleanest and enables the full throughput gain.

### Expected throughput

- 24 experiments × Tier 1 optimization × one 4-hr Colab session
- Full ν × φ parameter sweep in a single overnight run
- Iterate on science questions same-day rather than week-timescale

### Implementation cost

Significant but mechanical — not algorithmically complex, just propagating the `E`
dimension through every tensor in `system.py`, the force kernels, and the candidacy
manager. The hardest part is the manager batching. Everything else is broadcasting.

This is essentially turning the EPD codebase into the DEM equivalent of Brax —
a fully vectorized particle simulator where the experiment dimension is free.
That infrastructure is publishable in its own right.

---

## Summary Table

| Tier | Description | Code changes | Estimated gain | Effort |
|------|-------------|-------------|---------------|--------|
| 1a | Smaller candidate matrix | Manager + kernel | 2-3x per sim | Medium |
| 1b | Skin-distance rebuild threshold | Manager | included in 1a | Low |
| 1c | Multithread manager rebuilds | Manager (C++) | 5-10% total | Medium |
| 2 | Subprocess parallelism | None (wrapper only) | 3-6x throughput | Low |
| 3 | Vectorized E-dimension graph | Full system rewrite | 20-30x throughput | High |

Tiers are independent and composable. Tier 1 + Tier 2 is achievable near-term with
modest effort and gets to ~10-15x over the current baseline. Tier 3 is the long-term
architecture if the research demands high-throughput parameter sweeps at scale.
