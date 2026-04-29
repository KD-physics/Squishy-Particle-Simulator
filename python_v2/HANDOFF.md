# DPM — Session Handoff

_Update this file at the end of every session and at natural stopping points.
This file is the single source of truth for picking up work on any machine or in any new session._

---

## Environment
| Item | Value |
|------|-------|
| Machine | Local RTX 3080, Windows 11 + WSL2 Ubuntu |
| Conda env | `DPM` (to be created — Phase 0) |
| Python target | 3.10 |
| TF target | 2.15 (CUDA 12.x) |
| Project path (WSL) | `/mnt/c/Users/kend1/OneDrive/Desktop/DesktopFolder/Research/SoftParticles/Deformable-Particle-Model/python` |
| Project path (Windows) | `C:\Users\kend1\OneDrive\Desktop\DesktopFolder\Research\SoftParticles\Deformable-Particle-Model\python` |
| Claude Code | Installed via npm (`@anthropic-ai/claude-code` v2.1.119), authenticated via claude.ai subscription |

---

## Current State
**Active phase:** Phase 4a ✅ COMPLETE — parallel_run prototype + benchmarked
**Last completed:** `src/simulation/parallel.py` (parallel_run for N concurrent Systems via XLA-fused N step_full_tf calls); 4-system benchmark at P=300 on RTX 3080 confirmed correctness (drift < 3e-4 r_c, dominated by fp64 ordering in chaotic exp); ~1× speedup on RTX 3080 due to fp64 saturation. Section added to `notebooks/getting_started.ipynb`. 120k overnight run is paused at step 20000 (resumable from `full_N1200_Bo007_120k_state.npz`).
**Next action:** Run `getting_started.ipynb` parallel_run cell on Colab A100/H100 to validate predicted 3–7× speedup, then resume the overnight 120k.

---

## Last Session Summary (2026-04-29)
- Replaced absolute-displacement watchdog with symmetric per-contact **closing-rate** ρ_contact (`tf_sim.step_full_tf`):
  - Particle--particle: `|Δx_a − Δx_b| / (r_c_a + r_c_b)` per CapCandidate pair
  - Particle--wall: `(|Δx_a| + |v_wall|·dt) / (r_c_a + r_c_w_s)` per (node, primitive)
  - `make_prim_data` now precomputes `seg_v_max`, `arc_v_max` (max wall velocity per primitive: trans + |ω|·arm + |A·ω|).
- Added `dt`/`dt_factor` sync API on `System`: setting either property updates the other and propagates to `_dt_tf` and `params['_dt_tf']` (the JIT path reads from params, so both must be updated). `_dt_max` (CFL bound) stored in `initialize()` and `restore_state()` keeps `dt_factor` in sync.
- Re-ran all three dt sweeps with the new criterion (T=2.0 sim-time matched comparison, drift vs REF=0.4):
  - Bo=0.05 emul: ρ_contact ramps 0.002→0.026 from f=0.4→2.0; CRITICAL at f=3.0 (ρ=0.47).
  - Bo=0.10 emul: 0.007→0.025 from f=0.4→1.5; STRONG at f=2.0 (ρ=0.17); CRITICAL at f=3.0.
  - Elastic ν=0.5: 0.014→0.053 from f=0.4→1.5 (advisory); NaN at f=3.0.
  Findings agree with prior absolute-disp data because the hopper's max is wall-dominated (seg_r_c = 0 makes new = old in that limit). The new criterion is symmetric and will report differently in pure-PP regimes (bulk shear, rotating drum).
- Paper §11 updated: closing-rate equation, watchdog cover (PP + PW), dt/dt_factor sync API, refreshed tables, recommendations now state min-of-populations rule for mixed/polydisperse. Recompiled main.pdf (20 pages, 2.2 MB).

### Phase 0 Verification Logs
- `nvidia-smi`: Driver 576.02, CUDA 12.9 runtime, RTX 3080 10GB
- `tf.config.list_physical_devices('GPU')`: `[PhysicalDevice(name='/physical_device:GPU:0', device_type='GPU')]`
- WSL NUMA warnings present in TF startup output — harmless; standard for WSL2 (kernel without NUMA support)
- cuDNN/cuFFT/cuBLAS "factory already registered" warnings present — harmless; cosmetic TF 2.15 issue

---

## Known Issues / Gotchas

### WSL DNS
WSL DNS requires `nameserver 8.8.8.8` in `/etc/resolv.conf`.
Persists via `/etc/wsl.conf` with `generateResolvConf = false`.
`cdn.anthropic.com` does not resolve via WSL's default DNS; Claude Code installed via npm as workaround.
If DNS breaks again: `sudo rm /etc/resolv.conf && echo "nameserver 8.8.8.8" | sudo tee /etc/resolv.conf`

### reproduce_paper.ipynb
Cells were assembled from earlier pre-TF-API scripts. Some cells work, some don't.
Phase 1 works through them sequentially. Failures are fixed by patching classes, not cells.

### ClaudeArchive
Contains historical session docs from the original Hetzner development session.
Read for background/architecture only — not instructions for this session.

---

## Decisions Log
_Record any non-trivial decisions made during autonomous work sessions._

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-04-27 | Targeting TF 2.15 + Python 3.10 | Colab compatibility; avoids TF 2.16+ Keras 3 migration |
| 2026-04-27 | Claude Code installed via npm | WSL DNS cannot resolve cdn.anthropic.com |
| 2026-04-27 | TF installed via `tensorflow[and-cuda]==2.15.*` (pip-bundled CUDA) | Avoids dependency on system CUDA toolkit; matches Colab's pip-managed CUDA install pattern |
| 2026-04-29 | Stability metric is per-contact **closing rate**, not per-node displacement | Symmetric in who's moving: bulk translation (everyone falling together) leaves contact distance unchanged and so should not trigger the watchdog. Old criterion was over-sensitive to free-fall velocity. |
| 2026-04-29 | dt and dt_factor are equivalent entry points, kept in sync | Users alternate between absolute (`dt=0.002`) and CFL-relative (`dt_factor=1.5`) reasoning. Either setter updates the other and re-emits TF tensors; both `self._dt_tf` and `params['_dt_tf']` must be updated because the JIT path reads from params. |
| 2026-04-29 | Default `dt_factor = 1.5` for emulsion; `1.0` for elastic; min for mixed | From sweep data: Bo=0.05 has 5× advisory headroom at f=1.5; elastic at ν=0.5 hits advisory at f=1.5, so 1.0 leaves margin. Mixed/polydisperse use min over populations. |
| 2026-04-29 | Phase 4a uses N hardcoded `step_full_tf` calls inside one jit-compiled `tf.while_loop`, NOT a leading `N_exp` tensor dim | Same XLA-fused outcome with zero kernel-layer changes. User builds N independent Systems → `parallel_run(systems, ...)` → each system fully populated as before. Cost: one ~10s XLA recompile per distinct N (cached for the session). |
| 2026-04-29 | parallel_run speedup is hardware-dependent; ~1× on RTX 3080 in fp64, predicted 3–7× on A100/H100 | RTX 3080 has only 0.47 TFLOPS fp64 (1/64 of fp32). Single-experiment kernel at our N saturates that. A100 has 9.7 TFLOPS fp64, H100 has 30+; same kernel uses ~5% of fp64 capacity, leaving room for batched experiments to amortize. Framework is correct on RTX 3080 but only gives speedup on data-center GPUs. |

---

## Quick Resume (New Session)
```bash
# In WSL:
cd /mnt/c/Users/kend1/OneDrive/Desktop/DesktopFolder/Research/SoftParticles/Deformable-Particle-Model/python
conda activate DPM
claude
# Then: read HANDOFF.md → read PLAN.md → check in with Ken
```
