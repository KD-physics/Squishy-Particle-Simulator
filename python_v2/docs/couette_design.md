# Couette Cell Geometry — Design Recipe

Practical guide for picking `R_inner`, `R_outer`, and particle count `P` for a
2-D Couette shear experiment in EPD.  Used by `getting_started.ipynb`'s
Test C and by the T1-tracking workflow in `results/test_C_couette_T1/`.

---

## Physics motivation

In an annular Couette geometry, shear stress decays as **1/r²** with distance
from the inner (driven) wall — force per unit length is conserved across each
annular shell.  As a result:

- Strain rate is **highest at the inner wall** and decays exponentially into
  the bulk.
- A **shear band** localises near the moving wall.  Its thickness is set by
  the material (typically **5–10 droplet diameters** for 2-D dense
  emulsions / granular).
- The **outer wall** is mostly a passive boundary — making the gap wider
  doesn't widen the band, it just adds undeformed bulk beyond it.
- A **rough wall** (driven inner ring of frozen particles) is the standard
  way to inject strain without slip.  An equally rough **outer pinned ring**
  is needed because frictionless walls would let the whole annulus rotate
  rigidly under the inner drive.

---

## Design knobs

Given particle radius `R0` (the simulator normalises mean R0 to exactly 1.0),
two natural geometric scales:

- **`N_inner = R_inner / R0`** — inner radius in droplet radii.
  Sets curvature of the inner wall.
- **`N_width = (R_outer − R_inner) / (2·R0)`** — gap thickness in droplet
  **diameters**.  Sets shear-band coverage.

These plus the target packing fraction `phi` determine the particle count `P`,
or vice versa.

---

## Constraint formula

With effective droplet radius `r_eff = 1 + 2π/N_nodes` (≈ 1.105 at N=60 and
1.196 at N=32):

$$
\phi \; = \; \frac{P \cdot r_\text{eff}^2}{4 \cdot N_\text{width} \cdot (N_\text{inner} + N_\text{width})}
$$

(Derivation: annulus area = `π·[(N_inner + 2 N_width)² − N_inner²]·R0² = 4π·N_width·(N_inner + N_width)·R0²`,
particle area = `P·π·r_eff²·R0²`, ratio = phi.)

**Four knobs `{N_inner, N_width, phi, P}`, one constraint, so any three set
the fourth:**

| Fix three of | Solve for |
|---|---|
| `N_inner, N_width, phi` | **P** = 4 · phi · N_width · (N_inner + N_width) / r_eff² |
| `N_inner, N_width, P` | **phi** = P · r_eff² / [4 · N_width · (N_inner + N_width)] |
| `N_inner, P, phi` | **N_width** = ½·[−N_inner + √(N_inner² + P · r_eff²/phi)] |
| `N_width, P, phi` | **N_inner** = P · r_eff² / (4 · phi · N_width) − N_width |

---

## Choosing `N_inner` and `N_width`

### Inner-wall curvature (`N_inner`)

The number of droplets that fit tangentially around the inner ring is
roughly `π·(N_inner + 1)`.  If `N_inner` is too small the inner ring is
discrete and curves strongly into the bulk; too large and the geometry
approaches planar Couette.

| `N_inner` | Inner-wall droplets | Regime |
|---:|---:|---|
| 3 | ~13 | strong curvature, sparse ring |
| 5 | ~19 | conventional minimum for "rough wall" |
| **8–10** | **~28–35** | **standard for shear-band studies** |
| 15 | ~50 | approaching planar Couette |

### Gap thickness (`N_width`)

Determines what fraction of the gap is shear-band vs. bulk.

| `N_width` (diameters) | Use |
|---:|---|
| 1–2 | "Narrow-gap" rheometers (DIN 53019, ISO 3219) — uniform shear stress |
| 5–8 | Shear band only; no bulk margin |
| **10–15** | **Resolves shear band + a few diameters of bulk** |
| > 20 | Wide-gap; weak boundary effect — band well-isolated |

### Packing fraction (`phi`)

| `phi` | Regime |
|---:|---|
| ≤ 0.40 | Loose, gas-like; little contact |
| 0.50–0.70 | Marginal; intermittent contacts |
| 0.80–0.83 | Jammed but compressible (use `kappa` for emulsion stiffness) |
| **0.84–0.86** | **Random close packing (2D); standard for plasticity** |
| > 0.88 | Highly deformed; near hexagonal limit |

---

## Reference recipes

All at `phi = 0.84`, `N_nodes = 60` (`r_eff² = 1.221`):

| Target | `N_inner` | `N_width` | `R_ratio` | **P** | Notes |
|---|---:|---:|---:|---:|---|
| Minimal shear band | 5 | 10 | 5.0 | **412** | strong inner curvature; tight bulk margin |
| Balanced | 8 | 12 | 4.0 | **658** | reasonable curvature, full band + bulk |
| **Textbook shear band** | **10** | **15** | **4.0** | **1030** | flat-ish inner wall, clean band + bulk |
| Wide-gap luxury | 10 | 20 | 5.0 | **1646** | overkill for band; lots of bulk |
| Planar-Couette limit | 15 | 15 | 3.0 | **1238** | low curvature throughout |
| **Test C reference** | **8.1** | **9.5** | **3.3** | **500** | what the included notebook example uses |

Where `R_ratio = R_outer / R_inner = (N_inner + 2·N_width) / N_inner`.

---

## Code snippet

```python
import numpy as np

def couette_geometry(N_inner, N_width, phi, N_nodes=60, R0=1.0):
    """Return (R_inner, R_outer, P, R_ratio) for a Couette cell."""
    r_eff   = 1.0 + 2.0 * np.pi / N_nodes
    R_inner = N_inner * R0
    R_outer = R_inner + 2.0 * N_width * R0
    P       = int(round(4.0 * phi * N_width * (N_inner + N_width) / r_eff**2))
    R_ratio = R_outer / R_inner
    return R_inner, R_outer, P, R_ratio

# Textbook shear-band run
R_in, R_out, P, ratio = couette_geometry(N_inner=10, N_width=15, phi=0.84)
print(f"R_inner={R_in:.1f}  R_outer={R_out:.1f}  P={P}  R_ratio={ratio:.2f}")
# R_inner=10.0  R_outer=40.0  P=1030  R_ratio=4.00
```

---

## Practical notes

1. **Start with `N_inner` and `N_width`** for a shear-band study — they directly
   correspond to physical features (curvature, band coverage).  Compute `P`
   last.
2. **`R_ratio` is a derived diagnostic**, not a design knob.  Quote it for
   comparison with the rheometer literature (typical ranges: narrow-gap 1.05–1.15,
   wide-gap granular 2.0–5.0).
3. **For large P** (> ~2000) on consumer GPUs in fp64, runs become
   substantial (~10 hr for 500k driven steps at P=2000 on RTX 3080).  See
   the `parallel_run` workflow in `getting_started.ipynb` for sweeps and
   the dt/E tuning notes in `README.md` §14 for managing per-step cost.
4. **Outer wall must be pinned** for any driven run, otherwise the
   frictionless system co-rotates rigidly.  Identify outer-ring particles
   by `r_cm > R_outer − 2·R0` and pass them to `set_driven_particles` with
   a zero-velocity trajectory and `frozen=True`.

---

## Further reading

- DIN 53019 / ISO 3219 — narrow-gap rheometer geometry standards
- Pouliquen & Forterre (2002), *J. Fluid Mech.* — granular Couette flow
- Mueth et al. (2000), *Nature* — shear band thickness in 3D granular
- Bocquet, Colin & Ajdari (2009), *PRL* — soft-glass non-local rheology in Couette
