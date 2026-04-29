---
name: nu_as_user_knob
description: ν (Poisson ratio) is the single user-level material parameter for EPD simulations; all other params are derived from it via the calibration table
type: project
---

## ν is the only user-level material knob

The EPD model reduces to a single physically interpretable parameter: the effective
Poisson ratio ν, which controls particle "squishiness" from area-compressible
(ν ≈ 0.18) to near-incompressible (ν ≈ 0.94, limited by area conservation).

**Valid range:** ν ∈ [0.18, 0.94]  (N=32, b=0.2 working point)

**Why:** All other model parameters (q, τ, C, K_area, α) are either fixed at the
b=0.2 working point or derived from ν. The user specifies ν; the code looks up q
from the calibration table and computes everything else.

---

## How to go from ν_target to simulation parameters

**Step 1 — look up q:**
```python
import json, numpy as np
data = json.load(open('results/calibration_sweep/calibration_data.json'))
rows = [(r['q'], r['metrics']['0.08']['nu'])
        for r in data if r['N'] == 32]
rows.sort()
q_arr  = np.array([r[0] for r in rows])
nu_arr = np.array([r[1] for r in rows])
q = np.exp(np.interp(nu_target, nu_arr, np.log(q_arr)))   # log-interp
```

**Step 2 — derive all simulation parameters:**
```python
B_TARGET  = 0.2
TAU       = np.sqrt(12.0 * B_TARGET)     # ≈ 1.5492  (fixed working point)
El_t_base = 12.0 * S / TAU**2            # membrane stiffness at S=1, R0=1
K_area    = q * El_t_base                # area penalty
C         = 3000.0 * S * (1.0 + q)      # contact stiffness (τ-independent)
alpha_damp = ALPHA0 / T_wave             # where ALPHA0=2.0 gives Q≈3
```

**Step 3 — pass to CapsuleParticle:**
```python
p = CapsuleParticle(N=N, R0=R0, tau=TAU, S=S, C=C, rho_f=RHO_F, K_area=K_area)
```

---

## Second knob: damping (α / Q factor)

Damping α is the second user-level knob, set independently of ν.

Design formula: α = ω₂_bend / Q_target = (6/√5)·√(S/(ρ_f·τ)) / (Q_target·R₀)

In practice, use ALPHA0 = α·T_wave = 2.0 (Q≈3) as the default. This keeps the
damping ratio ζ = ALPHA0/(4π) ≈ 16% constant regardless of τ or q.

COR strongly depends on α because T_contact ≈ 1.5·T₂_bend — see
project_alpha_damp_calibration.md.

---

## Calibration data location

`results/calibration_sweep/calibration_data.json`
- 18 q-values: 0.05 → 50 (N=32, b=0.2, τ≈1.549)
- Metrics at ε_ref ∈ {0.04, 0.06, 0.08, 0.10, 0.12}
- Use ε_ref = 0.08 as the reference for q → ν lookup

ν table at ε_ref=0.08 (N=32):
  q=0.05  → ν≈0.18  |  q=0.10  → ν≈0.21  |  q=0.30  → ν≈0.33
  q=0.50  → ν≈0.42  |  q=1.00  → ν≈0.56  |  q=2.00  → ν≈0.69
  q=5.00  → ν≈0.83  |  q=10.0  → ν≈0.88  |  q=50.0  → ν≈0.94

**Why:** ν encapsulates the combined effect of membrane stiffness and area penalty
on the ratio of lateral to axial strain. It is the physically meaningful output that
connects EPD to real-material Poisson ratios and to other DEM models.
