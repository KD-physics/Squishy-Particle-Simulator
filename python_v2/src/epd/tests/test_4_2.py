"""
test_4_2.py — Phase 4.2: particles.py verification.

Tests:
  1. ParticleSpec(count=5, nu=0.5): q from calibration table, TAU, K_area, C
  2. ParticleSpec(count=5, q=X, tau_b=0.2): low-level override, same formulas
  3. Polydisperse (poly_dist=0.05): mean=1.0, std/mean≈0.05
  4. Monotone: nu=0.3 (stiff) has q < nu=0.8 (soft)
  5. Rigid type: frozen_shape=True, large C
  6. build() returns CapsuleParticle list
  7. extra_forces stored correctly

Run: python src/epd/tests/test_4_2.py
"""

import sys, os
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)

from src.epd.particles import ParticleSpec, nu_to_q, q_to_nu, _load_calibration

results = []

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════════════
# Test 1 — elastic by nu
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 1: ParticleSpec(count=5, type='elastic', nu=0.5) ───────────────")

spec1 = ParticleSpec(count=5, type='elastic', nu=0.5)
d1    = spec1.derived

# q from calibration table at nu=0.5
nu_arr, q_arr = _load_calibration()
q_expected = float(np.exp(np.interp(0.5, nu_arr, np.log(q_arr))))

check("1.1 q matches calibration table",
      abs(d1['q'] - q_expected) / q_expected < 0.01,
      f"expected={q_expected:.4f} got={d1['q']:.4f}")

TAU_expected = np.sqrt(12.0 * 0.2)
check("1.2 TAU = sqrt(12*0.2)",
      abs(d1['TAU'] - TAU_expected) < 1e-10,
      f"expected={TAU_expected:.6f} got={d1['TAU']:.6f}")

El_t_expected = 12.0 / TAU_expected**2
check("1.3 El_t = 12/TAU²",
      abs(d1['El_t'] - El_t_expected) < 1e-10,
      f"expected={El_t_expected:.6f} got={d1['El_t']:.6f}")

K_area_expected = d1['q'] * El_t_expected
check("1.4 K_area = q * El_t",
      abs(d1['K_area'] - K_area_expected) < 1e-10,
      f"expected={K_area_expected:.6f} got={d1['K_area']:.6f}")

C_expected = 3000.0 * (1.0 + d1['q'])
check("1.5 C = 3000*(1+q)",
      abs(d1['C'] - C_expected) < 1e-6,
      f"expected={C_expected:.2f} got={d1['C']:.2f}")

check("1.6 alpha = 2.0 (default)",
      abs(d1['alpha'] - 2.0) < 1e-10)

print(f"  (derived params: q={d1['q']:.4f}, TAU={d1['TAU']:.4f}, "
      f"K_area={d1['K_area']:.4f}, C={d1['C']:.1f})")


# ══════════════════════════════════════════════════════════════════════════════
# Test 2 — elastic by q directly
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 2: ParticleSpec(count=5, q=<same q>, tau_b=0.2) ─────────────────")

q_val = spec1.derived['q']
spec2 = ParticleSpec(count=5, q=q_val, tau_b=0.2)
d2    = spec2.derived

check("2.1 q passed directly = q from nu",
      abs(d2['q'] - q_val) < 1e-10,
      f"expected={q_val:.6f} got={d2['q']:.6f}")
check("2.2 Same K_area as nu-derived",
      abs(d2['K_area'] - d1['K_area']) < 1e-8)
check("2.3 Same C as nu-derived",
      abs(d2['C'] - d1['C']) < 1e-6)
check("2.4 nu stored as None when q given directly",
      spec2.nu is None)


# ══════════════════════════════════════════════════════════════════════════════
# Test 3 — polydisperse size distribution
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 3: ParticleSpec(count=100, poly_dist=0.05, seed=42) ────────────")

spec3  = ParticleSpec(count=500, poly_dist=0.05)
R0_arr = spec3.sample_R0(seed=42)

check("3.1 R0_arr mean = 1.0 exactly",
      abs(R0_arr.mean() - 1.0) < 1e-14,
      f"mean={R0_arr.mean():.15f}")

std_frac = R0_arr.std() / R0_arr.mean()
check("3.2 R0 std/mean ≈ 0.05 (within 30%)",
      abs(std_frac - 0.05) / 0.05 < 0.30,
      f"std/mean={std_frac:.4f}")

check("3.3 All R0 > 0",
      np.all(R0_arr > 0),
      f"min_R0={R0_arr.min():.4f}")

# Shorthand vs dict form should give same result
spec3b  = ParticleSpec(count=500, poly_dist={'type': 'gaussian', 'sigma': 0.05})
R0_dict = spec3b.sample_R0(seed=42)
check("3.4 Shorthand 0.05 == dict form",
      np.allclose(R0_arr, R0_dict),
      f"max_diff={np.max(np.abs(R0_arr-R0_dict)):.2e}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 4 — monotone q vs nu
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 4: Monotone q(nu) ──────────────────────────────────────────────")

spec_stiff = ParticleSpec(count=10, nu=0.3)
spec_soft  = ParticleSpec(count=10, nu=0.8)

q_stiff = spec_stiff.derived['q']
q_soft  = spec_soft.derived['q']

check("4.1 q(nu=0.3) < q(nu=0.8) (monotone)",
      q_stiff < q_soft,
      f"q_stiff={q_stiff:.4f} q_soft={q_soft:.4f}")

print(f"  nu=0.3 → q={q_stiff:.4f};  nu=0.8 → q={q_soft:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 5 — rigid type
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 5: ParticleSpec(count=3, type='rigid') ─────────────────────────")

spec_r = ParticleSpec(count=3, type='rigid')
check("5.1 rigid → frozen_shape=True", spec_r.frozen_shape)
check("5.2 rigid → C = 1e6",
      abs(spec_r.derived['C'] - 1e6) < 1e-6,
      f"C={spec_r.derived['C']:.2e}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 6 — build() produces CapsuleParticle list
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 6: build() produces CapsuleParticle list ───────────────────────")

from src.simulation.capsule_shell import CapsuleParticle

spec4 = ParticleSpec(count=4, nu=0.5, N_nodes=32)
centers = [(float(i) * 3.0, 0.0) for i in range(4)]
particles = spec4.build(seed=42, centers=centers)

check("6.1 build returns list of 4", len(particles) == 4)
check("6.2 All are CapsuleParticle", all(isinstance(p, CapsuleParticle) for p in particles))
check("6.3 Particles have correct N_nodes",
      all(p.N == 32 for p in particles),
      f"N={[p.N for p in particles]}")

# R0 mean ≈ 1.0 (monodisperse)
R0s = [p.R0 for p in particles]
check("6.4 Monodisperse R0 = 1.0",
      np.allclose(R0s, 1.0, atol=1e-10),
      f"R0s={R0s}")

# CM positions set correctly
cm_x = [float(p.x_cm[0]) for p in particles]
check("6.5 CMs at specified centers",
      np.allclose(cm_x, [0, 3, 6, 9], atol=1e-10),
      f"cm_x={cm_x}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 7 — extra_forces dict
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 7: extra_forces stored correctly ───────────────────────────────")

spec5 = ParticleSpec(count=5, extra_forces={'drag': 0.1, 'activity': 0.05})
check("7.1 extra_forces stored",
      spec5.extra_forces == {'drag': 0.1, 'activity': 0.05})
check("7.2 extra_forces defaults to empty dict",
      ParticleSpec(count=5).extra_forces == {})
check("7.3 accessing missing key returns 0.0",
      spec5.extra_forces.get('gravity', 0.0) == 0.0)


# ══════════════════════════════════════════════════════════════════════════════
# Test 8 — bimodal and explicit distributions
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 8: poly_dist bimodal and explicit ───────────────────────────────")

spec_bi = ParticleSpec(count=20, poly_dist={'type': 'bimodal', 'ratio': 0.5, 'delta': 0.1})
R0_bi   = spec_bi.sample_R0(seed=42)
check("8.1 Bimodal mean=1.0", abs(R0_bi.mean() - 1.0) < 1e-14)
unique  = np.unique(np.round(R0_bi, 6))
check("8.2 Bimodal has 2 unique sizes",
      len(unique) == 2, f"unique={unique}")

vals = list(np.linspace(0.8, 1.2, 10))
spec_ex = ParticleSpec(count=10, poly_dist={'type': 'explicit', 'values': vals})
R0_ex   = spec_ex.sample_R0(seed=0)
check("8.3 Explicit values mean=1.0", abs(R0_ex.mean() - 1.0) < 1e-14)
# Normalised vals match
raw = np.array(vals); raw /= raw.mean()
check("8.4 Explicit values correct",
      np.allclose(R0_ex, raw), f"max_diff={np.max(np.abs(R0_ex-raw)):.2e}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary table
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Parameter table ─────────────────────────────────────────────────────")
for nu_test in [0.3, 0.5, 0.7, 0.9]:
    try:
        sp = ParticleSpec(count=1, nu=nu_test)
        d  = sp.derived
        print(f"  nu={nu_test:.1f}: q={d['q']:7.4f}  K_area={d['K_area']:8.4f}  C={d['C']:8.1f}")
    except Exception:
        print(f"  nu={nu_test:.1f}: (out of range)")

n_pass = sum(r[1] for r in results)
n_fail = sum(not r[1] for r in results)
print(f"\n{'='*60}")
print(f"TOTAL: {n_pass}/{len(results)} PASS")
if n_fail:
    print("FAILED tests:")
    for name, passed, detail in results:
        if not passed:
            print(f"  - {name}  ({detail})")

sys.exit(0 if n_fail == 0 else 1)
