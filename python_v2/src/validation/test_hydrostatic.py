"""
Waypoint 2A.2 gate: analytic test of hydrostatic gravity implementation.

For a stationary circular droplet:
  Σ F_y^hydro = -rho_d * g * A = -m*g   (gravity, downward)
  Σ F_x^hydro = 0                         (no horizontal net force)

Both should hold to machine precision on a closed polygon.
Also tests: area penalty alone gives Σ F = 0 for circular shape at A = A0.
"""
import numpy as np
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from src.simulation.capsule_shell import CapsuleParticle

def run_test(N, R0, rho_d, g, tol=1e-8):
    p = CapsuleParticle(N=N, R0=R0, tau=0.1, S=1.0, C=100.0, rho_d=rho_d)
    # Disable elastic forces: set El_t=EI=0 so only hydrostatic acts
    p.El_t = 0.0
    p.EI   = 0.0
    p.zero_forces()
    p.accumulate_internal_forces(g=g)

    Fy_sum = p.f[:, 1].sum()
    Fx_sum = p.f[:, 0].sum()
    # Expected: -rho_d * g * A_polygon (Green's theorem on the actual polygon)
    A_polygon   = p.area()
    Fy_expected = -rho_d * g * A_polygon

    err_y = abs(Fy_sum - Fy_expected) / (abs(Fy_expected) + 1e-30)
    err_x = abs(Fx_sum)

    status_y = "PASS" if err_y < tol else "FAIL"
    status_x = "PASS" if err_x < tol else "FAIL"

    print(f"  N={N:4d}  R0={R0}  g={g:+.1f}  "
          f"Σ Fy={Fy_sum:+.6e}  expected={Fy_expected:+.6e}  "
          f"rel_err={err_y:.2e}  [{status_y}]  "
          f"Σ Fx={Fx_sum:.2e}  [{status_x}]")
    return status_y == "PASS" and status_x == "PASS"

print("=" * 80)
print("Waypoint 2A.2 — Hydrostatic gravity analytic gate")
print("=" * 80)

all_pass = True
for N in [32, 64, 120, 240]:
    for g in [1.0, -1.0, 9.81]:
        ok = run_test(N=N, R0=1.0, rho_d=1.0, g=g)
        all_pass = all_pass and ok

# Non-unit R0
for R0 in [0.5, 2.0, 5.0]:
    ok = run_test(N=120, R0=R0, rho_d=1.0, g=1.0)
    all_pass = all_pass and ok

# Non-unit rho_d
for rho in [0.1, 3.0, 10.0]:
    ok = run_test(N=120, R0=1.0, rho_d=rho, g=1.0)
    all_pass = all_pass and ok

print()
print("AREA PENALTY SELF-TEST (circular shape at A=A0 → Σ F = 0):")
p = CapsuleParticle(N=120, R0=1.0, tau=0.1, S=1.0, C=100.0, rho_d=1.0)
p.El_t = 0.0; p.EI = 0.0
p.zero_forces()
p.accumulate_internal_forces(g=0.0)
Fmag = np.linalg.norm(p.f, axis=1).max()
print(f"  Max nodal force magnitude at rest: {Fmag:.2e}  "
      f"[{'PASS' if Fmag < 1e-10 else 'FAIL'}]")

print()
print(f"Overall: {'ALL PASS ✅' if all_pass else 'FAILURES DETECTED ❌'}")
sys.exit(0 if all_pass else 1)
