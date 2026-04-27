"""
test_4_1.py — Phase 4.1: objects.py + motion.py verification.

Tests:
  1. Box resolved(t=0): 4 primitives, inward normals, region_polygon, point-in-polygon
  2. CouetteCell resolved(t=0): 2 arc primitives, contains() for annulus
  3. Box with set_motion(MotionSpec(vx=1.0)): origin shifts by dx=vx*t
  4. RegularPolygon (hexagon): 6 walls, inward normals, region_polygon
  5. CircleObstacle: 1 arc, exclusion=interior

Run: python src/epd/tests/test_4_1.py
Output: results/phase41_objects/objects_test.png
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MPoly
from matplotlib.collections import PatchCollection

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
sys.path.insert(0, ROOT)

from src.epd.objects import (Box, Channel, CouetteCell, CircleObstacle,
                               RegularPolygon, CustomObject, _point_in_polygon)
from src.epd.motion import MotionSpec
from src.simulation.contact_primitives import LineSegment, Arc

OUTDIR = os.path.join(ROOT, 'results', 'phase41_objects')
os.makedirs(OUTDIR, exist_ok=True)

results = []

def check(name, cond, detail=''):
    tag = 'PASS' if cond else 'FAIL'
    print(f"  [{tag}] {name}" + (f"  — {detail}" if detail else ''))
    results.append((name, bool(cond), detail))
    return bool(cond)


# ══════════════════════════════════════════════════════════════════════════════
# Test 1 — Box
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 1: Box(width=10, height=8, x0=5, y0=4) ─────────────────────────")

box = Box(width=10, height=8, x0=5.0, y0=4.0, exclusion='exterior')
prims = box.resolved(t=0.0)

check("1.1 Box has 4 primitives", len(prims) == 4, f"n={len(prims)}")
all_line = all(isinstance(d['prim'], LineSegment) for d in prims)
check("1.2 All Box primitives are LineSegment", all_line)

# Normals should point inward: dot(normal, centre - midpoint) > 0
cx, cy = 5.0, 4.0
inward_ok = True
for d in prims:
    seg = d['prim']
    mid = 0.5 * (seg.p0 + seg.p1)
    to_centre = np.array([cx, cy]) - mid
    if np.dot(seg.normal, to_centre) <= 0:
        inward_ok = False
        print(f"    normal {seg.normal} at mid {mid} not inward")
check("1.3 All normals inward", inward_ok)

poly = box.region_polygon(t=0.0)
check("1.4 region_polygon returns 4 vertices", poly is not None and len(poly['vertices']) == 4,
      f"n={len(poly['vertices']) if poly else 'None'}")
check("1.5 region_polygon exclusion='exterior'", poly['exclusion'] == 'exterior')

# Points inside box (x in [0,10], y in [0,8])
inside_pts = [[5.0, 4.0], [1.0, 1.0], [9.0, 7.0], [3.0, 6.0], [7.0, 2.0]]
outside_pts = [[12.0, 4.0], [5.0, 10.0], [-1.0, 4.0], [5.0, -1.0],
               [-0.5, 4.0], [10.5, 4.0]]     # clearly outside

in_ok  = all(box.contains(p) for p in inside_pts)
out_ok = all(not box.contains(p) for p in outside_pts)
check("1.6 interior points pass contains()", in_ok,
      f"failed: {[p for p in inside_pts if not box.contains(p)]}")
check("1.7 exterior points fail contains()", out_ok,
      f"failed: {[p for p in outside_pts if box.contains(p)]}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 2 — CouetteCell
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 2: CouetteCell(inner=2, outer=5, x0=0, y0=0) ───────────────────")

cc = CouetteCell(inner_radius=2.0, outer_radius=5.0, x0=0.0, y0=0.0)
prims_cc = cc.resolved(t=0.0)

check("2.1 CouetteCell has 2 primitives", len(prims_cc) == 2, f"n={len(prims_cc)}")
all_arc = all(isinstance(d['prim'], Arc) for d in prims_cc)
check("2.2 Both CouetteCell primitives are Arc", all_arc)

radii = sorted([d['prim'].radius for d in prims_cc])
check("2.3 Arc radii are inner=2, outer=5",
      abs(radii[0] - 2.0) < 1e-10 and abs(radii[1] - 5.0) < 1e-10,
      f"radii={radii}")

# Contains: annulus r in (2, 5)
check("2.4 Point at r=3 is inside annulus", cc.contains([3.0, 0.0]),
      f"contains([3,0])={cc.contains([3,0])}")
check("2.5 Point at r=1 (inside inner) is outside", not cc.contains([1.0, 0.0]),
      f"contains([1,0])={cc.contains([1,0])}")
check("2.6 Point at r=6 (outside outer) is outside", not cc.contains([6.0, 0.0]),
      f"contains([6,0])={cc.contains([6,0])}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 3 — Box with motion
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 3: Box with MotionSpec(vx=1.0).resolved(t=2.0) ─────────────────")

box2 = Box(width=4.0, height=4.0, x0=0.0, y0=0.0)
ms = MotionSpec(vx=1.0)
box2.set_motion(ms)

prims_t0 = box2.resolved(t=0.0)
prims_t2 = box2.resolved(t=2.0)

# At t=2, all x-coordinates should be shifted by +2
# Compute centroid of all endpoints at t=0 and t=2
def seg_cx(prims):
    xs = []
    for d in prims:
        seg = d['prim']
        if isinstance(seg, LineSegment):
            xs.extend([seg.p0[0], seg.p1[0]])
    return np.mean(xs)

cx0 = seg_cx(prims_t0)
cx2 = seg_cx(prims_t2)
check("3.1 Box shifted by dx=2.0 at t=2", abs(cx2 - cx0 - 2.0) < 1e-10,
      f"cx0={cx0:.3f} cx2={cx2:.3f} Δ={cx2-cx0:.4f}")

# Velocity at t=0 should be vx=1
vx, vy, om = ms.velocity(t=0.0)
check("3.2 MotionSpec velocity vx=1 at t=0", abs(vx - 1.0) < 1e-12)

# Displacement at t=3
dx, dy, dth = ms.displacement(t=3.0)
check("3.3 MotionSpec displacement dx=3 at t=3", abs(dx - 3.0) < 1e-12)

# Wall velocity in resolved prim = (1, 0)
vel0 = prims_t2[0]['vel']
check("3.4 Wall velocity from resolved = (1.0, 0.0)",
      np.allclose(vel0, [1.0, 0.0]), f"vel={vel0}")

# region_polygon at t=2 should be centred at x=2
poly2 = box2.region_polygon(t=2.0)
cx_poly = np.mean(poly2['vertices'][:, 0])
check("3.5 region_polygon centred at x=2 at t=2",
      abs(cx_poly - 2.0) < 1e-10, f"cx_poly={cx_poly:.4f}")


# ══════════════════════════════════════════════════════════════════════════════
# Test 4 — RegularPolygon (hexagon)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 4: RegularPolygon(sides=6, radius=3) ────────────────────────────")

hex_obj = RegularPolygon(sides=6, radius=3.0, x0=0.0, y0=0.0, exclusion='interior')
prims_h = hex_obj.resolved(t=0.0)

check("4.1 Hexagon has 6 primitives", len(prims_h) == 6, f"n={len(prims_h)}")

# All normals point inward (toward origin)
inward_h = all(np.dot(d['prim'].normal, -0.5*(d['prim'].p0+d['prim'].p1)) > 0
               for d in prims_h)
check("4.2 Hexagon normals inward", inward_h)

poly_h = hex_obj.region_polygon(t=0.0)
check("4.3 Hexagon region_polygon has 6 vertices",
      poly_h is not None and len(poly_h['vertices']) == 6)

# Centre (0,0) is inside hexagon
check("4.4 Centre is inside hexagon polygon",
      _point_in_polygon([0, 0], poly_h['vertices']))
# exclusion='interior' → contains() should return False for interior point
check("4.5 contains([0,0]) = False (interior excluded)",
      not hex_obj.contains([0, 0]))
# Outside point
check("4.6 contains([10,0]) = True (exterior allowed)",
      hex_obj.contains([10, 0]))


# ══════════════════════════════════════════════════════════════════════════════
# Test 5 — CircleObstacle
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 5: CircleObstacle(radius=2, x0=3, y0=3) ────────────────────────")

circ = CircleObstacle(radius=2.0, x0=3.0, y0=3.0, exclusion='interior')
prims_c = circ.resolved(t=0.0)

check("5.1 CircleObstacle has 1 primitive", len(prims_c) == 1)
check("5.2 CircleObstacle primitive is Arc",
      isinstance(prims_c[0]['prim'], Arc))
check("5.3 Arc center at (3,3)",
      np.allclose(prims_c[0]['prim'].center, [3, 3], atol=1e-10))
check("5.4 Arc radius=2",
      abs(prims_c[0]['prim'].radius - 2.0) < 1e-10)
check("5.5 contains([3,3]) = False (interior excluded)",
      not circ.contains([3, 3]))
check("5.6 contains([0,0]) = True (exterior allowed)",
      circ.contains([0, 0]))


# ══════════════════════════════════════════════════════════════════════════════
# Test 6 — to_make_prim_list round-trip
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Test 6: to_make_prim_list → make_prim_data round-trip ───────────────")

import src.simulation.tf_sim as tf_sim_mod
import tensorflow as tf
tf_sim_mod.set_dtype(tf.float64)
from src.simulation.tf_sim import make_prim_data, DTYPE

box3 = Box(width=20.0, height=20.0, x0=10.0, y0=10.0)
prim_list = box3.to_make_prim_list(t=0.0)
prim_data = make_prim_data(prim_list)

check("6.1 make_prim_data succeeds with Box prim_list", prim_data is not None)
n_segs = prim_data['seg_p0'].shape[0]
check("6.2 prim_data has 4 segments", n_segs == 4, f"n_segs={n_segs}")


# ══════════════════════════════════════════════════════════════════════════════
# Summary figure
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot 1: Box with test points
ax = axes[0]
poly_box = box.region_polygon()['vertices']
ax.add_patch(MPoly(poly_box, fill=False, edgecolor='navy', lw=2))
for p in inside_pts:
    ax.plot(p[0], p[1], 'go', ms=8)
for p in outside_pts:
    ax.plot(p[0], p[1], 'rx', ms=8, mew=2)
ax.set_xlim(-2, 14); ax.set_ylim(-2, 10)
ax.set_aspect('equal'); ax.set_title('Box: green=inside, red=outside')
ax.grid(True, alpha=0.3)

# Plot 2: CouetteCell
ax = axes[1]
theta_arr = np.linspace(0, 2*np.pi, 200)
ax.plot(5*np.cos(theta_arr), 5*np.sin(theta_arr), 'navy', lw=2, label='outer')
ax.plot(2*np.cos(theta_arr), 2*np.sin(theta_arr), 'navy', lw=2, ls='--', label='inner')
test_pts = [(3, 0, True), (1, 0, False), (6, 0, False), (0, 3, True), (0, 4.9, True)]
for (px, py, exp) in test_pts:
    c = cc.contains([px, py])
    col = 'go' if c else 'rx'
    ax.plot(px, py, col, ms=8, mew=2)
ax.set_xlim(-7, 7); ax.set_ylim(-7, 7)
ax.set_aspect('equal'); ax.set_title('CouetteCell: green=annulus, red=outside')
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Plot 3: Hexagon + moved box
ax = axes[2]
poly_h_v = hex_obj.region_polygon()['vertices']
ax.add_patch(MPoly(poly_h_v, fill=False, edgecolor='darkorange', lw=2, label='hexagon'))
ax.plot(0, 0, 'rx', ms=8, mew=2, label='inside→excluded')
ax.plot(4, 0, 'go', ms=8, label='outside→allowed')

# Show moved box at t=2
poly2_v = box2.region_polygon(t=2.0)['vertices']
ax.add_patch(MPoly(poly2_v, fill=False, edgecolor='teal', lw=2, ls='--', label='Box@t=2'))
ax.plot(0, 0, 'k+', ms=10, mew=2, label='box origin @t=0')

ax.set_xlim(-5, 8); ax.set_ylim(-5, 6)
ax.set_aspect('equal'); ax.set_title('Hexagon + Box with motion')
ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

plt.tight_layout()
outpath = os.path.join(OUTDIR, 'objects_test.png')
plt.savefig(outpath, dpi=120)
print(f"\nSummary plot: {outpath}")

# Final
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
