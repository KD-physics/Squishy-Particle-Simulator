"""
contact_primitives.py — Rigid contact primitive objects for disk-boundary contact.

Each primitive exposes:
    gap, normal = primitive.gap_and_normal(xy)
        xy     : (2,) array — query point (disk perimeter node position)
        gap    : scalar — signed distance; gap < 0 means penetrating
        normal : (2,) unit vector pointing FROM primitive surface TOWARD disk interior
                 (i.e. the direction a repulsive force should push the disk node)

Penalty force on node = k_pen * max(0, -gap) * normal

Primitive types
---------------
LineSegment  : finite flat wall segment defined by two endpoints
               inward_normal must be supplied (which side the disk is on)
Arc          : circular arc; convex=True means the arc bulges toward the disk
               (like a rigid circular punch / indenter pressing in)
               convex=False means the disk is inside the circle (enclosure / pipe wall)
Polygon      : ordered list of vertices forming a closed boundary
               corner handling: closest-point projection — each query point resolves
               against exactly one feature (segment or vertex), no double counting
               Convenience subclass: Box(cx, cy, w, h)
"""

import numpy as np


# ── helpers ───────────────────────────────────────────────────────────────────

def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-15 else v


def _closest_point_on_segment(p, a, b):
    """
    Closest point on segment a→b to query point p.
    Returns (closest_point, t) where t ∈ [0,1].
    """
    ab = b - a
    ab2 = np.dot(ab, ab)
    if ab2 < 1e-30:
        return a.copy(), 0.0
    t = np.clip(np.dot(p - a, ab) / ab2, 0.0, 1.0)
    return a + t * ab, t


# ── LineSegment ───────────────────────────────────────────────────────────────

class LineSegment:
    """
    Finite flat wall segment.

    Parameters
    ----------
    p0, p1 : array-like (2,)
        Endpoints of the segment.
    inward_normal : array-like (2,)
        Unit normal pointing TOWARD the disk (the interior side).
        Must be perpendicular to (p1 - p0); caller is responsible.
        If None, computed automatically as the left-hand normal of p0→p1
        (assumes disk is to the left when walking from p0 to p1).
    """

    def __init__(self, p0, p1, inward_normal=None):
        self.p0 = np.asarray(p0, dtype=float)
        self.p1 = np.asarray(p1, dtype=float)
        seg = self.p1 - self.p0
        if inward_normal is None:
            # left-hand normal of p0→p1
            self.normal = _unit(np.array([-seg[1], seg[0]]))
        else:
            self.normal = _unit(np.asarray(inward_normal, dtype=float))
        # Capsule radius of this primitive (0 = zero-thickness line, set >0 for thick walls)
        self.r_c = 0.0

    def gap_and_normal(self, xy):
        """
        gap  = signed distance from xy to segment (positive = outside / safe side)
        normal = inward_normal (constant for a flat wall)
        """
        xy = np.asarray(xy, dtype=float)
        cp, _ = _closest_point_on_segment(xy, self.p0, self.p1)
        diff = xy - cp                         # vector from surface to query point
        gap = np.dot(diff, self.normal)        # positive if on the safe (inward) side
        return gap, self.normal.copy()

    def gap_and_normal_batch(self, points):
        """
        Vectorized gap_and_normal for (M, 2) array of query points.
        Returns gaps (M,) and normals (M, 2).
        """
        points = np.asarray(points, dtype=float)  # (M, 2)
        M = len(points)
        # Closest point on segment for each query point
        ab  = self.p1 - self.p0
        ab2 = np.dot(ab, ab)
        if ab2 < 1e-30:
            cp = np.tile(self.p0, (M, 1))
        else:
            t   = np.clip(np.dot(points - self.p0, ab) / ab2, 0.0, 1.0)  # (M,)
            cp  = self.p0 + t[:, np.newaxis] * ab                          # (M, 2)
        diff = points - cp                                                   # (M, 2)
        gaps = np.dot(diff, self.normal)                                    # (M,)
        normals = np.tile(self.normal, (M, 1))                             # (M, 2)
        return gaps, normals

    def __repr__(self):
        return f"LineSegment({self.p0}, {self.p1})"


# ── Arc ───────────────────────────────────────────────────────────────────────

class Arc:
    """
    Circular arc primitive.

    Parameters
    ----------
    center : array-like (2,)
        Centre of the circle.
    radius : float
        Radius of the arc.
    angle_range : (float, float) or None
        (theta_min, theta_max) in radians; None means full circle.
        Angles measured CCW from +x axis.
    convex : bool
        True  → arc surface faces outward; disk is OUTSIDE the circle.
                 Gap = dist(xy, center) - radius   (positive when disk is outside)
                 Normal = unit vector from center toward xy (pushes disk outward)
                 Use for: rigid circular punch / indenter pressing in from outside.
        False → arc surface faces inward; disk is INSIDE the circle.
                 Gap = radius - dist(xy, center)   (positive when disk is inside)
                 Normal = unit vector from xy toward center (pushes disk inward)
                 Use for: pipe walls, cylindrical enclosures.
    """

    def __init__(self, center, radius, angle_range=None, convex=True):
        self.center = np.asarray(center, dtype=float)
        self.radius = float(radius)
        self.angle_range = angle_range   # None = full circle
        self.convex = convex

    def _in_arc(self, xy):
        """Return True if the angle from center to xy falls within angle_range."""
        if self.angle_range is None:
            return True
        theta = np.arctan2(xy[1] - self.center[1], xy[0] - self.center[0])
        t0, t1 = self.angle_range
        # Normalise to [t0, t0+2π)
        while theta < t0:
            theta += 2 * np.pi
        return theta <= t1

    def gap_and_normal(self, xy):
        xy = np.asarray(xy, dtype=float)
        diff = xy - self.center
        dist = np.linalg.norm(diff)

        if dist < 1e-15:
            # Degenerate: node exactly at center
            n = np.array([1.0, 0.0])
            return (self.radius if not self.convex else -self.radius), n

        radial = diff / dist   # unit vector from center to xy

        if self.convex:
            # Disk is outside; arc faces outward
            # gap > 0 when node is beyond the arc (safe); < 0 when penetrating
            gap = dist - self.radius
            normal = radial          # push node away from center
        else:
            # Disk is inside; arc faces inward
            gap = self.radius - dist
            normal = -radial         # push node toward center

        return gap, normal

    def __repr__(self):
        kind = "convex" if self.convex else "concave"
        return f"Arc(center={self.center}, R={self.radius}, {kind})"


# ── Polygon ───────────────────────────────────────────────────────────────────

class Polygon:
    """
    Closed polygon boundary (rigid wall).

    Vertices should be ordered so that the disk interior is to the LEFT
    when traversing the boundary counter-clockwise (standard CCW convention).
    Equivalently: the inward normal of each edge points toward the disk.

    Corner handling: closest-point projection across all segments and vertices.
    Each query point maps to exactly one nearest feature → no double-counting.

    Parameters
    ----------
    vertices : list of array-like (2,)
        Ordered polygon vertices (CCW with disk inside).
    """

    def __init__(self, vertices):
        verts = [np.asarray(v, dtype=float) for v in vertices]
        n = len(verts)
        assert n >= 3, "Polygon requires at least 3 vertices"

        self._verts = verts
        self._segments = []   # list of (p0, p1, inward_normal)

        for i in range(n):
            p0 = verts[i]
            p1 = verts[(i + 1) % n]
            seg = p1 - p0
            # Left-hand normal (CCW convention → inward when disk is inside CCW polygon)
            normal = _unit(np.array([-seg[1], seg[0]]))
            self._segments.append((p0, p1, normal))

    def gap_and_normal(self, xy):
        """
        Find the active wall for convex-polygon SDF and return its gap and normal.

        For a convex polygon, the signed-distance function equals the MAX gap over
        all half-planes:
          SDF > 0  →  node outside (no force)
          SDF ≤ 0  →  node inside  (force from the argmax wall)

        Using argmax(gap) instead of argmin(dist) prevents spurious forces on
        outside nodes near corners, where the nearest-by-distance segment can
        have a negative gap even though the node is geometrically outside.

        gap   = max signed gap across all segments (positive = outside/safe)
        normal = wall normal at the argmax segment
        """
        xy = np.asarray(xy, dtype=float)
        best_gap = -np.inf
        best_normal = np.array([1.0, 0.0])

        for (p0, p1, seg_normal) in self._segments:
            cp, t = _closest_point_on_segment(xy, p0, p1)
            diff = xy - cp
            gap = np.dot(diff, seg_normal)   # signed: positive = outside for CW polygon

            if gap > best_gap:
                best_gap = gap
                best_normal = seg_normal

        return best_gap, best_normal.copy()

    def __repr__(self):
        return f"Polygon({len(self._segments)} sides)"


# ── Convenience: Box ──────────────────────────────────────────────────────────

class Box(Polygon):
    """
    Axis-aligned rectangular box.
    Disk is assumed to be INSIDE the box.

    Parameters
    ----------
    cx, cy : float  — centre of the box
    w, h   : float  — full width and height
    """

    def __init__(self, cx, cy, w, h):
        hw, hh = w / 2.0, h / 2.0
        # CCW order (disk inside): bottom-left → bottom-right → top-right → top-left
        vertices = [
            [cx - hw, cy - hh],
            [cx + hw, cy - hh],
            [cx + hw, cy + hh],
            [cx - hw, cy + hh],
        ]
        super().__init__(vertices)
        self.cx, self.cy, self.w, self.h = cx, cy, w, h

    def __repr__(self):
        return f"Box(cx={self.cx}, cy={self.cy}, w={self.w}, h={self.h})"


# ── Convenience: Hopper ───────────────────────────────────────────────────────

class Hopper:
    """
    Two inward-facing line segments forming a V-shaped hopper (open at top).

    The hopper apex is at (apex_x, apex_y).  Two walls rise at ±half_angle
    from vertical (i.e. half_angle=0 means vertical walls, =45° means 90° V).

    The disk is assumed to be between the two walls.

    Parameters
    ----------
    apex     : (float, float) — tip of the V
    half_angle : float — half-opening angle in radians (0 = vertical walls)
    height   : float — height of the walls above apex
    """

    def __init__(self, apex, half_angle, height):
        ax, ay = float(apex[0]), float(apex[1])
        ha = float(half_angle)
        h  = float(height)

        # Left wall: runs from apex upward-left
        dx_l = -np.sin(ha)
        dy_l =  np.cos(ha)
        self._left  = LineSegment(
            [ax, ay],
            [ax + dx_l * h, ay + dy_l * h],
            inward_normal=_unit(np.array([np.cos(ha), np.sin(ha)]))   # points right
        )
        # Right wall: runs from apex upward-right
        dx_r =  np.sin(ha)
        dy_r =  np.cos(ha)
        self._right = LineSegment(
            [ax, ay],
            [ax + dx_r * h, ay + dy_r * h],
            inward_normal=_unit(np.array([-np.cos(ha), np.sin(ha)]))  # points left
        )
        self._apex = np.array([ax, ay])

    def gap_and_normal(self, xy):
        """
        Returns gap/normal for the closer of the two walls.
        The apex vertex is the Voronoi tie-breaker: at the apex the closest
        wall is determined by which side of the bisector the node is on.
        """
        xy = np.asarray(xy, dtype=float)
        g_l, n_l = self._left.gap_and_normal(xy)
        g_r, n_r = self._right.gap_and_normal(xy)

        cp_l, _ = _closest_point_on_segment(xy, self._left.p0,  self._left.p1)
        cp_r, _ = _closest_point_on_segment(xy, self._right.p0, self._right.p1)
        d_l = np.linalg.norm(xy - cp_l)
        d_r = np.linalg.norm(xy - cp_r)

        if d_l <= d_r:
            return g_l, n_l
        else:
            return g_r, n_r

    def __repr__(self):
        return f"Hopper(apex={self._apex})"


# ── Unit tests ────────────────────────────────────────────────────────────────

def _run_unit_tests():
    import sys
    fails = 0

    def check(name, cond, detail=''):
        nonlocal fails
        if cond:
            print(f"  PASS  {name}")
        else:
            print(f"  FAIL  {name}  {detail}")
            fails += 1

    print("contact_primitives unit tests")
    print("=" * 50)

    # ── LineSegment ──
    print("\nLineSegment:")
    seg = LineSegment([0, 0], [1, 0], inward_normal=[0, 1])
    g, n = seg.gap_and_normal([0.5, 0.3])
    check("point above segment: gap > 0", g > 0, f"gap={g:.4f}")
    check("point above segment: normal=[0,1]", np.allclose(n, [0, 1], atol=1e-10))

    g, n = seg.gap_and_normal([0.5, -0.1])
    check("point below segment: gap < 0", g < 0, f"gap={g:.4f}")

    g, n = seg.gap_and_normal([0.5, 0.0])
    check("point on segment: gap=0", abs(g) < 1e-12, f"gap={g:.2e}")

    # Point off end (closest point is p1)
    g, n = seg.gap_and_normal([1.5, 0.2])
    check("point off end: gap finite", np.isfinite(g))

    # ── Arc (convex) ──
    print("\nArc (convex=True, indenter):")
    arc = Arc(center=[2, 0], radius=1.0, convex=True)
    # Point at distance 1.5 from center → gap = 0.5 (outside)
    g, n = arc.gap_and_normal([0.5, 0])   # dist = 1.5 from center [2,0]
    check("outside arc: gap > 0", g > 0, f"gap={g:.4f}")
    # normal should point from center toward node (radial outward = away from center)
    # node=[0.5,0], center=[2,0] → (node-center)=[-1.5,0], normal=[-1,0]
    # dot product of normal with (node-center) should be > 0 (same direction)
    check("outside arc: normal points away from center",
          np.dot(n, np.array([0.5, 0]) - np.array([2, 0])) > 0,
          f"n={n}")

    # Point at distance 0.5 from center → gap = -0.5 (penetrating)
    g, n = arc.gap_and_normal([1.5, 0])   # dist = 0.5 from center [2,0]
    check("inside arc: gap < 0", g < 0, f"gap={g:.4f}")

    # Point on arc surface
    g, n = arc.gap_and_normal([1.0, 0])   # dist = 1.0 from center [2,0]
    check("on arc surface: gap=0", abs(g) < 1e-12, f"gap={g:.2e}")

    # ── Arc (concave) ──
    print("\nArc (convex=False, enclosure):")
    enc = Arc(center=[0, 0], radius=2.0, convex=False)
    g, n = enc.gap_and_normal([1.0, 0])
    check("inside enclosure: gap > 0", g > 0, f"gap={g:.4f}")

    g, n = enc.gap_and_normal([2.5, 0])
    check("outside enclosure: gap < 0", g < 0, f"gap={g:.4f}")

    g, n = enc.gap_and_normal([2.0, 0])
    check("on enclosure wall: gap=0", abs(g) < 1e-12, f"gap={g:.2e}")

    # ── Polygon (Box) ──
    print("\nBox (Polygon):")
    box = Box(cx=0, cy=0, w=4, h=4)

    g, n = box.gap_and_normal([0, 0])
    check("interior center: gap > 0", g > 0, f"gap={g:.4f}")

    g, n = box.gap_and_normal([1.5, 0])
    check("near right wall: gap > 0", g > 0, f"gap={g:.4f}")

    g, n = box.gap_and_normal([2.5, 0])
    check("outside right wall: gap < 0", g < 0, f"gap={g:.4f}")

    g, n = box.gap_and_normal([2.0, 0])
    check("on right wall: gap=0", abs(g) < 1e-12, f"gap={g:.2e}")

    # Corner: point outside at corner
    g, n = box.gap_and_normal([2.5, 2.5])
    check("outside corner: gap < 0", g < 0, f"gap={g:.4f}")

    # Corner continuity: points on either side of corner apex should give
    # smoothly varying (no discontinuous jump) normals
    pts = [[2.01, y] for y in np.linspace(-2.0, 2.0, 5)]
    gaps = [box.gap_and_normal(p)[0] for p in pts]
    check("right wall gap monotone along wall",
          all(abs(g - gaps[2]) < 0.02 for g in gaps[1:-1]),
          f"gaps={[f'{g:.3f}' for g in gaps]}")

    # ── Hopper ──
    print("\nHopper:")
    hop = Hopper(apex=(0, 0), half_angle=np.radians(30), height=2.0)
    # Point above apex on centreline: should be inside (gap > 0 from both walls)
    g, n = hop.gap_and_normal([0, 1.0])
    check("centreline above apex: gap > 0", g > 0, f"gap={g:.4f}")

    # Point at apex: gap near 0
    g, n = hop.gap_and_normal([0.0, 0.0])
    check("at apex: gap ≈ 0", abs(g) < 1e-10, f"gap={g:.2e}")

    print(f"\n{'All tests PASS' if fails == 0 else f'{fails} test(s) FAILED'}")
    return fails == 0


if __name__ == '__main__':
    import sys
    ok = _run_unit_tests()
    sys.exit(0 if ok else 1)
