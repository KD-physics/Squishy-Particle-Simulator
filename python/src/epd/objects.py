"""
objects.py — Simulation object hierarchy for the EPD user API.

Classes
-------
SimulationObject          base class: exclusion semantics, motion, resolved()
  Wall                    single line-segment primitive
  ArcWall                 single arc primitive
  CompositeObject         list of children + origin/rotation transform
    Box                   4 walls (fully closed), exclusion='exterior'
    Channel               top+bottom walls only
    CouetteCell           outer arc (container) + inner arc (obstacle)
    CircleObstacle        1 arc, exclusion='interior'
    RegularPolygon        N walls
    CustomObject          user-assembled composite

resolved(t) contract
--------------------
Returns a flat list of primitive dicts.  Each dict:
    {
      'prim'    : LineSegment or Arc  (world-frame position at time t),
      'k_pen'   : float,
      'vel'     : (2,) ndarray  — translational wall velocity at time t,
      'omega'   : float         — spin of wall about r_ref,
      'r_ref'   : (2,) ndarray,
      'r_c_wall': float,
      'exclusion': str or None,
    }

Use to_make_prim_list(resolved_list) to convert to the tuples expected by
make_prim_data() in tf_sim.py.
"""

import numpy as np
from src.simulation.contact_primitives import LineSegment, Arc


# ── point-in-polygon test ────────────────────────────────────────────────────

def _point_in_polygon(pt, vertices):
    """
    Ray-casting point-in-polygon test.
    vertices : (N, 2) array of polygon vertex coordinates (closed implicitly).
    Returns True if pt is strictly inside the polygon.
    """
    x, y  = float(pt[0]), float(pt[1])
    verts = np.asarray(vertices, dtype=float)
    n     = len(verts)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = verts[i]
        xj, yj = verts[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-300) + xi):
            inside = not inside
        j = i
    return inside


# ── base class ────────────────────────────────────────────────────────────────

class SimulationObject:
    """
    Base class for all simulation boundary / obstacle objects.

    Parameters
    ----------
    kind : str — descriptive label
    """

    def __init__(self, kind):
        self._kind       = kind
        self._exclusion  = None    # 'interior', 'exterior', or None
        self._motion     = None    # MotionSpec or None
        self._k_pen      = 1.0
        self._r_c_wall   = 0.0

    def set_exclusion(self, mode):
        """'interior', 'exterior', or None."""
        assert mode in ('interior', 'exterior', None)
        self._exclusion = mode
        return self

    def set_motion(self, motion_spec):
        """Attach a MotionSpec (from motion.py)."""
        self._motion = motion_spec
        return self

    def set_k_pen(self, k):
        """Wall penalty multiplier (default 1.0)."""
        self._k_pen = float(k)
        return self

    def set_r_c_wall(self, r):
        """Wall capsule radius (default 0)."""
        self._r_c_wall = float(r)
        return self

    def exclusion_area(self, t=0.0):
        """
        Area (m²) excluded from particle placement by this object.
        Used to correct the phi denominator: accessible_area = Lx*Ly - sum(exclusion_area).
        Default 0; override in obstacle subclasses with exclusion='interior'.
        """
        return 0.0

    def set_render(self, color=None, linewidth=None, alpha=None,
                   fill=False, zorder=None):
        """Attach rendering style to this object. Returns self for chaining."""
        if not hasattr(self, '_render_style'):
            self._render_style = {}
        if color     is not None: self._render_style['color']     = color
        if linewidth is not None: self._render_style['linewidth'] = linewidth
        if alpha     is not None: self._render_style['alpha']     = alpha
        if zorder    is not None: self._render_style['zorder']    = zorder
        self._render_style['fill'] = bool(fill)
        return self

    def rescale(self, f):
        """
        Scale geometry by factor f (called after each box compression during swell).
        Base implementation is a no-op; subclasses override.
        """
        pass

    def resolved(self, t=0.0):
        """
        Return flat list of primitive dicts at time t.
        Subclasses override this.
        """
        raise NotImplementedError

    def to_make_prim_list(self, t=0.0):
        """
        Convert resolved(t) to list of 6-tuples for make_prim_data():
            (primitive_obj, k_pen, vel, r_c_wall, omega, r_ref)
        """
        out = []
        for d in self.resolved(t):
            out.append((d['prim'], d['k_pen'], d['vel'],
                        d['r_c_wall'], d['omega'], d['r_ref']))
        return out

    def is_time_varying(self):
        """
        Return True if this object's geometry changes with time.
        System._has_moving_objects() uses this to decide whether to rebuild
        prim_data every step. Override in subclasses whose resolved() is
        time-dependent but that don't use a MotionSpec (e.g. OscillatingWall).
        Default: True iff a non-static MotionSpec is attached.
        """
        ms = getattr(self, '_motion', None)
        return ms is not None and not ms.is_static()

    def region_polygon(self, t=0.0):
        """
        Return {'vertices': (N,2), 'exclusion': str} or None.
        Used by seeder for inside/outside point tests.
        Subclasses that have a closed boundary override this.
        """
        return None

    def contains(self, pt, t=0.0):
        """
        Test if point pt is in the *accessible* region at time t.
        - exclusion='exterior': accessible = inside polygon (seeder places particles inside)
        - exclusion='interior': accessible = outside polygon (seeder avoids interior)
        - None: no restriction
        """
        poly = self.region_polygon(t)
        if poly is None:
            return True
        inside = _point_in_polygon(pt, poly['vertices'])
        excl   = poly.get('exclusion') or self._exclusion
        if excl == 'exterior':
            return inside
        elif excl == 'interior':
            return not inside
        return True


# ── leaf primitives ───────────────────────────────────────────────────────────

class Wall(SimulationObject):
    """
    Single line-segment wall primitive.

    Parameters
    ----------
    p0, p1   : (2,) — endpoints in local frame (or world if not inside composite)
    normal   : (2,) — inward normal (pointing toward particles)
    """

    def __init__(self, p0, p1, normal=None):
        super().__init__('Wall')
        self._p0     = np.asarray(p0, dtype=float)
        self._p1     = np.asarray(p1, dtype=float)
        if normal is None:
            d    = self._p1 - self._p0
            d_len = np.linalg.norm(d)
            # Default inward normal: left of p0→p1
            self._normal = np.array([-d[1], d[0]]) / d_len
        else:
            n = np.asarray(normal, dtype=float)
            self._normal = n / np.linalg.norm(n)

    def resolved(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, _ = self._motion.displacement(t)
            p0 = self._p0 + np.array([dx, dy])
            p1 = self._p1 + np.array([dx, dy])
        else:
            p0, p1 = self._p0.copy(), self._p1.copy()
        return [{
            'prim'    : LineSegment(p0, p1, self._normal.copy()),
            'k_pen'   : self._k_pen,
            'vel'     : np.zeros(2),
            'omega'   : 0.0,
            'r_ref'   : np.zeros(2),
            'r_c_wall': self._r_c_wall,
            'exclusion': self._exclusion,
        }]


class OscillatingWall(Wall):
    """
    Horizontal wall that oscillates as y(t) = y0 + sign * A * sin(omega * t).

    Uses the same formula as the TF test reference scripts, guaranteeing
    bit-identical wall positions when cross-checking API vs TF direct.

    Parameters
    ----------
    y0       : float — rest position
    half_w   : float — half-width (wall spans [-half_w, half_w])
    A        : float — oscillation amplitude
    omega    : float — angular frequency (rad/s)
    is_top   : bool  — True = top wall (sign=-1, normal=[0,-1]);
                       False = bottom wall (sign=+1, normal=[0,+1])
    r_c_wall : float — wall capsule radius (default 0)
    """

    def __init__(self, y0, half_w, A, omega, is_top=True, r_c_wall=0.0):
        sign   = -1.0 if is_top else +1.0
        normal = [0, -1] if is_top else [0, +1]
        super().__init__([-half_w, float(y0)], [half_w, float(y0)], normal=normal)
        self.set_r_c_wall(r_c_wall)
        self._y0    = float(y0)
        self._A     = float(A)
        self._omega = float(omega)
        self._sign  = float(sign)

    def y_at(self, t):
        """Wall centre y position at time t."""
        return self._y0 + self._sign * self._A * np.sin(self._omega * float(t))

    def wall_strain_at(self, t, R0=1.0):
        """|y(t) - y0| / R0."""
        return abs(self.y_at(t) - self._y0) / float(R0)

    def is_time_varying(self):
        return True

    def to_make_prim_list(self, t=0.0):
        """
        Returns the wall at REST POSITION with oscillation encoded as parameters.
        primitive_forces_tf computes y(t) = y0 + sign*A*sin(omega*t) inside the
        TF graph — prim_data never needs rebuilding for moving walls.
        """
        p0_rest = np.array([self._p0[0], self._y0])
        p1_rest = np.array([self._p1[0], self._y0])
        seg = LineSegment(p0_rest, p1_rest, self._normal.copy())
        osc = (self._A, self._omega, self._sign)
        return [(seg, self._k_pen, np.zeros(2), self._r_c_wall, 0.0, np.zeros(2), osc)]

    def resolved(self, t=0.0):
        y  = self.y_at(t)
        p0 = self._p0.copy(); p0[1] = y
        p1 = self._p1.copy(); p1[1] = y
        return [{
            'prim'     : LineSegment(p0, p1, self._normal.copy()),
            'k_pen'    : self._k_pen,
            'vel'      : np.zeros(2),
            'omega'    : 0.0,
            'r_ref'    : np.zeros(2),
            'r_c_wall' : self._r_c_wall,
            'exclusion': self._exclusion,
        }]


class ArcWall(SimulationObject):
    """
    Single arc primitive (circle boundary).

    Parameters
    ----------
    center      : (2,) — arc center in local frame
    radius      : float
    convex      : bool — True = concave container (particles inside); False = convex obstacle
    angle_range : (float, float) or None — (theta_min, theta_max) in radians;
                  None means full circle.  Used for rendering only (physics sees full circle).
    """

    def __init__(self, center, radius, convex=True, angle_range=None):
        super().__init__('ArcWall')
        self._center      = np.asarray(center, dtype=float)
        self._radius      = float(radius)
        self._convex      = bool(convex)
        self._angle_range = angle_range   # (t0, t1) radians or None

    def resolved(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, _ = self._motion.displacement(t)
            center = self._center + np.array([dx, dy])
        else:
            center = self._center.copy()
        return [{
            'prim'    : Arc(center, self._radius,
                            angle_range=self._angle_range, convex=self._convex),
            'k_pen'   : self._k_pen,
            'vel'     : np.zeros(2),
            'omega'   : 0.0,
            'r_ref'   : np.zeros(2),
            'r_c_wall': self._r_c_wall,
            'exclusion': self._exclusion,
        }]


# ── composite object ──────────────────────────────────────────────────────────

class CompositeObject(SimulationObject):
    """
    Ordered list of child primitives (Wall or ArcWall) with a shared
    origin/rotation transform and optional motion.

    The world-frame position of a child with local position (cx, cy):
        effective_origin = origin + displacement(t)
        effective_theta  = theta  + dtheta(t)
        world = R(effective_theta) @ (cx, cy) + effective_origin
    """

    def __init__(self, kind='Composite'):
        super().__init__(kind)
        self._children      = []    # list of (child_obj, local_pos, local_theta)
        self._origin        = np.zeros(2)
        self._theta         = 0.0
        self._r_ref_offset  = np.zeros(2)   # body-frame offset from first prim p0 to r_ref

    def add_primitive(self, child_obj, local_pos=(0.0, 0.0), local_theta=0.0):
        """
        Append child to the composite.  local_pos and local_theta describe
        the child's position/orientation in the composite's local frame.
        """
        self._children.append((child_obj,
                                np.asarray(local_pos, dtype=float),
                                float(local_theta)))
        return self

    def set_origin(self, x, y):
        self._origin = np.array([float(x), float(y)])
        return self

    def set_rotation(self, theta):
        """Static orientation of the composite (radians)."""
        self._theta = float(theta)
        return self

    def set_r_ref_offset(self, dx, dy):
        """
        Body-frame offset from the first primitive's p0 to the rotation reference point.
        r_ref_abs = p0_first_world + R(theta) @ (dx, dy).
        Pre-built shapes set this automatically; CustomObject users call it manually.
        """
        self._r_ref_offset = np.array([float(dx), float(dy)])
        return self

    def to_make_prim_list(self, t=0.0):
        """
        Build TF prim tuples with correct motion parameters.
        Always uses rest positions (t=0 displacement) + encodes motion as parameters.
        """
        prims = self.resolved(0.0)

        if self._motion is not None:
            if self._motion.is_parametric():
                vx     = self._motion.vx_dc
                vy     = self._motion.vy_dc
                omega  = self._motion.omega_dc
            else:
                vx, vy, omega = self._motion.velocity(0.0)
        else:
            vx = vy = omega = 0.0

        # Compute r_ref_abs: origin (center) by default; or first-prim p0 + rotated offset
        if np.any(self._r_ref_offset != 0.0) and prims:
            p0_first = prims[0]['prim'].p0
            c, s = np.cos(self._theta), np.sin(self._theta)
            off = self._r_ref_offset
            r_ref_abs = p0_first + np.array([c * off[0] - s * off[1],
                                              s * off[0] + c * off[1]])
        else:
            r_ref_abs = self._origin.copy()

        vel = np.array([vx, vy])
        out = []
        for d in prims:
            out.append((d['prim'], d['k_pen'], vel, d['r_c_wall'], omega, r_ref_abs))
        return out

    def resolved(self, t=0.0):
        """
        Flatten all children into world-frame primitive dicts.
        Motion is applied to compute effective origin and theta at time t.
        """
        t = float(t)
        if self._motion is not None:
            dx, dy, dtheta = self._motion.displacement(t)
            vx, vy, om_parent = self._motion.velocity(t)
        else:
            dx, dy, dtheta = 0.0, 0.0, 0.0
            vx, vy, om_parent = 0.0, 0.0, 0.0

        ox    = self._origin[0] + dx
        oy    = self._origin[1] + dy
        theta = self._theta + dtheta

        c, s = np.cos(theta), np.sin(theta)
        r_ref = np.array([ox, oy])

        out = []
        for (child, local_pos, local_theta_child) in self._children:
            # World position of this child's local origin
            wx = ox + c * local_pos[0] - s * local_pos[1]
            wy = oy + s * local_pos[0] + c * local_pos[1]

            # Effective child rotation in world frame
            child_theta_world = theta + local_theta_child

            # Get child's own resolved list (in local frame) and transform to world
            child_prims = child.resolved(0.0)   # get in local frame (t=0)
            for d in child_prims:
                prim = d['prim']
                if isinstance(prim, LineSegment):
                    # Rotate and translate endpoints and normal
                    p0_w = _rotate_translate(prim.p0, child_theta_world, wx, wy)
                    p1_w = _rotate_translate(prim.p1, child_theta_world, wx, wy)
                    n_w  = _rotate(prim.normal, child_theta_world)
                    world_prim = LineSegment(p0_w, p1_w, n_w)
                elif isinstance(prim, Arc):
                    c_w = _rotate_translate(prim.center, child_theta_world, wx, wy)
                    world_prim = Arc(c_w, prim.radius, convex=prim.convex)
                else:
                    world_prim = prim   # unknown, pass through

                out.append({
                    'prim'    : world_prim,
                    'k_pen'   : d['k_pen'],
                    'vel'     : np.zeros(2),
                    'omega'   : 0.0,
                    'r_ref'   : np.zeros(2),
                    'r_c_wall': d['r_c_wall'],
                    'exclusion': d['exclusion'] or self._exclusion,
                })
        return out

    def region_polygon(self, t=0.0):
        """
        Build region polygon by extracting wall endpoints in order.
        Returns None if no closed polygon can be formed.
        Subclasses override for special shapes.
        """
        return None


# ── helpers ───────────────────────────────────────────────────────────────────

def _rotate(v, theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([c * v[0] - s * v[1], s * v[0] + c * v[1]])


def _rotate_translate(v, theta, ox, oy):
    w = _rotate(v, theta)
    return w + np.array([ox, oy])


# ── pre-built composite objects ───────────────────────────────────────────────

class Box(CompositeObject):
    """
    Fully closed rectangular box with 4 wall primitives and inward normals.

    Parameters
    ----------
    width, height : float — box dimensions
    x0, y0        : float — center position (default 0,0)
    theta         : float — rotation angle (default 0)
    exclusion     : 'exterior' = particles live inside (default)
    """

    def __init__(self, width, height, x0=0.0, y0=0.0, theta=0.0,
                 exclusion='exterior'):
        super().__init__('Box')
        self.set_origin(x0, y0)
        self.set_rotation(theta)
        self.set_exclusion(exclusion)
        self._width  = float(width)
        self._height = float(height)

        w2, h2 = width / 2.0, height / 2.0

        # 4 walls in local frame, normals pointing inward (+x, +y inside)
        # Bottom: p0=(-w2,-h2) p1=(+w2,-h2), normal=(0,+1)
        # Top:    p0=(+w2,+h2) p1=(-w2,+h2), normal=(0,-1)
        # Right:  p0=(+w2,-h2) p1=(+w2,+h2), normal=(-1, 0)
        # Left:   p0=(-w2,+h2) p1=(-w2,-h2), normal=(+1, 0)
        walls_local = [
            ((-w2, -h2), ( w2, -h2), ( 0,  1)),   # bottom
            (( w2,  h2), (-w2,  h2), ( 0, -1)),   # top
            (( w2, -h2), ( w2,  h2), (-1,  0)),   # right
            ((-w2,  h2), (-w2, -h2), ( 1,  0)),   # left
        ]
        for (p0, p1, n) in walls_local:
            child_wall = Wall(p0, p1, n)
            child_wall.set_exclusion(exclusion)
            self.add_primitive(child_wall)

    def region_polygon(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, dtheta = self._motion.displacement(t)
        else:
            dx, dy, dtheta = 0.0, 0.0, 0.0
        ox    = self._origin[0] + dx
        oy    = self._origin[1] + dy
        theta = self._theta + dtheta

        w2, h2 = self._width / 2.0, self._height / 2.0
        corners_local = [(-w2, -h2), (w2, -h2), (w2, h2), (-w2, h2)]
        verts = np.array([_rotate_translate(np.array(c), theta, ox, oy)
                          for c in corners_local])
        return {'vertices': verts, 'exclusion': self._exclusion}


class Channel(CompositeObject):
    """
    Two horizontal walls (top and bottom). Periodic in x assumed by user.

    Parameters
    ----------
    width, height : float — channel dimensions (width determines wall length)
    x0, y0        : float — center position
    exclusion     : 'exterior' = particles between the walls
    """

    def __init__(self, width, height, x0=0.0, y0=0.0, theta=0.0,
                 exclusion='exterior'):
        super().__init__('Channel')
        self.set_origin(x0, y0)
        self.set_rotation(theta)
        self.set_exclusion(exclusion)
        self._width  = float(width)
        self._height = float(height)

        w2, h2 = width / 2.0, height / 2.0
        # Bottom wall: normal upward
        bot = Wall((-w2, -h2), (w2, -h2), (0, 1))
        bot.set_exclusion(exclusion)
        self.add_primitive(bot)
        # Top wall: normal downward
        top = Wall((w2, h2), (-w2, h2), (0, -1))
        top.set_exclusion(exclusion)
        self.add_primitive(top)

    def region_polygon(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, dtheta = self._motion.displacement(t)
        else:
            dx, dy, dtheta = 0.0, 0.0, 0.0
        ox    = self._origin[0] + dx
        oy    = self._origin[1] + dy
        theta = self._theta + dtheta
        w2, h2 = self._width / 2.0, self._height / 2.0
        corners = [(-w2, -h2), (w2, -h2), (w2, h2), (-w2, h2)]
        verts   = np.array([_rotate_translate(np.array(c), theta, ox, oy)
                            for c in corners])
        return {'vertices': verts, 'exclusion': self._exclusion}

    def rescale(self, f):
        self._width  *= float(f)
        self._height *= float(f)
        self._origin  = self._origin * float(f)
        if self._motion is not None:
            self._motion.rescale(f)
        self._primitives = []
        w2, h2 = self._width / 2.0, self._height / 2.0
        bot = Wall((-w2, -h2), (w2, -h2), (0, 1))
        bot.set_exclusion(self._exclusion)
        self.add_primitive(bot)
        top = Wall((w2, h2), (-w2, h2), (0, -1))
        top.set_exclusion(self._exclusion)
        self.add_primitive(top)


class CouetteCell(CompositeObject):
    """
    Annular cell: outer arc (container) + inner arc (cylindrical obstacle).

    Particles live in the annulus (inner_radius < r < outer_radius).

    Parameters
    ----------
    inner_radius, outer_radius : float
    x0, y0                     : float — center
    exclusion                  : 'exterior' (particles between arcs; default)
    """

    def __init__(self, inner_radius, outer_radius, x0=0.0, y0=0.0,
                 exclusion='exterior'):
        super().__init__('CouetteCell')
        self.set_origin(x0, y0)
        self.set_exclusion(exclusion)
        self._inner_r = float(inner_radius)
        self._outer_r = float(outer_radius)

        # Outer arc: convex=False = concave container; particles inside, pushed inward
        outer = ArcWall((0, 0), outer_radius, convex=False)
        outer.set_exclusion('exterior')
        self.add_primitive(outer)

        # Inner arc: convex=True = convex obstacle; particles outside, pushed outward
        inner = ArcWall((0, 0), inner_radius, convex=True)
        inner.set_exclusion('interior')
        self.add_primitive(inner)

    def region_polygon(self, t=0.0):
        """Annulus approximated as two polygons (outer - inner), encoded as outer."""
        t = float(t)
        if self._motion is not None:
            dx, dy, _ = self._motion.displacement(t)
        else:
            dx, dy = 0.0, 0.0
        cx, cy = self._origin[0] + dx, self._origin[1] + dy

        n_pts = 64
        angles = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
        # Return outer polygon; seeder will also exclude inner disk separately
        verts = np.column_stack([
            cx + self._outer_r * np.cos(angles),
            cy + self._outer_r * np.sin(angles),
        ])
        return {'vertices': verts, 'exclusion': 'exterior',
                'inner_radius': self._inner_r,
                'center': np.array([cx, cy])}

    def contains(self, pt, t=0.0):
        """True if pt is in the annulus at time t."""
        t = float(t)
        if self._motion is not None:
            dx, dy, _ = self._motion.displacement(t)
        else:
            dx, dy = 0.0, 0.0
        cx, cy = self._origin[0] + dx, self._origin[1] + dy
        r = np.sqrt((float(pt[0]) - cx)**2 + (float(pt[1]) - cy)**2)
        return self._inner_r < r < self._outer_r

    def exclusion_area(self, t=0.0):
        """Area excluded from the enclosing (2*R_outer)² bounding box.

        = corner area (outside outer circle) + inner disk area
        = 4*R_outer² − π*(R_outer² − R_inner²)

        Assumes the System box is exactly 2*R_outer × 2*R_outer, so that
        System._accessible_area() = π*(R_outer² − R_inner²).
        """
        return 4.0 * self._outer_r**2 - np.pi * (self._outer_r**2 - self._inner_r**2)

    def rescale(self, f):
        """Scale radii and origin by f (called at each swell compression step)."""
        f = float(f)
        self._inner_r *= f
        self._outer_r *= f
        self._origin   = self._origin * f
        # Rebuild arc children with updated radii
        self._children = []
        outer = ArcWall((0, 0), self._outer_r, convex=False)
        outer.set_exclusion('exterior')
        self.add_primitive(outer)
        inner = ArcWall((0, 0), self._inner_r, convex=True)
        inner.set_exclusion('interior')
        self.add_primitive(inner)


class CircleObstacle(CompositeObject):
    """
    Circular obstacle — particles excluded from interior.

    Parameters
    ----------
    radius : float
    x0, y0 : float — center
    """

    def __init__(self, radius, x0=0.0, y0=0.0, exclusion='interior'):
        super().__init__('CircleObstacle')
        self.set_origin(x0, y0)
        self.set_exclusion(exclusion)
        self._radius = float(radius)
        arc = ArcWall((0, 0), radius, convex=False)
        arc.set_exclusion(exclusion)
        self.add_primitive(arc)

    def region_polygon(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, _ = self._motion.displacement(t)
        else:
            dx, dy = 0.0, 0.0
        cx, cy = self._origin[0] + dx, self._origin[1] + dy
        n_pts  = 64
        angles = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
        verts  = np.column_stack([
            cx + self._radius * np.cos(angles),
            cy + self._radius * np.sin(angles),
        ])
        return {'vertices': verts, 'exclusion': self._exclusion}


class RegularPolygon(CompositeObject):
    """
    Regular N-sided polygon — N wall primitives with inward normals.

    Parameters
    ----------
    sides  : int   — number of sides
    radius : float — circumradius
    x0, y0 : float — center
    theta  : float — rotation of first vertex
    exclusion : 'interior' = particles excluded from inside (obstacle)
                'exterior' = particles inside polygon (container)
    """

    def __init__(self, sides, radius, x0=0.0, y0=0.0, theta=0.0,
                 exclusion='interior'):
        super().__init__('RegularPolygon')
        self.set_origin(x0, y0)
        self.set_rotation(theta)
        self.set_exclusion(exclusion)
        self._sides  = int(sides)
        self._radius = float(radius)

        angles = [2 * np.pi * k / sides for k in range(sides)]
        verts  = [(radius * np.cos(a), radius * np.sin(a)) for a in angles]

        for k in range(sides):
            p0 = verts[k]
            p1 = verts[(k + 1) % sides]
            # Inward normal: perpendicular pointing toward centre
            d  = np.array(p1) - np.array(p0)
            n  = np.array([-d[1], d[0]])
            n /= np.linalg.norm(n)
            # Ensure inward: midpoint of edge → centre direction
            mid = 0.5 * (np.array(p0) + np.array(p1))
            if np.dot(n, -mid) < 0:
                n = -n
            wall = Wall(p0, p1, n)
            wall.set_exclusion(exclusion)
            self.add_primitive(wall)

    def region_polygon(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, dtheta = self._motion.displacement(t)
        else:
            dx, dy, dtheta = 0.0, 0.0, 0.0
        ox    = self._origin[0] + dx
        oy    = self._origin[1] + dy
        theta_eff = self._theta + dtheta
        angles = [2 * np.pi * k / self._sides for k in range(self._sides)]
        local  = [(self._radius * np.cos(a), self._radius * np.sin(a)) for a in angles]
        verts  = np.array([_rotate_translate(np.array(c), theta_eff, ox, oy)
                           for c in local])
        return {'vertices': verts, 'exclusion': self._exclusion}


class SquareObstacle(CompositeObject):
    """
    Square obstacle — particles excluded from the interior.

    The obstacle is placed by specifying its center (x0, y0).
    The r_ref for rigid-body spin is set to the center automatically.

    Parameters
    ----------
    side      : float — side length
    x0, y0   : float — center position (default 0, 0)
    theta     : float — initial rotation (radians, default 0)
    exclusion : 'interior' (default) — seeder keeps particles outside
    """

    def __init__(self, side, x0=0.0, y0=0.0, theta=0.0, exclusion='interior'):
        super().__init__('SquareObstacle')
        self._side = float(side)
        self.set_origin(x0, y0)
        self.set_rotation(theta)
        self.set_exclusion(exclusion)

        w2 = self._side / 2.0
        # r_ref_offset: from first corner (-w2,-w2) to center in body frame = (w2, w2)
        self.set_r_ref_offset(w2, w2)

        # Normals point OUTWARD (toward particles outside the obstacle)
        walls_local = [
            ((-w2, -w2), ( w2, -w2), ( 0, -1)),   # bottom  → normal down
            (( w2,  w2), (-w2,  w2), ( 0,  1)),   # top     → normal up
            (( w2, -w2), ( w2,  w2), ( 1,  0)),   # right   → normal right
            ((-w2,  w2), (-w2, -w2), (-1,  0)),   # left    → normal left
        ]
        for (p0, p1, n) in walls_local:
            child_wall = Wall(p0, p1, n)
            child_wall.set_exclusion(exclusion)
            self.add_primitive(child_wall)

    def region_polygon(self, t=0.0):
        t = float(t)
        if self._motion is not None:
            dx, dy, dtheta = self._motion.displacement(t)
        else:
            dx, dy, dtheta = 0.0, 0.0, 0.0
        ox    = self._origin[0] + dx
        oy    = self._origin[1] + dy
        theta = self._theta + dtheta
        w2    = self._side / 2.0
        corners_local = [(-w2, -w2), (w2, -w2), (w2, w2), (-w2, w2)]
        verts = np.array([_rotate_translate(np.array(c), theta, ox, oy)
                          for c in corners_local])
        return {'vertices': verts, 'exclusion': self._exclusion}

    def exclusion_area(self, t=0.0):
        """Interior area excluded from the accessible domain."""
        return self._side ** 2

    def rescale(self, f):
        self._side   *= float(f)
        self._origin  = self._origin * float(f)
        if self._motion is not None:
            self._motion.rescale(f)
        w2 = self._side / 2.0
        self.set_r_ref_offset(w2, w2)
        self._primitives = []
        walls_local = [
            ((-w2, -w2), ( w2, -w2), ( 0, -1)),
            (( w2,  w2), (-w2,  w2), ( 0,  1)),
            (( w2, -w2), ( w2,  w2), ( 1,  0)),
            ((-w2,  w2), (-w2, -w2), (-1,  0)),
        ]
        for (p0, p1, n) in walls_local:
            child_wall = Wall(p0, p1, n)
            child_wall.set_exclusion(self._exclusion)
            self.add_primitive(child_wall)

    def to_make_prim_list(self, t=0.0):
        """
        Submit walls as a single Polygon primitive so that make_prim_data assigns
        them a shared group ID.  Only the closest wall is active for any given
        particle node, which prevents forces from walls on the 'wrong' side.

        Vertices are listed in CW order so that Polygon's left-hand normals point
        OUTWARD — the correct direction for an obstacle (particles outside).
        """
        from src.simulation.contact_primitives import Polygon as CPoly

        poly  = self.region_polygon(0.0)           # CCW world-frame corners
        verts_cw = poly['vertices'][::-1]           # reverse → CW

        prim_obj = CPoly(verts_cw)

        if self._motion is not None:
            if self._motion.is_parametric():
                vx     = self._motion.vx_dc
                vy     = self._motion.vy_dc
                omega  = self._motion.omega_dc
            else:
                vx, vy, omega = self._motion.velocity(0.0)
        else:
            vx = vy = omega = 0.0

        vel     = np.array([vx, vy])
        r_ref   = self._origin.copy()   # center of the obstacle
        return [(prim_obj, self._k_pen, vel, self._r_c_wall, omega, r_ref)]


class HopperRegion(SimulationObject):
    """
    2-D hopper: funnel + vertical reservoir, open at top and bottom.

    The accessible region is INSIDE the hopper polygon (exclusion='exterior').
    Four wall primitives are created (left funnel, right funnel, left vertical,
    right vertical); there is no floor or ceiling.

    Parameters
    ----------
    x_c       : float — center x of outlet (and hopper symmetry axis)
    W_out     : float — outlet width (at y = y_bot)
    W_res     : float — reservoir width (full width of straight section)
    theta_deg : float — funnel wall angle from horizontal in degrees (default 30)
    h_res     : float — reservoir height above the funnel top
    y_bot     : float — y-coordinate of the outlet (default 0)
    """

    def __init__(self, x_c, W_out, W_res, theta_deg=30.0, h_res=42.0, y_bot=0.0,
                 corner_radius=0.4):
        super().__init__('HopperRegion')
        self.set_exclusion('exterior')
        self._x_c           = float(x_c)
        self._W_out         = float(W_out)
        self._W_res         = float(W_res)
        self._theta_deg     = float(theta_deg)
        self._h_res         = float(h_res)
        self._y_bot         = float(y_bot)
        self._corner_radius = float(corner_radius)
        self._rebuild()

    def _rebuild(self):
        import math
        xc = self._x_c
        Wo = self._W_out
        Wr = self._W_res
        th = math.radians(self._theta_deg)
        hf = (Wr - Wo) / 2.0 * math.tan(th)
        self._h_funnel = hf
        yb = self._y_bot
        hr = self._h_res
        r  = self._corner_radius

        h_exit = max(Wo * 0.5, 3.0)

        # Named corners of the hopper polygon
        BL     = np.array([xc - Wo/2,  yb])
        BR     = np.array([xc + Wo/2,  yb])
        BL_ext = np.array([xc - Wo/2,  yb - h_exit])
        BR_ext = np.array([xc + Wo/2,  yb - h_exit])
        ULF    = np.array([xc - Wr/2,  yb + hf])
        ULL    = np.array([xc - Wr/2,  yb + hf + hr])
        ULR    = np.array([xc + Wr/2,  yb + hf + hr])
        URF    = np.array([xc + Wr/2,  yb + hf])

        # 6-vertex polygon for exclusion/seeding (CCW, original outlet corners)
        self._verts = np.array([BL, ULF, ULL, ULR, URF, BR])

        # Center point for inward-normal computation
        center = np.array([xc, yb + hf + hr / 2.0])

        def _inward_wall(p0, p1):
            d     = np.asarray(p1) - np.asarray(p0)
            d_len = np.linalg.norm(d)
            n     = np.array([-d[1], d[0]]) / d_len
            mid   = 0.5 * (np.asarray(p0) + np.asarray(p1))
            if np.dot(n, center - mid) < 0:
                n = -n
            return Wall(np.asarray(p0, dtype=float), np.asarray(p1, dtype=float), n)

        # ── Arc corner rounding ───────────────────────────────────────────────
        # At each outlet corner the funnel wall (angled at theta from horizontal)
        # meets the vertical guide wall.  The sharp corner endpoint is replaced by
        # a circular arc of radius r whose center sits in the wall material at
        # distance r from both walls.
        #
        # Geometry (left corner, right is symmetric):
        #   C_L  = (xc - Wo/2 - r,  yb - r*(1-sin(th))/cos(th))
        #   T1_L = tangent on funnel wall = C_L + r*(sin(th),  cos(th))
        #   T2_L = tangent on guide wall  = (xc - Wo/2, C_L_y)
        #   arc angle range: 0  → pi/2 - th  (counterclockwise, T2 → T1)
        f_sin   = math.sin(th)
        f_cos   = math.cos(th)
        dy_arc  = r * (1.0 - f_sin) / f_cos

        C_L  = np.array([xc - Wo/2 - r,  yb - dy_arc])
        T1_L = C_L + r * np.array([ f_sin,  f_cos])   # tangent on left funnel wall
        T2_L = np.array([xc - Wo/2,  C_L[1]])          # tangent on left guide wall

        C_R  = np.array([xc + Wo/2 + r,  yb - dy_arc])
        T1_R = C_R + r * np.array([-f_sin,  f_cos])   # tangent on right funnel wall
        T2_R = np.array([xc + Wo/2,  C_R[1]])          # tangent on right guide wall

        ang_lo_L = 0.0
        ang_hi_L = math.pi / 2.0 - th    # e.g. 60° for theta=30°
        ang_lo_R = math.pi / 2.0 + th    # e.g. 120° for theta=30°
        ang_hi_R = math.pi

        self._walls = [
            _inward_wall(T1_L,   ULF),      # left funnel wall, ends at arc tangent
            _inward_wall(T1_R,   URF),      # right funnel wall
            _inward_wall(ULF,    ULL),      # left vertical reservoir wall
            _inward_wall(URF,    ULR),      # right vertical reservoir wall
            _inward_wall(BL_ext, T2_L),    # left outlet guide (below arc tangent)
            _inward_wall(BR_ext, T2_R),    # right outlet guide
        ]

        aw_L = ArcWall(C_L, r, convex=True, angle_range=(ang_lo_L, ang_hi_L))
        aw_L.set_exclusion('interior')
        aw_R = ArcWall(C_R, r, convex=True, angle_range=(ang_lo_R, ang_hi_R))
        aw_R.set_exclusion('interior')
        self._arc_walls = [aw_L, aw_R]

    def region_polygon(self, t=0.0):
        return {'vertices': self._verts, 'exclusion': 'exterior'}

    def exclusion_area(self, t=0.0):
        return 0.0

    def polygon_interior_area(self):
        """Shoelace area of the accessible hopper polygon."""
        v = self._verts
        x, y = v[:, 0], v[:, 1]
        return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))

    def resolved(self, t=0.0):
        out = []
        for w in self._walls:
            out.extend(w.resolved(t))
        for aw in self._arc_walls:
            out.extend(aw.resolved(t))
        return out

    def rescale(self, f):
        self._x_c           *= float(f)
        self._W_out         *= float(f)
        self._W_res         *= float(f)
        self._h_res         *= float(f)
        self._y_bot         *= float(f)
        self._corner_radius *= float(f)
        self._rebuild()

    @property
    def h_funnel(self):
        return self._h_funnel

    @property
    def h_total(self):
        return self._h_funnel + self._h_res


class CustomObject(CompositeObject):
    """
    User-assembled composite. Add primitives manually via add_primitive().
    """

    def __init__(self, x0=0.0, y0=0.0, theta=0.0, exclusion=None):
        super().__init__('Custom')
        self.set_origin(x0, y0)
        self.set_rotation(theta)
        self.set_exclusion(exclusion)

    def region_polygon(self, t=0.0):
        """
        Build region polygon from the p0 endpoints of assembled Wall primitives.
        Assumes walls were added in a consistent order forming a closed polygon.
        Returns None if fewer than 3 wall primitives exist.
        """
        prims = self.resolved(t)
        verts = []
        for d in prims:
            p = d['prim']
            if isinstance(p, LineSegment):
                verts.append(p.p0.copy())
        if len(verts) < 3:
            return None
        return {'vertices': np.array(verts), 'exclusion': self._exclusion}
