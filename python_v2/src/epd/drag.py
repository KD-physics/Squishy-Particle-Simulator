"""
drag.py — Stokes drag convenience helpers.

The primary parameter is the Ohnesorge number Oh = xi / v_ref, where
v_ref is the particle's intrinsic velocity scale:
  emulsion : v_ref = sqrt(gamma / (rho_d * R0))
  elastic  : v_ref = sqrt(El_t  / (rho_d * R0))

At the working point (gamma=rho_d=R0=1 for emulsion, rho_d=R0=1 for elastic):
  emulsion : xi = Oh
  elastic  : xi = Oh * sqrt(El_t)

Terminal velocity under gravity g (2D disk, radius R0, density rho_d):
  v_t = rho_d * R0 * g / (2 * xi)
      = g / (2 * Oh)    [emulsion working point]
"""

import numpy as np


def oh_from_terminal_velocity(v_t, g, particle_type='emulsion',
                               gamma=1.0, El_t=None, rho_d=1.0, R0=1.0):
    """
    Compute the Ohnesorge number Oh that gives terminal velocity v_t
    under gravity g for a single particle.

    Parameters
    ----------
    v_t           : float — target terminal velocity (same units as simulation)
    g             : float — gravitational acceleration (must be > 0)
    particle_type : str   — 'emulsion' or 'elastic'
    gamma         : float — surface tension (emulsion only; default 1.0)
    El_t          : float — membrane stiffness (elastic only; required if elastic)
    rho_d         : float — particle density (default 1.0)
    R0            : float — particle radius (default 1.0)

    Returns
    -------
    Oh : float — Ohnesorge number

    Notes
    -----
    Force balance: rho_d * pi * R0^2 * g = xi * 2*pi*R0 * v_t
    =>  xi = rho_d * R0 * g / (2 * v_t)
    =>  Oh = xi / v_ref
    """
    if g <= 0.0:
        raise ValueError("g must be > 0 to define terminal velocity.")
    if v_t <= 0.0:
        raise ValueError("v_t must be > 0.")

    xi = float(rho_d * R0 * g / (2.0 * v_t))

    if particle_type == 'emulsion':
        v_ref = float(np.sqrt(gamma / (rho_d * R0)))
    elif particle_type in ('elastic', 'rigid'):
        if El_t is None:
            raise ValueError("El_t must be provided for elastic particles.")
        v_ref = float(np.sqrt(El_t / (rho_d * R0)))
        if v_ref == 0.0:
            raise ValueError("El_t=0 — elastic reference velocity undefined.")
    else:
        raise ValueError(f"Unknown particle_type '{particle_type}'. "
                         f"Use 'emulsion' or 'elastic'.")

    return xi / v_ref


def terminal_velocity(Oh, g, particle_type='emulsion',
                      gamma=1.0, El_t=None, rho_d=1.0, R0=1.0):
    """
    Compute terminal velocity for given Oh and gravity g.
    Inverse of oh_from_terminal_velocity.
    """
    if g <= 0.0:
        raise ValueError("g must be > 0.")
    if Oh <= 0.0:
        raise ValueError("Oh must be > 0.")

    if particle_type == 'emulsion':
        v_ref = float(np.sqrt(gamma / (rho_d * R0)))
    elif particle_type in ('elastic', 'rigid'):
        if El_t is None:
            raise ValueError("El_t required for elastic particles.")
        v_ref = float(np.sqrt(El_t / (rho_d * R0)))
    else:
        raise ValueError(f"Unknown particle_type '{particle_type}'.")

    xi = Oh * v_ref
    return float(rho_d * R0 * g / (2.0 * xi))
