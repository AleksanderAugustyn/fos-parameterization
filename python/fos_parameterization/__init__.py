"""Python bindings for the Fourier-over-Spheroid (FoS) nuclear-shape parameterization.

Two tiers: one-shot module-level functions that build and discard their own
tables per call, and :class:`Cache`, a read-only cache built once and shared
freely between threads, for many shapes at one resolution.
Validation failures are returned as a :class:`Status` on the result object;
only usage errors raise :class:`FosParamError`.
"""
from ._cdefs import MAX_PARAMS, N_POINTS_FLOOR
from ._libloader import load_library
from .api import (
    Cache,
    FosParamError,
    FosShape,
    NeckResult,
    RadiusDerivativeResult,
    RadiusGridResult,
    RhoZGridResult,
    StarConvexityResult,
    Status,
    a2,
    neck,
    radius_and_derivative,
    radius_grid,
    rho_at_z,
    rho_z_grid,
    rho_z_grid_unchecked,
    shape,
    star_convexity_optimum,
    status_message,
    theta_grid,
    z_shift,
)

__all__ = [
    "Cache", "FosParamError", "FosShape", "NeckResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoZGridResult",
    "StarConvexityResult", "Status",
    "a2", "neck", "radius_and_derivative", "radius_grid", "rho_at_z",
    "rho_z_grid", "rho_z_grid_unchecked", "shape", "star_convexity_optimum",
    "status_message", "theta_grid", "z_shift",
    "MAX_PARAMS", "N_POINTS_FLOOR",
    "load_library",
]
