"""High-level API: Status, result dataclasses, standalone functions.

Hot-path functions return a Status instead of raising, so parameter sweeps
and plotters can render invalid shapes. Only the library loader raises.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional

import numpy as np
import numpy.typing as npt

from ._cdefs import MESSAGE_BUFFER_SIZE, c_dbl_p, configure
from ._libloader import load_library

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class Status(IntEnum):
    """Mirrors FOS_* codes in include/fos_parameterization.h."""
    VALID = 0
    ERROR_RHO_NEGATIVE = 1
    ERROR_NOT_STAR_CONVEX = 2
    ERROR_INVALID_C = 3
    ERROR_BEAK_SINGULARITY = 4
    ERROR_INVALID_ARGUMENTS = 5


@dataclass(frozen=True)
class RadiusGridResult:
    """One R(θ) evaluation. z_shift is the total shift baked into the grid."""
    radii: npt.NDArray[np.float64]
    z_shift: float
    status: Status
    message: str

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


@dataclass(frozen=True)
class RhoProfileResult:
    """Cylindrical profile ρ(z) with exact dρ/dz, COM frame."""
    z: npt.NDArray[np.float64]
    rho: npt.NDArray[np.float64]
    drho_dz: npt.NDArray[np.float64]
    status: Status
    message: str

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


@dataclass(frozen=True)
class NeckResult:
    """Neck position/radius in the COM frame; found=False for neckless shapes."""
    z_neck: float
    rho_neck: float
    found: bool
    status: Status

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


@dataclass(frozen=True)
class ShapeResult:
    """Shape validity, total z-shift, and analytic pole radii R(0), R(π)."""
    z_shift: float
    r_north: float
    r_south: float
    status: Status
    message: str

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


@dataclass(frozen=True)
class RadiusDerivativeResult:
    """Batch R(θ) and analytic dR/dθ at caller-supplied θ nodes."""
    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status

    @property
    def ok(self) -> bool:
        return self.status == Status.VALID


def theta_grid(n_grid: int) -> npt.NDArray[np.float64]:
    """The θ grid the library evaluates on: n_grid uniform points over [0, π]."""
    return np.linspace(0.0, np.pi, n_grid)


def _as_params(params: npt.ArrayLike) -> npt.NDArray[np.float64]:
    arr = np.ascontiguousarray(params, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"params must be 1-D, got shape {arr.shape}")
    return arr


def radius_grid(params: npt.ArrayLike, n_grid: int) -> RadiusGridResult:
    """Validate the shape and compute R(θ) on n_grid uniform points over [0, π]."""
    arr = _as_params(params)
    radii = np.zeros(int(n_grid), dtype=np.float64)
    z_shift = ctypes.c_double(0.0)
    buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
    status = _get_lib().fos_compute_radius_grid(
        arr.ctypes.data_as(c_dbl_p), arr.size, int(n_grid),
        radii.ctypes.data_as(c_dbl_p), ctypes.byref(z_shift),
        MESSAGE_BUFFER_SIZE, buf)
    return RadiusGridResult(
        radii=radii, z_shift=z_shift.value, status=Status(status),
        message=buf.value.decode(errors="replace"))


def rho_profile(params: npt.ArrayLike, n_z: int) -> RhoProfileResult:
    """ρ(z) profile with exact dρ/dz in the COM frame (no star-convexity gate)."""
    arr = _as_params(params)
    n = int(n_z)
    z = np.zeros(n, dtype=np.float64)
    rho = np.zeros(n, dtype=np.float64)
    drho_dz = np.zeros(n, dtype=np.float64)
    buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
    status = _get_lib().fos_compute_rho_profile(
        arr.ctypes.data_as(c_dbl_p), arr.size, n,
        z.ctypes.data_as(c_dbl_p), rho.ctypes.data_as(c_dbl_p),
        drho_dz.ctypes.data_as(c_dbl_p),
        MESSAGE_BUFFER_SIZE, buf)
    return RhoProfileResult(
        z=z, rho=rho, drho_dz=drho_dz, status=Status(status),
        message=buf.value.decode(errors="replace"))


def neck(params: npt.ArrayLike) -> NeckResult:
    """Neck (interior ρ minimum), Newton-refined, COM frame."""
    arr = _as_params(params)
    z_neck = ctypes.c_double(0.0)
    rho_neck = ctypes.c_double(0.0)
    found = ctypes.c_int(0)
    status = _get_lib().fos_compute_neck(
        arr.ctypes.data_as(c_dbl_p), arr.size,
        ctypes.byref(z_neck), ctypes.byref(rho_neck), ctypes.byref(found))
    return NeckResult(
        z_neck=z_neck.value, rho_neck=rho_neck.value,
        found=bool(found.value), status=Status(status))


def shape(params: npt.ArrayLike, n_rho_grid: int) -> ShapeResult:
    """Resolve shape validity, total z-shift, and analytic pole radii.

    Parameters
    ----------
    params : array_like
        FoS parameters: params[0] = c, params[k-2] = a_k for k >= 3.
    n_rho_grid : int
        Internal rho(z) grid size (used verbatim; >= 2).
    """
    arr = _as_params(params)
    z_shift = ctypes.c_double(0.0)
    r_north = ctypes.c_double(0.0)
    r_south = ctypes.c_double(0.0)
    buf = ctypes.create_string_buffer(MESSAGE_BUFFER_SIZE)
    status = _get_lib().fos_compute_shape(
        arr.ctypes.data_as(c_dbl_p), arr.size, int(n_rho_grid),
        ctypes.byref(z_shift), ctypes.byref(r_north), ctypes.byref(r_south),
        MESSAGE_BUFFER_SIZE, buf)
    return ShapeResult(
        z_shift=z_shift.value, r_north=r_north.value, r_south=r_south.value,
        status=Status(status), message=buf.value.decode(errors="replace"))


def radius_and_derivative(params: npt.ArrayLike, thetas: npt.ArrayLike,
                          z_shift: float) -> RadiusDerivativeResult:
    """Batch R(theta) and analytic dR/dtheta at caller-supplied thetas.

    z_shift must come from shape(); degenerate params yield the unit-sphere
    fallback (r = 1, dr_dtheta = 0) — a library guarantee, not an error.
    """
    arr = _as_params(params)
    t = np.ascontiguousarray(thetas, dtype=np.float64)
    radii = np.zeros(t.size, dtype=np.float64)
    dr = np.zeros(t.size, dtype=np.float64)
    status = _get_lib().fos_compute_radius_and_derivative_at_thetas(
        arr.ctypes.data_as(c_dbl_p), arr.size,
        t.ctypes.data_as(c_dbl_p), t.size, float(z_shift),
        radii.ctypes.data_as(c_dbl_p), dr.ctypes.data_as(c_dbl_p))
    return RadiusDerivativeResult(radii=radii, dr_dtheta=dr, status=Status(status))


def z_shift(params: npt.ArrayLike) -> float:
    """Intrinsic COM z-shift (closed form). 0.0 for invalid params."""
    arr = _as_params(params)
    return float(_get_lib().fos_z_shift(arr.ctypes.data_as(c_dbl_p), arr.size))


def a2(params: npt.ArrayLike) -> float:
    """a2 from the volume-conservation constraint. 0.0 for empty params."""
    arr = _as_params(params)
    return float(_get_lib().fos_a2(arr.ctypes.data_as(c_dbl_p), arr.size))
