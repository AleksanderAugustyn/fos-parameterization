"""High-level API: Status, result objects, Cache, flat tier-1 functions.

Two tiers, matching the C API. The module-level functions build, use and
discard their own tables per call; :class:`Cache` keeps the per-shape working
state, so a sweep over many shapes only recomputes what the changed parameters
invalidated.

A cache is THREAD-CONFINED: every compute call mutates it, so give each thread
its own.

Shape-validation failures come back as result objects carrying a
:class:`Status` with every numeric output zero-filled; :class:`FosParamError`
is raised only for usage errors — a failed create, a closed handle, or non-1-D
input.

Frames: the R(theta) results and :func:`shape` /
:func:`star_convexity_optimum` report the TOTAL z-shift (intrinsic
centre-of-mass shift plus any extra star-convexity shift). The rho(z) grid and
the neck report in the COM frame, which is also what :func:`z_shift` returns.

Precondition — every ``params`` and ``thetas`` array must be finite. Non-finite
input is undefined behavior: the library cannot detect NaN under fast-math, so
a call returns ``Status.valid`` with NaN outputs instead of an error. Screen
inputs before calling.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass
from enum import IntEnum
from typing import Optional, Tuple

import numpy as np
import numpy.typing as npt

from ._cdefs import (CACHE_MAX_PARAMS, MAX_PARAMS, N_POINTS_FLOOR, c_dbl_p,
                     configure)
from ._libloader import load_library

_lib: Optional[ctypes.CDLL] = None


def _get_lib() -> ctypes.CDLL:
    """Load and configure the shared library once per process."""
    global _lib
    if _lib is None:
        _lib = configure(load_library())
    return _lib


class FosParamError(RuntimeError):
    """Raised for usage errors: failed create, closed handle, non-1-D input."""


class Status(IntEnum):
    """Status codes of the C API (``FOS_*`` in the header).

    0-6 are the shared shape-parameterization contract codes, identical in
    every library of the family; 100+ are FoS's own, append-only after 2.0.0.
    """

    valid = 0
    too_many_params = 1
    cache_not_initialized = 2
    invalid_grid = 3
    wrong_param_count = 4
    invalid_init = 5
    tables_not_initialized = 6
    rho_negative = 100
    not_star_convex = 101
    invalid_c = 102
    beak_singularity = 103
    convergence = 104
    buffer_mismatch = 105


def status_message(status: int) -> str:
    """Fixed description of a status code.

    Parameters
    ----------
    status : int
        A :class:`Status` value or its integer code.

    Returns
    -------
    str
        The library's static message; ``"unknown status code"`` for codes it
        does not know.
    """
    return _get_lib().fos_param_status_message(int(status)).decode()


# --- result objects --------------------------------------------------------

@dataclass(frozen=True)
class RadiusGridResult:
    """R(theta) on a theta set. Radii are zero-filled when ``status`` is nonzero."""

    radii: npt.NDArray[np.float64]
    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class RadiusDerivativeResult:
    """R(theta) and dR/dtheta. Both buffers are zero-filled on failure."""

    radii: npt.NDArray[np.float64]
    dr_dtheta: npt.NDArray[np.float64]
    status: Status

    @property
    def ok(self) -> bool:
        """True when the evaluation succeeded."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class FosShape:
    """A resolved shape: the TOTAL z-shift and the analytic pole radii.

    ``z_shift`` is the intrinsic COM shift plus any extra shift the
    star-convexity search applied, i.e. the frame the R(theta) outputs live in.
    """

    z_shift: float
    r_north: float
    r_south: float
    status: Status

    @property
    def ok(self) -> bool:
        """True when the shape resolved successfully."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class RhoZGridResult:
    """Cylindrical profile rho(z) with exact drho/dz, in the COM frame.

    ``z_shift`` is the intrinsic COM shift only — never the star-convexity
    shift. rho is exactly 0 at both tips.
    """

    z: npt.NDArray[np.float64]
    rho: npt.NDArray[np.float64]
    drho_dz: npt.NDArray[np.float64]
    z_shift: float
    status: Status

    @property
    def ok(self) -> bool:
        """True when the profile was computed."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class NeckResult:
    """Interior minimum of rho between two maxima, in the COM frame.

    A valid shape without a neck is ``ok`` with ``found`` False and both
    coordinates zero.
    """

    z_neck: float
    rho_neck: float
    found: bool
    status: Status

    @property
    def ok(self) -> bool:
        """True when the search ran to completion."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


@dataclass(frozen=True)
class StarConvexityResult:
    """Optimum of the star-convexity test over the origin shift.

    ``z_shift_total`` is the intrinsic COM shift plus the optimum s*, and
    ``g_opt`` is the test value there (negative for a star-convex shape).
    """

    z_shift_total: float
    g_opt: float
    status: Status

    @property
    def ok(self) -> bool:
        """True when the search ran to completion."""
        return self.status == Status.valid

    @property
    def message(self) -> str:
        """Description of :attr:`status`."""
        return status_message(self.status)


# --- helpers ---------------------------------------------------------------

def theta_grid(n: int) -> npt.NDArray[np.float64]:
    """The open uniform theta grid ``theta_i = i*pi/(n+1)``, ``i = 1..n``.

    Open: neither pole is a node. Pole radii come analytically from
    :func:`shape` / :meth:`Cache.shape`.

    Every theta handed to the library must lie inside ``[0, pi]``; anything
    outside is rejected with ``Status.invalid_grid``. This construction cannot
    overshoot — the largest node is ``n*pi/(n+1) < pi``. An INCLUSIVE grid
    built as ``i*pi/(n-1)`` can round one ulp ABOVE pi at its last node, so pin
    both endpoints if you build one yourself.

    Parameters
    ----------
    n : int
        Number of interior points.

    Returns
    -------
    numpy.ndarray
        ``n`` thetas in radians, strictly inside ``(0, pi)``.
    """
    return np.arange(1, int(n) + 1, dtype=np.float64) * np.pi / (int(n) + 1.0)


def _as_1d(values: npt.ArrayLike, name: str) -> npt.NDArray[np.float64]:
    """Contiguous 1-D float64 view of ``values``; FosParamError otherwise."""
    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise FosParamError(f"{name} must be 1-D, got shape {arr.shape}")
    return arr


def _ptr(arr: npt.NDArray[np.float64]) -> ctypes.POINTER(ctypes.c_double):  # type: ignore[valid-type]
    """Raw double pointer to a contiguous float64 array."""
    return arr.ctypes.data_as(c_dbl_p)


# --- flat tier-1 functions -------------------------------------------------

def radius_grid(params: npt.ArrayLike, thetas: npt.ArrayLike,
                n_points: int) -> RadiusGridResult:
    """R(theta) for one shape, tables built and discarded per call.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries (a2 follows from the
        volume constraint).
    thetas : array_like
        Evaluation angles in radians, all inside ``[0, pi]``.
    n_points : int
        u-grid resolution of the internal tables, at least ``N_POINTS_FLOOR``.

    Returns
    -------
    RadiusGridResult
        Radii in the total-shift frame, and the status; radii are zero-filled
        on failure.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = np.zeros(t.size, dtype=np.float64)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_radius_grid(
        _ptr(p), p.size, _ptr(t), t.size, int(n_points), _ptr(radii),
        ctypes.byref(status))
    return RadiusGridResult(radii=radii, status=Status(status.value))


def radius_and_derivative(params: npt.ArrayLike, thetas: npt.ArrayLike,
                          n_points: int) -> RadiusDerivativeResult:
    """R(theta) and dR/dtheta for one shape, tables discarded per call.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.
    thetas : array_like
        Evaluation angles in radians, all inside ``[0, pi]``.
    n_points : int
        u-grid resolution of the internal tables.

    Returns
    -------
    RadiusDerivativeResult
        Radii, derivatives and status; both buffers zero-filled on failure.
    """
    p = _as_1d(params, "params")
    t = _as_1d(thetas, "thetas")
    radii = np.zeros(t.size, dtype=np.float64)
    dr_dtheta = np.zeros(t.size, dtype=np.float64)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_radius_and_derivative(
        _ptr(p), p.size, _ptr(t), t.size, int(n_points), _ptr(radii),
        _ptr(dr_dtheta), ctypes.byref(status))
    return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta,
                                  status=Status(status.value))


def shape(params: npt.ArrayLike, n_points: int) -> FosShape:
    """Resolve one shape without evaluating a grid.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.
    n_points : int
        u-grid resolution of the internal tables.

    Returns
    -------
    FosShape
        Total z-shift, ``R(0)``, ``R(pi)`` and the status.
    """
    p = _as_1d(params, "params")
    z_shift_out = ctypes.c_double(0.0)
    r_north = ctypes.c_double(0.0)
    r_south = ctypes.c_double(0.0)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_shape(
        _ptr(p), p.size, int(n_points), ctypes.byref(z_shift_out),
        ctypes.byref(r_north), ctypes.byref(r_south), ctypes.byref(status))
    return FosShape(z_shift=z_shift_out.value, r_north=r_north.value,
                    r_south=r_south.value, status=Status(status.value))


def rho_z_grid(params: npt.ArrayLike, n_points: int) -> RhoZGridResult:
    """Cylindrical profile rho(z) on the internal u-grid, in the COM frame.

    Star-convexity and the beak gate do NOT apply here, so a shape the
    R(theta) conversion rejects still yields a profile — this is the plotting
    path.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.
    n_points : int
        Number of profile nodes, at least ``N_POINTS_FLOOR``.

    Returns
    -------
    RhoZGridResult
        ``n_points``-long z, rho and drho/dz, the intrinsic COM shift, and the
        status; all buffers zero-filled on failure.
    """
    p = _as_1d(params, "params")
    n = max(int(n_points), 0)
    z = np.zeros(n, dtype=np.float64)
    rho = np.zeros(n, dtype=np.float64)
    drho_dz = np.zeros(n, dtype=np.float64)
    z_shift_out = ctypes.c_double(0.0)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_rho_z_grid(
        _ptr(p), p.size, int(n_points), _ptr(z), _ptr(rho), _ptr(drho_dz),
        ctypes.byref(z_shift_out), ctypes.byref(status))
    return RhoZGridResult(z=z, rho=rho, drho_dz=drho_dz,
                          z_shift=z_shift_out.value, status=Status(status.value))


def neck(params: npt.ArrayLike, n_points: int) -> NeckResult:
    """Newton-refined neck of one shape, in the COM frame.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.
    n_points : int
        u-grid resolution the search is seeded from.

    Returns
    -------
    NeckResult
        Neck coordinates, whether one exists, and the status.
    """
    p = _as_1d(params, "params")
    z_neck = ctypes.c_double(0.0)
    rho_neck = ctypes.c_double(0.0)
    found = ctypes.c_int32(0)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_neck(
        _ptr(p), p.size, int(n_points), ctypes.byref(z_neck),
        ctypes.byref(rho_neck), ctypes.byref(found), ctypes.byref(status))
    return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value,
                      found=bool(found.value), status=Status(status.value))


def star_convexity_optimum(params: npt.ArrayLike,
                           n_points: int) -> StarConvexityResult:
    """Best origin shift for the star-convexity test, and the value there.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.
    n_points : int
        u-grid resolution of the internal tables.

    Returns
    -------
    StarConvexityResult
        Total shift at the optimum, ``g(s*)``, and the status.
    """
    p = _as_1d(params, "params")
    z_shift_total = ctypes.c_double(0.0)
    g_opt = ctypes.c_double(0.0)
    status = ctypes.c_int32(0)
    _get_lib().fos_param_star_convexity_optimum(
        _ptr(p), p.size, int(n_points), ctypes.byref(z_shift_total),
        ctypes.byref(g_opt), ctypes.byref(status))
    return StarConvexityResult(z_shift_total=z_shift_total.value,
                               g_opt=g_opt.value, status=Status(status.value))


def z_shift(params: npt.ArrayLike) -> float:
    """Intrinsic COM z-shift in closed form.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.

    Returns
    -------
    float
        The shift, or 0.0 if the library rejected the vector (too many
        parameters is the only rejection this closed form has). Use
        :func:`shape` when you need a status.
    """
    p = _as_1d(params, "params")
    out = ctypes.c_double(0.0)
    _get_lib().fos_param_z_shift(_ptr(p), p.size, ctypes.byref(out), None)
    return out.value


def a2(params: npt.ArrayLike) -> float:
    """``a2`` from the volume-conservation constraint.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``, 1 .. ``MAX_PARAMS`` entries.

    Returns
    -------
    float
        The constrained coefficient, or 0.0 if the library rejected the vector.
    """
    p = _as_1d(params, "params")
    out = ctypes.c_double(0.0)
    _get_lib().fos_param_a2(_ptr(p), p.size, ctypes.byref(out), None)
    return out.value


def rho_at_z(params: npt.ArrayLike, z: float,
             z_shift_value: float) -> Tuple[float, float]:
    """rho and drho/dz at one z — raw, unvalidated, status-free.

    Outside the tier rules: no tables, no gates. ``z`` is in the SHIFTED frame,
    so the FoS coordinate is ``u = (z - z_shift_value) / c``. Degenerate input
    (empty vector, c at or below the minimum, |u| at a tip within roundoff)
    yields the documented fallback ``(0.0, 0.0)``.

    Parameters
    ----------
    params : array_like
        ``[c, a3, a4, ...]``; zero-extended past its end, no length cap.
    z : float
        Axial coordinate in the shifted frame.
    z_shift_value : float
        The shift that defines that frame, e.g. from :func:`z_shift`.

    Returns
    -------
    tuple of float
        ``(rho, drho_dz)``.
    """
    p = _as_1d(params, "params")
    rho = ctypes.c_double(0.0)
    drho_dz = ctypes.c_double(0.0)
    _get_lib().fos_param_rho_at_z(
        _ptr(p), p.size, float(z), float(z_shift_value),
        ctypes.byref(rho), ctypes.byref(drho_dz))
    return rho.value, drho_dz.value


# --- cached tier -----------------------------------------------------------

class Cache:
    """Per-shape working state over a fixed resolution and theta set.

    Create once, evaluate many shapes: only the intermediates invalidated by
    the changed parameters are recomputed. THREAD-CONFINED — every compute
    call mutates the cache, so give each thread its own.

    The handle is released by :meth:`close`, by leaving a ``with`` block, or at
    garbage collection.

    Parameters
    ----------
    n_params : int
        Exact length every later ``params`` array must have,
        1 .. ``CACHE_MAX_PARAMS``. A different length comes back as
        ``Status.wrong_param_count``.
    n_points : int
        u-grid resolution, at least ``N_POINTS_FLOOR``.
    thetas : array_like
        Primary evaluation angles in radians, at least one, all inside
        ``[0, pi]``. See :func:`theta_grid`.

    Raises
    ------
    FosParamError
        If the library refuses the arguments and returns no handle.
    """

    def __init__(self, n_params: int, n_points: int,
                 thetas: npt.ArrayLike) -> None:
        lib = _get_lib()
        t = _as_1d(thetas, "thetas")
        handle = lib.fos_param_cache_create(
            int(n_params), int(n_points), _ptr(t), t.size)
        if not handle:
            raise FosParamError(
                f"Cache creation failed (n_params={n_params}, "
                f"n_points={n_points}, n_thetas={t.size}): need "
                f"1 <= n_params <= {CACHE_MAX_PARAMS}, "
                f"n_points >= {N_POINTS_FLOOR}, at least one theta in [0, pi]")
        self._handle: Optional[ctypes.c_void_p] = ctypes.c_void_p(handle)
        # Bound here so __del__ never needs module globals during shutdown.
        self._destroy = lib.fos_param_cache_destroy
        self.n_params = int(n_params)
        self.n_points = int(n_points)
        self.n_thetas = int(t.size)

    def close(self) -> None:
        """Destroy the underlying handle. Idempotent."""
        handle = getattr(self, "_handle", None)
        if handle is not None:
            self._handle = None
            self._destroy(handle)

    def __enter__(self) -> "Cache":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    def _require_handle(self) -> ctypes.c_void_p:
        handle = getattr(self, "_handle", None)
        if handle is None:
            raise FosParamError("Cache is closed")
        return handle

    def radius_grid(self, params: npt.ArrayLike) -> RadiusGridResult:
        """R(theta) at the cache's own thetas, in the total-shift frame."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = np.zeros(self.n_thetas, dtype=np.float64)
        status = _get_lib().fos_param_cache_radius_grid(
            handle, _ptr(p), p.size, _ptr(radii), radii.size)
        return RadiusGridResult(radii=radii, status=Status(status))

    def radius_and_derivative(self, params: npt.ArrayLike) -> RadiusDerivativeResult:
        """R(theta) and dR/dtheta at the cache's own thetas."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        radii = np.zeros(self.n_thetas, dtype=np.float64)
        dr_dtheta = np.zeros(self.n_thetas, dtype=np.float64)
        status = _get_lib().fos_param_cache_radius_and_derivative(
            handle, _ptr(p), p.size, _ptr(radii), _ptr(dr_dtheta), radii.size)
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta,
                                      status=Status(status))

    def radius_and_derivative_at_thetas(
            self, params: npt.ArrayLike,
            thetas: npt.ArrayLike) -> RadiusDerivativeResult:
        """R and dR/dtheta at caller-supplied thetas.

        Uncached evaluation against the cache's resolved shape: the thetas need
        not be the cache's own, but must lie inside ``[0, pi]``.

        Parameters
        ----------
        params : array_like
            Exactly ``n_params`` entries.
        thetas : array_like
            Evaluation angles in radians.

        Returns
        -------
        RadiusDerivativeResult
            Buffers as long as ``thetas``; zero-filled on failure.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        t = _as_1d(thetas, "thetas")
        radii = np.zeros(t.size, dtype=np.float64)
        dr_dtheta = np.zeros(t.size, dtype=np.float64)
        status = _get_lib().fos_param_cache_radius_and_derivative_at_thetas(
            handle, _ptr(p), p.size, _ptr(t), t.size, _ptr(radii), _ptr(dr_dtheta))
        return RadiusDerivativeResult(radii=radii, dr_dtheta=dr_dtheta,
                                      status=Status(status))

    def shape(self, params: npt.ArrayLike) -> FosShape:
        """Resolve a shape: total z-shift and the analytic pole radii."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z_shift_out = ctypes.c_double(0.0)
        r_north = ctypes.c_double(0.0)
        r_south = ctypes.c_double(0.0)
        status = _get_lib().fos_param_cache_shape(
            handle, _ptr(p), p.size, ctypes.byref(z_shift_out),
            ctypes.byref(r_north), ctypes.byref(r_south))
        return FosShape(z_shift=z_shift_out.value, r_north=r_north.value,
                        r_south=r_south.value, status=Status(status))

    def rho_z_grid(self, params: npt.ArrayLike) -> RhoZGridResult:
        """Cylindrical profile rho(z) on the cache's u-grid, in the COM frame.

        Neither the beak gate nor star-convexity applies, so this path renders
        shapes :meth:`radius_grid` rejects.
        """
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z = np.zeros(self.n_points, dtype=np.float64)
        rho = np.zeros(self.n_points, dtype=np.float64)
        drho_dz = np.zeros(self.n_points, dtype=np.float64)
        z_shift_out = ctypes.c_double(0.0)
        status = _get_lib().fos_param_cache_rho_z_grid(
            handle, _ptr(p), p.size, _ptr(z), _ptr(rho), _ptr(drho_dz),
            z.size, ctypes.byref(z_shift_out))
        return RhoZGridResult(z=z, rho=rho, drho_dz=drho_dz,
                              z_shift=z_shift_out.value, status=Status(status))

    def neck(self, params: npt.ArrayLike) -> NeckResult:
        """Newton-refined neck in the COM frame."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z_neck = ctypes.c_double(0.0)
        rho_neck = ctypes.c_double(0.0)
        found = ctypes.c_int32(0)
        status = _get_lib().fos_param_cache_neck(
            handle, _ptr(p), p.size, ctypes.byref(z_neck),
            ctypes.byref(rho_neck), ctypes.byref(found))
        return NeckResult(z_neck=z_neck.value, rho_neck=rho_neck.value,
                          found=bool(found.value), status=Status(status))

    def star_convexity_optimum(self, params: npt.ArrayLike) -> StarConvexityResult:
        """Best origin shift for the star-convexity test, and ``g(s*)``."""
        handle = self._require_handle()
        p = _as_1d(params, "params")
        z_shift_total = ctypes.c_double(0.0)
        g_opt = ctypes.c_double(0.0)
        status = _get_lib().fos_param_cache_star_convexity_optimum(
            handle, _ptr(p), p.size, ctypes.byref(z_shift_total),
            ctypes.byref(g_opt))
        return StarConvexityResult(z_shift_total=z_shift_total.value,
                                   g_opt=g_opt.value, status=Status(status))


__all__ = [
    "Cache", "FosParamError", "FosShape", "NeckResult",
    "RadiusDerivativeResult", "RadiusGridResult", "RhoZGridResult",
    "StarConvexityResult", "Status",
    "a2", "neck", "radius_and_derivative", "radius_grid", "rho_at_z",
    "rho_z_grid", "shape", "star_convexity_optimum", "status_message",
    "theta_grid", "z_shift",
    "CACHE_MAX_PARAMS", "MAX_PARAMS", "N_POINTS_FLOOR",
]
