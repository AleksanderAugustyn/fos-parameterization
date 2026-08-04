"""ctypes signatures for the fos-parameterization 2.0.0 C API.

One entry per exported symbol in ``include/fos_parameterization.h``. Every
integer the C API takes or returns is ``int32_t``; every buffer is ``double*``.
"""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_i32 = ctypes.c_int32
c_i32_p = ctypes.POINTER(ctypes.c_int32)
c_void = ctypes.c_void_p

#: Highest ``n_params`` the flat tier-1 calls accept (``FOS_PARAM_MAX_PARAMS``).
MAX_PARAMS = 50
#: Highest ``n_params`` a :class:`~fos_parameterization.api.Cache` accepts.
CACHE_MAX_PARAMS = 8
#: Lowest u-grid resolution a tables/cache handle may be built for.
N_POINTS_FLOOR = 100

# name -> (argtypes, restype)
_SIGNATURES = {
    # --- diagnostics ---
    "fos_param_status_message": ([c_i32], ctypes.c_char_p),

    # --- handle lifecycle ---
    "fos_param_tables_create": ([c_i32, c_dbl_p, c_i32], c_void),
    "fos_param_tables_destroy": ([c_void], None),
    "fos_param_cache_create": ([c_i32, c_i32, c_dbl_p, c_i32], c_void),
    "fos_param_cache_create_shared": ([c_void, c_i32], c_void),
    "fos_param_cache_destroy": ([c_void], None),

    # --- cached computes (tier 2), status is the return value ---
    "fos_param_cache_radius_grid": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_i32], c_i32),
    "fos_param_cache_radius_and_derivative": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32], c_i32),
    "fos_param_cache_radius_and_derivative_at_thetas": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_i32, c_dbl_p, c_dbl_p], c_i32),
    "fos_param_cache_shape": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_dbl_p], c_i32),
    "fos_param_cache_rho_z_grid": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32, c_dbl_p], c_i32),
    "fos_param_cache_neck": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p, c_i32_p], c_i32),
    "fos_param_cache_star_convexity_optimum": (
        [c_void, c_dbl_p, c_i32, c_dbl_p, c_dbl_p], c_i32),

    # --- flat tier-1 computes, trailing nullable status out-parameter ---
    "fos_param_radius_grid": (
        [c_dbl_p, c_i32, c_dbl_p, c_i32, c_i32, c_dbl_p, c_i32_p], None),
    "fos_param_radius_and_derivative": (
        [c_dbl_p, c_i32, c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_i32_p], None),
    "fos_param_shape": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "fos_param_rho_z_grid": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_dbl_p, c_dbl_p, c_i32_p], None),
    "fos_param_neck": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_i32_p, c_i32_p], None),
    "fos_param_star_convexity_optimum": (
        [c_dbl_p, c_i32, c_i32, c_dbl_p, c_dbl_p, c_i32_p], None),
    "fos_param_z_shift": ([c_dbl_p, c_i32, c_dbl_p, c_i32_p], None),
    "fos_param_a2": ([c_dbl_p, c_i32, c_dbl_p, c_i32_p], None),

    # --- raw evaluator, no status ---
    "fos_param_rho_at_z": (
        [c_dbl_p, c_i32, ctypes.c_double, ctypes.c_double, c_dbl_p, c_dbl_p], None),
}


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restype on every exported symbol (idempotent).

    Parameters
    ----------
    lib : ctypes.CDLL
        The loaded shared library.

    Returns
    -------
    ctypes.CDLL
        The same object, with every signature attached.

    Raises
    ------
    AttributeError
        If the library predates the 2.0.0 symbol set.
    """
    for name, (argtypes, restype) in _SIGNATURES.items():
        func = getattr(lib, name)
        func.argtypes = argtypes
        func.restype = restype
    return lib
