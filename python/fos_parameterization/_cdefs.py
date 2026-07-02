"""ctypes signatures for the fos-parameterization C API."""
from __future__ import annotations

import ctypes

c_dbl_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int)

MESSAGE_BUFFER_SIZE = 256


def configure(lib: ctypes.CDLL) -> ctypes.CDLL:
    """Set argtypes/restypes on the loaded library (idempotent)."""
    lib.fos_compute_radius_grid.argtypes = [
        c_dbl_p, ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.fos_compute_radius_grid.restype = ctypes.c_int

    lib.fos_compute_rho_profile.argtypes = [
        c_dbl_p, ctypes.c_int, ctypes.c_int, c_dbl_p, c_dbl_p, c_dbl_p,
        ctypes.c_int, ctypes.c_char_p]
    lib.fos_compute_rho_profile.restype = ctypes.c_int

    lib.fos_compute_neck.argtypes = [
        c_dbl_p, ctypes.c_int, c_dbl_p, c_dbl_p, c_int_p]
    lib.fos_compute_neck.restype = ctypes.c_int

    lib.fos_z_shift.argtypes = [c_dbl_p, ctypes.c_int]
    lib.fos_z_shift.restype = ctypes.c_double

    lib.fos_a2.argtypes = [c_dbl_p, ctypes.c_int]
    lib.fos_a2.restype = ctypes.c_double
    return lib
