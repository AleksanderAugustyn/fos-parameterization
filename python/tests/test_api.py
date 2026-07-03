"""Python == Fortran: bindings must reproduce the library bit-for-bit."""
from __future__ import annotations

import math

import numpy as np
import pytest

import fos_parameterization as fp

# Golden shapes captured by tests/golden_capture.f08 (17 sig digits).
# n_grid = 181; indices below are 0-based equivalents of Fortran [1,31,...,181].
GOLDEN_IDX = [0, 30, 60, 90, 120, 150, 180]

F2_PARAMS = [1.80, 0.20, 0.30, 0.01, -0.02, 0.0, 0.0]
F2_EXPECTED = [
    1.9675901550757660E+000,
    1.2229761015937455E+000,
    6.5003108135859899E-001,
    6.2134923161743805E-001,
    8.4183035874080525E-001,
    1.3592368734106997E+000,
    1.6324098449242341E+000,
]
F2_Z_SHIFT = 1.6759015507576580E-001

F3_PARAMS = [1.50, 0.10, 0.20, 0.0, 0.0, 0.0, 0.0]
F3_EXPECTED = [
    1.5716197243913530E+000,
    1.2756928073264713E+000,
    8.0180080738517090E-001,
    7.0768449215772455E-001,
    9.0700239126619542E-001,
    1.2690634114059940E+000,
    1.4283802756086470E+000,
]
F3_Z_SHIFT = 7.1619724391352904E-002

F4_PARAMS = [2.00, 0.00, 0.50, 0.0, 0.0, 0.0, 0.0]
F4_EXPECTED = [
    2.0000000000000000E+000,
    1.4225493540227652E+000,
    5.2408493336619644E-001,
    4.0824829046386313E-001,
    5.2408493336619610E-001,
    1.4225493540227652E+000,
    2.0000000000000000E+000,
]
F4_Z_SHIFT = 0.0000000000000000E+000


def test_sphere_exact() -> None:
    res = fp.radius_grid([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], n_grid=181)
    assert res.status == fp.Status.VALID
    assert res.ok
    np.testing.assert_allclose(res.radii, np.ones(181), rtol=0.0, atol=5e-12)
    assert res.z_shift == 0.0


@pytest.mark.parametrize("params,expected,z_shift", [
    (F2_PARAMS, F2_EXPECTED, F2_Z_SHIFT),
    (F3_PARAMS, F3_EXPECTED, F3_Z_SHIFT),
    (F4_PARAMS, F4_EXPECTED, F4_Z_SHIFT),
])
def test_goldens_match_fortran(params, expected, z_shift) -> None:
    res = fp.radius_grid(params, n_grid=181)
    assert res.status == fp.Status.VALID
    np.testing.assert_allclose(res.radii[GOLDEN_IDX], expected, rtol=1e-15, atol=0.0)
    assert res.z_shift == pytest.approx(z_shift, rel=1e-15, abs=1e-300)


def test_error_codes_surface_as_status() -> None:
    res = fp.radius_grid([0.0, 0.0, 0.0], n_grid=41)          # c = 0
    assert res.status == fp.Status.ERROR_INVALID_C
    assert not res.ok
    assert res.message                                         # non-empty
    assert np.all(res.radii == 0.0) and res.z_shift == 0.0

    res = fp.radius_grid([2.0, 0.0, 0.7497], n_grid=41)        # beak
    assert res.status == fp.Status.ERROR_BEAK_SINGULARITY

    res = fp.radius_grid([1.0], n_grid=1)                      # bad n_grid
    assert res.status == fp.Status.ERROR_INVALID_ARGUMENTS


def test_rho_profile_sphere() -> None:
    res = fp.rho_profile([1.0, 0.0, 0.0], n_z=201)
    assert res.status == fp.Status.VALID
    mid = 100
    assert res.z[mid] == pytest.approx(0.0, abs=1e-14)
    assert res.rho[mid] == pytest.approx(1.0, abs=1e-12)
    assert res.rho[0] == 0.0 and res.rho[-1] == 0.0
    assert res.z.shape == res.rho.shape == res.drho_dz.shape == (201,)


def test_neck_analytic() -> None:
    res = fp.neck([2.0, 0.0, 0.5])
    assert res.status == fp.Status.VALID
    assert res.found
    assert res.rho_neck == pytest.approx(math.sqrt(1.0 / 6.0), abs=1e-9)
    assert res.z_neck == pytest.approx(0.0, abs=1e-9)

    res = fp.neck([1.0, 0.0, 0.0])                             # sphere: no neck
    assert res.status == fp.Status.VALID
    assert not res.found


def test_scalar_helpers() -> None:
    assert fp.z_shift([2.0, 0.0, 0.5]) == 0.0                  # symmetric
    assert fp.a2([2.0, 0.0, 0.5]) == pytest.approx(0.5 / 3.0, rel=1e-15)
    assert fp.z_shift([1.80, 0.20, 0.30, 0.01, -0.02]) != 0.0  # asymmetric


def test_theta_grid() -> None:
    t = fp.theta_grid(181)
    assert t.shape == (181,)
    assert t[0] == 0.0
    assert t[-1] == pytest.approx(np.pi, rel=1e-15)


def test_params_must_be_1d() -> None:
    with pytest.raises(ValueError):
        fp.radius_grid(np.zeros((2, 3)), n_grid=41)


@pytest.mark.parametrize("params", [
    [1.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0],  # z_shift != 0: left tip lands 1 ulp inside
    [1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # dz accumulation: right tip 1 ulp inside
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # sphere: clean before the fix, must stay clean
])
def test_profile_tips_are_exact_zero(params) -> None:
    """f(±1) = 0 by construction, so tips are rho = 0 with the drho_dz = 0 convention.

    Regression (2026-07-03): the exact abs(u) >= 1 tip guard let ~1-ulp endpoint
    residue through; f evaluated at roundoff scale and drho_dz = f'/(2c*sqrt(c*f))
    exploded to ~1e7 on whichever tip the residue landed inside.
    """
    prof = fp.rho_profile(params, 721)
    assert prof.ok
    for i in (0, -1):
        assert prof.rho[i] == 0.0
        assert prof.drho_dz[i] == 0.0
