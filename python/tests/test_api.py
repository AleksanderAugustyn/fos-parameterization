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

# F5: marginal star-convex shape recovered by the golden-section shift search
# (best origin gives max-T = -0.103 R0, just past the -0.1 margin; the old coarse
# shift search rejected it). z_shift = -0.16905 R0 is the true optimum.
F5_PARAMS = [2.00, 0.40, 0.66, 0.0, 0.0, 0.0, 0.0]
F5_EXPECTED = [
    1.8309533690419022E+000,
    1.2362050765905319E+000,
    2.2194559887019227E-001,
    1.7533696725021786E-001,
    2.5349878018568273E-001,
    1.7128255699392803E+000,
    2.1690466309580976E+000,
]
F5_Z_SHIFT = -1.6904663095809774E-001

# Shape-split + derivative goldens (same capture program, n_rho_grid = 181).
# Thetas are pre-rounded literals for pi/8, pi/2, 7pi/8 — bit-identical to the
# literals in golden_capture.f08, so both sides feed the same doubles.
DERIV_THETAS = [0.39269908169872414, 1.5707963267948966, 2.748893571891069]

# dR/dtheta is a direct arithmetic expression of (r, theta). Unlike R — a
# Newton fixed point that self-corrects to the same root — it inherits
# FMA/contraction differences between the capture binary and the shared
# library (Release is -ffast-math -flto: two separately optimized codegens),
# leaving a few-ulp floor (~1.2e-15 observed) between the two.
DERIV_DR_ATOL = 5e-15

F2_R_NORTH = 1.9675901550757660E+000
F2_R_SOUTH = 1.6324098449242341E+000
F2_DERIV_R = [
    1.5387379214194603E+000,
    6.2134923161743805E-001,
    1.4646885007634520E+000,
]
F2_DERIV_DR = [
    -2.0821776359437809E+000,
    1.5235663870749072E-001,
    7.3011136010047017E-001,
]

F3_R_NORTH = 1.5716197243913530E+000
F3_R_SOUTH = 1.4283802756086470E+000
F3_DERIV_R = [
    1.3915669254733638E+000,
    7.0768449215772455E-001,
    1.3317901145714259E+000,
]
F3_DERIV_DR = [
    -7.9984191593394949E-001,
    9.6279565037092382E-002,
    4.3222840642498583E-001,
]

F4_R_NORTH = 2.0000000000000000E+000
F4_R_SOUTH = 2.0000000000000000E+000
F4_DERIV_R = [
    1.6553976256663301E+000,
    4.0824829046386313E-001,
    1.6553976256663301E+000,
]
F4_DERIV_DR = [
    -1.5307124866118125E+000,
    -3.9863274022861999E-017,
    1.5307124866118125E+000,
]

F5_R_NORTH = 1.8309533690419022E+000
F5_R_SOUTH = 2.1690466309580976E+000
F5_DERIV_R = [
    1.4857936540921055E+000,
    1.7533696725021786E-001,
    1.9105555374098624E+000,
]
F5_DERIV_DR = [
    -1.5693450137431226E+000,
    1.6198246115502524E-002,
    1.2205515080911753E+000,
]


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
    (F5_PARAMS, F5_EXPECTED, F5_Z_SHIFT),
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


def test_shape_and_derivative_grid_parity() -> None:
    params = [1.5, 0.08, 0.05]
    n_grid = 91
    ref = fp.radius_grid(params, n_grid)
    assert ref.ok

    shp = fp.shape(params, n_rho_grid=n_grid)
    assert shp.ok
    assert shp.z_shift == pytest.approx(ref.z_shift, abs=1e-15)
    assert shp.r_north > 0.0 and shp.r_south > 0.0

    thetas = fp.theta_grid(n_grid)
    result = fp.radius_and_derivative(params, thetas, shp.z_shift)
    assert result.ok
    np.testing.assert_allclose(result.radii, ref.radii, rtol=0.0, atol=1e-15)
    assert result.dr_dtheta[0] == 0.0 and result.dr_dtheta[-1] == 0.0


def test_derivative_vs_fd() -> None:
    params = [1.5, 0.08, 0.05]
    shp = fp.shape(params, n_rho_grid=1000)
    thetas = np.linspace(0.2, np.pi - 0.2, 50)
    h = 1e-3

    r = fp.radius_and_derivative(params, thetas, shp.z_shift)
    fd_stencils = [
        fp.radius_and_derivative(params, thetas + k * h, shp.z_shift).radii
        for k in (-2, -1, 1, 2)
    ]
    fd = (fd_stencils[0] - 8 * fd_stencils[1] + 8 * fd_stencils[2] - fd_stencils[3]) / (12 * h)
    np.testing.assert_allclose(r.dr_dtheta, fd, rtol=0.0, atol=1e-9)


def test_shape_invalid_c() -> None:
    shp = fp.shape([-1.0], n_rho_grid=1000)
    assert not shp.ok
    assert shp.status == fp.Status.ERROR_INVALID_C
    assert shp.z_shift == 0.0 and shp.r_north == 0.0 and shp.r_south == 0.0


@pytest.mark.parametrize("params,z_shift,r_north,r_south,deriv_r,deriv_dr", [
    (F2_PARAMS, F2_Z_SHIFT, F2_R_NORTH, F2_R_SOUTH, F2_DERIV_R, F2_DERIV_DR),
    (F3_PARAMS, F3_Z_SHIFT, F3_R_NORTH, F3_R_SOUTH, F3_DERIV_R, F3_DERIV_DR),
    (F4_PARAMS, F4_Z_SHIFT, F4_R_NORTH, F4_R_SOUTH, F4_DERIV_R, F4_DERIV_DR),
    (F5_PARAMS, F5_Z_SHIFT, F5_R_NORTH, F5_R_SOUTH, F5_DERIV_R, F5_DERIV_DR),
])
def test_derivative_goldens(params, z_shift, r_north, r_south, deriv_r, deriv_dr) -> None:
    shp = fp.shape(params, n_rho_grid=181)
    assert shp.ok
    assert shp.z_shift == pytest.approx(z_shift, rel=1e-15, abs=1e-300)
    np.testing.assert_allclose([shp.r_north, shp.r_south], [r_north, r_south],
                               rtol=0.0, atol=1e-15)

    res = fp.radius_and_derivative(params, DERIV_THETAS, shp.z_shift)
    assert res.ok
    np.testing.assert_allclose(res.radii, deriv_r, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(res.dr_dtheta, deriv_dr, rtol=0.0, atol=DERIV_DR_ATOL)


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
