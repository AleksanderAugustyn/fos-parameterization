"""Contract tests for the 3.0.0 Python surface.

Shape-validation failures are results, not exceptions: the buffers come back
zero-filled with a nonzero :class:`Status`. Only usage errors — a failed
create, a closed handle, non-1-D input — raise ``FosParamError``.

The golden blocks are bit-comparisons against the Fortran library, captured by
``tests/golden_capture.f08`` through the same flat tier-1 entries these tests
call. Recapture procedure lives in that file's header.
"""
from __future__ import annotations

import math
import threading

import numpy as np
import pytest

import fos_parameterization as fp

# --- fixtures --------------------------------------------------------------

SPHERE = [1.0, 0.0, 0.0]
N_POINTS = 181           # the golden-capture resolution
ASYMMETRIC = [1.80, 0.20, 0.30, 0.01, -0.02, 0.0, 0.0]
# Beak-rejected by the R(theta) conversion, still a drawable rho(z) profile.
# f(0) = 1 - 4*a4/3 = 2.0e-4 < F_MIN_THRESHOLD = 5.0e-4 (2026-08-11 retune)
BEAK = [2.0, 0.0, 0.74985]
# Two fragments: f(0) = 1 - 4*a4/3 = -0.2.
SPLIT = [2.0, 0.0, 0.9]
TOO_MANY = [1.0] + [0.0] * 50   # 51 entries, one past FOS_PARAM_MAX_PARAMS


def _inclusive_grid(n: int) -> np.ndarray:
    """The golden-capture grid: ``n`` uniform nodes on the CLOSED ``[0, pi]``.

    Bit-identical to the loop in ``tests/golden_capture.f08`` — multiply then
    divide, both endpoints pinned so no rounding can push a node past pi.

    Parameters
    ----------
    n : int
        Node count, at least 2.

    Returns
    -------
    numpy.ndarray
        ``n`` thetas with ``grid[0] == 0.0`` and ``grid[-1] == pi`` exactly.
    """
    grid = np.arange(n, dtype=np.float64) * np.pi / (n - 1.0)
    grid[0] = 0.0
    grid[-1] = np.pi
    return grid


def _same_bits(a: np.ndarray, b: np.ndarray) -> bool:
    """True when ``a`` and ``b`` hold the same IEEE754 bit patterns.

    ``np.array_equal`` compares values, so it calls ``+0.0`` and ``-0.0``
    equal; the contract's bitwise claims need the bytes. In the wheel job this
    suite is the only bitwise check on the compiler that ships.
    """
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and a.dtype == b.dtype and a.tobytes() == b.tobytes()


# Golden shapes captured by tests/golden_capture.f08 (Debug build, 17 sig
# digits). n_points = 181; the indices are 0-based equivalents of the Fortran
# [1, 31, ..., 181].
GOLDEN_IDX = [0, 30, 60, 90, 120, 150, 180]

# Thetas as pre-rounded literals for pi/8, pi/2, 7pi/8 — bit-identical to the
# literals in golden_capture.f08, so both sides feed the same doubles.
DERIV_THETAS = [0.39269908169872414, 1.5707963267948966, 2.748893571891069]

# dR/dtheta is a direct arithmetic expression of (r, theta). Unlike R — a
# Newton fixed point that self-corrects to the same root — it inherits
# FMA/contraction differences between the capture binary and the shared
# library, and between Debug and Release codegen (Release is -ffast-math
# -flto), leaving a few-ulp floor (~2.2e-15 observed) between the two.
DERIV_DR_ATOL = 5e-15

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
F2_R_NORTH = 1.9675901550757660E+000
F2_R_SOUTH = 1.6324098449242341E+000
F2_DERIV_R = [
    1.5387379214194605E+000,
    6.2134923161743805E-001,
    1.4646885007634518E+000,
]
F2_DERIV_DR = [
    -2.0821776359437787E+000,
    1.5235663870749072E-001,
    7.3011136010047140E-001,
]

F3_PARAMS = [1.50, 0.10, 0.20, 0.0, 0.0, 0.0, 0.0]
F3_EXPECTED = [
    1.5716197243913530E+000,
    1.2756928073264715E+000,
    8.0180080738517101E-001,
    7.0768449215772455E-001,
    9.0700239126619531E-001,
    1.2690634114059938E+000,
    1.4283802756086470E+000,
]
F3_Z_SHIFT = 7.1619724391352904E-002
F3_R_NORTH = 1.5716197243913530E+000
F3_R_SOUTH = 1.4283802756086470E+000
F3_DERIV_R = [
    1.3915669254733640E+000,
    7.0768449215772455E-001,
    1.3317901145714262E+000,
]
F3_DERIV_DR = [
    -7.9984191593394893E-001,
    9.6279565037092382E-002,
    4.3222840642498511E-001,
]

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

# F5: marginal star-convex shape recovered by the golden-section shift search
# (best origin gives max-T = -0.103 R0, just past the -0.1 margin).
F5_PARAMS = [2.00, 0.40, 0.66, 0.0, 0.0, 0.0, 0.0]
F5_EXPECTED = [
    1.8309533690419022E+000,
    1.2362050765905319E+000,
    2.2194559887019230E-001,
    1.7533696725021786E-001,
    2.5349878018568284E-001,
    1.7128255699392803E+000,
    2.1690466309580976E+000,
]
F5_Z_SHIFT = -1.6904663095809774E-001
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


# --- theta_grid ------------------------------------------------------------

def test_theta_grid_is_open_and_uniform() -> None:
    thetas = fp.theta_grid(64)
    assert thetas.size == 64
    assert thetas[0] > 0.0
    assert thetas[-1] < np.pi
    assert np.allclose(np.diff(thetas), np.pi / 65.0, rtol=0.0, atol=1e-14)


def test_theta_grid_small_n_is_exact() -> None:
    assert np.allclose(fp.theta_grid(3), np.array([1.0, 2.0, 3.0]) * np.pi / 4.0)


@pytest.mark.parametrize("n", [1, 2, 3, 7, 64, 181, 721, 2001, 100001])
def test_theta_grid_never_leaves_the_domain(n: int) -> None:
    """Endpoint safety: no node may round past pi (the 1.x ulp trap).

    The library rejects a theta outside [0, pi] with Status.invalid_grid, and
    an inclusive grid built as ``(n-1)*pi/(n-1)`` can land one ulp above pi.
    theta_grid is open by construction, so this must hold for every n.
    """
    thetas = fp.theta_grid(n)
    assert thetas.size == n
    assert np.all(thetas > 0.0)
    assert np.all(thetas < np.pi)


def test_theta_grid_is_accepted_by_the_library() -> None:
    res = fp.radius_grid(SPHERE, fp.theta_grid(2001), N_POINTS)
    assert res.ok


# --- tier 1: happy path ----------------------------------------------------

def test_sphere_radii_are_unit() -> None:
    thetas = _inclusive_grid(N_POINTS)
    res = fp.radius_grid(SPHERE, thetas, N_POINTS)
    assert res.ok
    assert res.status == fp.Status.valid
    assert np.allclose(res.radii, 1.0, rtol=0.0, atol=5e-12)


def test_sphere_shape_is_unit_and_centred() -> None:
    shp = fp.shape(SPHERE, N_POINTS)
    assert shp.ok
    assert shp.z_shift == 0.0
    assert shp.r_north == pytest.approx(1.0, abs=1e-12)
    assert shp.r_south == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("params,expected,z_shift", [
    (F2_PARAMS, F2_EXPECTED, F2_Z_SHIFT),
    (F3_PARAMS, F3_EXPECTED, F3_Z_SHIFT),
    (F4_PARAMS, F4_EXPECTED, F4_Z_SHIFT),
    (F5_PARAMS, F5_EXPECTED, F5_Z_SHIFT),
])
def test_radius_goldens_match_fortran(params, expected, z_shift) -> None:
    res = fp.radius_grid(params, _inclusive_grid(N_POINTS), N_POINTS)
    assert res.ok
    np.testing.assert_allclose(res.radii[GOLDEN_IDX], expected, rtol=1e-15, atol=0.0)
    shp = fp.shape(params, N_POINTS)
    assert shp.z_shift == pytest.approx(z_shift, rel=1e-15, abs=1e-300)


@pytest.mark.parametrize("params,r_north,r_south,deriv_r,deriv_dr", [
    (F2_PARAMS, F2_R_NORTH, F2_R_SOUTH, F2_DERIV_R, F2_DERIV_DR),
    (F3_PARAMS, F3_R_NORTH, F3_R_SOUTH, F3_DERIV_R, F3_DERIV_DR),
    (F4_PARAMS, F4_R_NORTH, F4_R_SOUTH, F4_DERIV_R, F4_DERIV_DR),
    (F5_PARAMS, F5_R_NORTH, F5_R_SOUTH, F5_DERIV_R, F5_DERIV_DR),
])
def test_derivative_goldens_match_fortran(params, r_north, r_south,
                                          deriv_r, deriv_dr) -> None:
    shp = fp.shape(params, N_POINTS)
    assert shp.ok
    np.testing.assert_allclose([shp.r_north, shp.r_south], [r_north, r_south],
                               rtol=0.0, atol=1e-15)
    res = fp.radius_and_derivative(params, DERIV_THETAS, N_POINTS)
    assert res.ok
    np.testing.assert_allclose(res.radii, deriv_r, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(res.dr_dtheta, deriv_dr, rtol=0.0, atol=DERIV_DR_ATOL)


def test_derivative_matches_finite_difference() -> None:
    params = [1.5, 0.08, 0.05]
    thetas = np.linspace(0.2, np.pi - 0.2, 50)
    h = 1e-3
    res = fp.radius_and_derivative(params, thetas, 501)
    stencils = [fp.radius_and_derivative(params, thetas + k * h, 501).radii
                for k in (-2, -1, 1, 2)]
    fd = (stencils[0] - 8 * stencils[1] + 8 * stencils[2] - stencils[3]) / (12 * h)
    np.testing.assert_allclose(res.dr_dtheta, fd, rtol=0.0, atol=1e-9)


def test_rho_z_grid_sphere() -> None:
    res = fp.rho_z_grid(SPHERE, 201)
    assert res.ok
    assert res.z.shape == res.rho.shape == res.drho_dz.shape == (201,)
    assert res.z_shift == 0.0
    mid = 100
    assert res.z[mid] == pytest.approx(0.0, abs=1e-14)
    assert res.rho[mid] == pytest.approx(1.0, abs=1e-12)
    assert res.rho[0] == 0.0 and res.rho[-1] == 0.0


@pytest.mark.parametrize("params", [
    [1.0, 0.0, 0.0, 0.05, 0.0, 0.0, 0.0],  # z_shift != 0: left tip 1 ulp inside
    [1.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # dz accumulation: right tip 1 ulp inside
    [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],   # sphere: clean before the fix, must stay clean
])
def test_profile_tips_are_exact_zero(params) -> None:
    """f(+-1) = 0 by construction, so tips are rho = 0 with drho_dz = 0."""
    res = fp.rho_z_grid(params, 721)
    assert res.ok
    for i in (0, -1):
        assert res.rho[i] == 0.0
        assert res.drho_dz[i] == 0.0


def test_neck_analytic() -> None:
    res = fp.neck(F4_PARAMS, 501)
    assert res.ok
    assert res.found
    assert res.rho_neck == pytest.approx(math.sqrt(1.0 / 6.0), abs=1e-9)
    assert res.z_neck == pytest.approx(0.0, abs=1e-9)

    res = fp.neck(SPHERE, 501)
    assert res.ok
    assert not res.found
    assert res.z_neck == 0.0 and res.rho_neck == 0.0


def test_star_convexity_optimum_reports_the_applied_shift() -> None:
    """F5 only resolves because of the extra shift, so shape() reports it."""
    opt = fp.star_convexity_optimum(F5_PARAMS, N_POINTS)
    assert opt.ok
    assert opt.g_opt < 0.0
    assert opt.z_shift_total == pytest.approx(F5_Z_SHIFT, rel=1e-14)
    assert fp.shape(F5_PARAMS, N_POINTS).z_shift == pytest.approx(
        opt.z_shift_total, abs=1e-14)


def test_frames_differ_between_radius_and_rho_grids() -> None:
    """R(theta) is in the total-shift frame, rho(z) in the COM frame."""
    intrinsic = fp.z_shift(F5_PARAMS)
    assert fp.rho_z_grid(F5_PARAMS, N_POINTS).z_shift == pytest.approx(
        intrinsic, abs=1e-15)
    assert fp.shape(F5_PARAMS, N_POINTS).z_shift != pytest.approx(intrinsic)


def test_scalar_helpers() -> None:
    assert fp.z_shift(F4_PARAMS) == 0.0                  # symmetric
    assert fp.a2(F4_PARAMS) == pytest.approx(0.5 / 3.0, rel=1e-15)
    assert fp.z_shift(ASYMMETRIC) != 0.0


def test_rho_at_z_raw_evaluator() -> None:
    rho, drho_dz = fp.rho_at_z(SPHERE, 0.0, 0.0)
    assert rho == pytest.approx(1.0, abs=1e-15)
    assert drho_dz == pytest.approx(0.0, abs=1e-15)


# --- statuses, not exceptions ----------------------------------------------

def test_invalid_c_is_a_status() -> None:
    res = fp.radius_grid([1e-11], fp.theta_grid(16), N_POINTS)
    assert res.status == fp.Status.invalid_c
    assert not res.ok
    assert res.message
    assert np.all(res.radii == 0.0)


@pytest.mark.parametrize("call", [
    lambda: fp.radius_grid(TOO_MANY, fp.theta_grid(16), N_POINTS),
    lambda: fp.rho_z_grid(TOO_MANY, N_POINTS),
    lambda: fp.shape(TOO_MANY, N_POINTS),
    lambda: fp.neck(TOO_MANY, N_POINTS),
])
def test_too_many_params_is_a_status_not_a_raise(call) -> None:
    """51 params exceeds the one-shot cap: Status(1) must exist and be returned.

    Regression: the 1.x enum had no member for 1, so Status(1) raised
    ValueError instead of reporting the rejection.
    """
    res = call()
    assert res.status == fp.Status.too_many_params
    assert not res.ok
    assert res.message


@pytest.mark.parametrize("call", [
    lambda: fp.radius_grid([], fp.theta_grid(16), N_POINTS),
    lambda: fp.radius_and_derivative([], fp.theta_grid(16), N_POINTS),
    lambda: fp.rho_z_grid([], N_POINTS),
    lambda: fp.rho_z_grid_unchecked([], N_POINTS),
    lambda: fp.shape([], N_POINTS),
    lambda: fp.neck([], N_POINTS),
    lambda: fp.star_convexity_optimum([], N_POINTS),
])
def test_empty_params_is_wrong_param_count(call) -> None:
    """An empty vector is code 4 in every one-shot output (2.x returned 102)."""
    res = call()
    assert res.status == fp.Status.wrong_param_count
    assert not res.ok


def test_invalid_grid_is_a_status() -> None:
    res = fp.radius_grid(SPHERE, fp.theta_grid(16), 50)   # below the n_points floor
    assert res.status == fp.Status.invalid_grid
    assert not res.ok
    assert np.all(res.radii == 0.0)


def test_failures_zero_fill_every_output() -> None:
    shp = fp.shape([1e-11], N_POINTS)
    assert shp.status == fp.Status.invalid_c
    assert shp.z_shift == 0.0 and shp.r_north == 0.0 and shp.r_south == 0.0

    res = fp.rho_z_grid([1e-11], N_POINTS)
    assert not res.ok
    assert np.all(res.z == 0.0) and np.all(res.rho == 0.0)
    assert np.all(res.drho_dz == 0.0) and res.z_shift == 0.0

    nck = fp.neck([1e-11], N_POINTS)
    assert not nck.ok
    assert nck.z_neck == 0.0 and nck.rho_neck == 0.0 and not nck.found

    opt = fp.star_convexity_optimum([1e-11], N_POINTS)
    assert not opt.ok
    assert opt.z_shift_total == 0.0 and opt.g_opt == 0.0


def test_beak_gates_radius_grid_but_not_rho_z_grid() -> None:
    """Gating asymmetry: the plotter path renders what R(theta) rejects."""
    res = fp.radius_grid(BEAK, fp.theta_grid(64), N_POINTS)
    assert res.status == fp.Status.beak_singularity
    assert np.all(res.radii == 0.0)

    profile = fp.rho_z_grid(BEAK, N_POINTS)
    assert profile.ok
    assert profile.rho.max() > 0.0


# --- unchecked cylindrical profile -----------------------------------------

def test_unchecked_profile_draws_a_separated_shape() -> None:
    checked = fp.rho_z_grid(SPLIT, N_POINTS)
    assert checked.status == fp.Status.rho_negative
    assert np.all(checked.rho == 0.0)

    raw = fp.rho_z_grid_unchecked(SPLIT, N_POINTS)
    assert raw.ok
    mid = N_POINTS // 2
    assert raw.rho[mid] == 0.0 and raw.drho_dz[mid] == 0.0      # the void
    assert raw.rho[:mid].max() > 0.1 and raw.rho[mid + 1:].max() > 0.1

    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        assert cache.rho_z_grid(SPLIT).status == fp.Status.rho_negative
        cached = cache.rho_z_grid_unchecked(SPLIT)
    assert cached.ok
    assert _same_bits(cached.z, raw.z)
    assert _same_bits(cached.rho, raw.rho)
    assert _same_bits(cached.drho_dz, raw.drho_dz)
    assert cached.z_shift == raw.z_shift


def test_unchecked_profile_equals_checked_on_a_connected_shape() -> None:
    checked = fp.rho_z_grid(ASYMMETRIC, N_POINTS)
    raw = fp.rho_z_grid_unchecked(ASYMMETRIC, N_POINTS)
    assert checked.ok and raw.ok
    assert _same_bits(raw.z, checked.z)
    assert _same_bits(raw.rho, checked.rho)
    assert _same_bits(raw.drho_dz, checked.drho_dz)


def test_unchecked_profile_keeps_the_c_gate() -> None:
    res = fp.rho_z_grid_unchecked([1e-11], N_POINTS)
    assert res.status == fp.Status.invalid_c
    assert np.all(res.rho == 0.0) and res.z_shift == 0.0


# --- cached tier -----------------------------------------------------------

def test_same_bits_distinguishes_signed_zero() -> None:
    # The hole _same_bits closes: value equality calls +0.0 and -0.0 equal.
    assert np.array_equal([0.0], [-0.0])
    assert not _same_bits(np.array([0.0]), np.array([-0.0]))
    assert _same_bits(np.array([0.0, 1.5]), np.array([0.0, 1.5]))


def test_cache_matches_tier1_bitwise() -> None:
    thetas = fp.theta_grid(64)
    flat = fp.radius_grid(ASYMMETRIC, thetas, N_POINTS)
    with fp.Cache(7, N_POINTS, thetas) as cache:
        cached = cache.radius_grid(ASYMMETRIC)
    assert flat.ok and cached.ok
    assert _same_bits(flat.radii, cached.radii)


def test_cache_reuse_is_stable_across_shapes() -> None:
    thetas = fp.theta_grid(32)
    with fp.Cache(3, N_POINTS, thetas) as cache:
        first = cache.radius_grid(SPHERE)
        cache.radius_grid([1.5, 0.1, 0.2])
        again = cache.radius_grid(SPHERE)
    assert _same_bits(first.radii, again.radii)


def test_cache_derivative_and_at_thetas_agree() -> None:
    thetas = fp.theta_grid(32)
    with fp.Cache(3, N_POINTS, thetas) as cache:
        grid = cache.radius_and_derivative([1.5, 0.08, 0.05])
        at = cache.radius_and_derivative_at_thetas([1.5, 0.08, 0.05], thetas)
    assert grid.ok and at.ok
    # Exact, on every compiler: both outputs call one solve loop that no tier
    # routine can inline. This is the wheel's bitwise gate (GCC 10).
    assert _same_bits(at.radii, grid.radii)
    assert _same_bits(at.dr_dtheta, grid.dr_dtheta)


def test_cache_of_larger_max_params_gives_the_same_bits() -> None:
    thetas = fp.theta_grid(32)
    with fp.Cache(7, N_POINTS, thetas) as exact, \
            fp.Cache(fp.MAX_PARAMS, N_POINTS, thetas) as wide:
        a = exact.radius_and_derivative(ASYMMETRIC)
        b = wide.radius_and_derivative(ASYMMETRIC)
    assert a.ok and b.ok
    assert _same_bits(a.radii, b.radii)
    assert _same_bits(a.dr_dtheta, b.dr_dtheta)


def test_short_vector_equals_its_zero_padded_form() -> None:
    thetas = fp.theta_grid(32)
    short = [1.5, 0.08, 0.05]
    with fp.Cache(8, N_POINTS, thetas) as cache:
        a = cache.radius_and_derivative(short)
        b = cache.radius_and_derivative(short + [0.0] * 5)
        c = cache.radius_and_derivative(short + [-0.0] * 2)
        shape_a = cache.shape(short)
        shape_b = cache.shape(short + [0.0] * 5)
    assert a.ok and b.ok and c.ok
    assert _same_bits(a.radii, b.radii) and _same_bits(a.radii, c.radii)
    assert _same_bits(a.dr_dtheta, b.dr_dtheta)
    assert _same_bits(a.dr_dtheta, c.dr_dtheta)
    assert shape_a == shape_b

    flat = fp.radius_and_derivative(short, thetas, N_POINTS)
    flat_padded = fp.radius_and_derivative(short + [0.0] * 5, thetas, N_POINTS)
    assert _same_bits(flat.radii, a.radii)
    assert _same_bits(flat_padded.radii, a.radii)
    assert _same_bits(flat_padded.dr_dtheta, a.dr_dtheta)


def test_cache_is_stateless_across_a_rejection() -> None:
    thetas = fp.theta_grid(32)
    with fp.Cache(7, N_POINTS, thetas) as cache:
        first = cache.radius_and_derivative(ASYMMETRIC)
        assert cache.radius_grid(BEAK).status == fp.Status.beak_singularity
        assert cache.radius_grid([1e-11]).status == fp.Status.invalid_c
        again = cache.radius_and_derivative(ASYMMETRIC)
    assert first.ok and again.ok
    assert _same_bits(first.radii, again.radii)
    assert _same_bits(first.dr_dtheta, again.dr_dtheta)


def test_cache_shared_between_threads_reproduces_serial_results() -> None:
    """Contract family 3: many threads computing on ONE cache."""
    thetas = fp.theta_grid(48)
    n_threads, n_steps = 8, 40

    def shape_for(tid: int, step: int) -> list:
        t = 0.01 * (tid + 1)
        u = 0.005 * (step % 41)
        full = [1.30 + u + t, 0.05 * t, 0.10 - 0.5 * u, 0.02, 0.01, 0.005, 0.002]
        return full[:3 + 2 * (step % 3)]

    def walk(cache: "fp.Cache", tid: int) -> list:
        out = []
        for step in range(n_steps):
            p = shape_for(tid, step)
            rd = cache.radius_and_derivative(p)
            prof = cache.rho_z_grid(p)
            shp = cache.shape(p)
            out.append((rd.status, rd.radii, rd.dr_dtheta, prof.status, prof.rho,
                        shp.status, shp.z_shift, shp.r_north, shp.r_south))
        return out

    with fp.Cache(7, N_POINTS, thetas) as cache:
        serial = [walk(cache, tid) for tid in range(n_threads)]
        threaded: list = [None] * n_threads

        def worker(tid: int) -> None:
            threaded[tid] = walk(cache, tid)

        pool = [threading.Thread(target=worker, args=(tid,)) for tid in range(n_threads)]
        for th in pool:
            th.start()
        for th in pool:
            th.join()

    n_valid = 0
    for tid in range(n_threads):
        assert threaded[tid] is not None
        for a, b in zip(threaded[tid], serial[tid]):
            assert a[0] == b[0] and a[3] == b[3] and a[5] == b[5]
            assert _same_bits(a[1], b[1]) and _same_bits(a[2], b[2])
            assert _same_bits(a[4], b[4])
            assert a[6:] == b[6:]
            n_valid += int(b[0] == fp.Status.valid and b[3] == fp.Status.valid
                           and b[5] == fp.Status.valid)
    # Zero-filled rejections would compare equal too: the test must not be
    # able to pass vacuously.
    assert n_valid == n_threads * n_steps


def test_cache_shape_rho_neck_and_star_convexity() -> None:
    with fp.Cache(7, N_POINTS, fp.theta_grid(32)) as cache:
        shp = cache.shape(F5_PARAMS)
        prof = cache.rho_z_grid(F5_PARAMS)
        nck = cache.neck(F5_PARAMS)
        opt = cache.star_convexity_optimum(F5_PARAMS)
    assert shp.ok and prof.ok and nck.ok and opt.ok
    assert shp.z_shift == pytest.approx(F5_Z_SHIFT, rel=1e-14)
    assert prof.z.size == N_POINTS
    assert nck.found
    assert opt.z_shift_total == pytest.approx(shp.z_shift, abs=1e-14)


def test_cache_wrong_param_count_is_a_status() -> None:
    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        short = cache.radius_grid([1.0, 0.0])            # 1 .. max_params: accepted
        long = cache.radius_grid([1.0, 0.0, 0.0, 0.01])  # max_params + 1
        empty = cache.radius_grid([])
    assert short.ok
    assert long.status == fp.Status.wrong_param_count
    assert not long.ok
    assert np.all(long.radii == 0.0)
    assert empty.status == fp.Status.wrong_param_count


def test_cache_beak_gating_asymmetry() -> None:
    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        res = cache.radius_grid(BEAK)
        prof = cache.rho_z_grid(BEAK)
    assert res.status == fp.Status.beak_singularity
    assert prof.ok and prof.rho.max() > 0.0


# --- usage errors ----------------------------------------------------------

@pytest.mark.parametrize("max_params,n_points,thetas,code", [
    (51, N_POINTS, fp.theta_grid(16), fp.Status.too_many_params),  # above MAX_PARAMS
    (0, N_POINTS, fp.theta_grid(16), fp.Status.invalid_init),      # below 1
    (3, 50, fp.theta_grid(16), fp.Status.invalid_grid),            # below the n_points floor
    (3, N_POINTS, [], fp.Status.invalid_grid),                     # no thetas
    (3, N_POINTS, [0.5, 4.0], fp.Status.invalid_grid),             # theta outside [0, pi]
    (0, 50, [], fp.Status.invalid_init),                           # count is judged first
])
def test_bad_cache_arguments_raise_with_the_status(max_params, n_points, thetas,
                                                   code) -> None:
    with pytest.raises(fp.FosParamError) as excinfo:
        fp.Cache(max_params, n_points, thetas)
    assert excinfo.value.status == code
    assert code.name in str(excinfo.value)


def test_cache_accepts_both_limits_of_max_params() -> None:
    for max_params in (1, fp.MAX_PARAMS):
        with fp.Cache(max_params, N_POINTS, fp.theta_grid(16)) as cache:
            assert cache.max_params == max_params
            assert cache.n_points == N_POINTS and cache.n_thetas == 16
            assert cache.radius_grid(SPHERE[:1]).ok


def test_python_side_errors_carry_no_status() -> None:
    with pytest.raises(fp.FosParamError) as excinfo:
        fp.radius_grid([[1.0, 0.0]], fp.theta_grid(16), N_POINTS)
    assert excinfo.value.status is None


def test_use_after_close_raises() -> None:
    cache = fp.Cache(3, N_POINTS, fp.theta_grid(16))
    cache.close()
    cache.close()   # idempotent
    with pytest.raises(fp.FosParamError, match="closed"):
        cache.radius_grid(SPHERE)


def test_context_manager_closes() -> None:
    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        assert cache.radius_grid(SPHERE).ok
    with pytest.raises(fp.FosParamError, match="closed"):
        cache.radius_grid(SPHERE)


def test_non_1d_input_raises() -> None:
    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        with pytest.raises(fp.FosParamError):
            cache.radius_grid(np.zeros((2, 3)))
    with pytest.raises(fp.FosParamError):
        fp.radius_grid(np.zeros((2, 3)), fp.theta_grid(16), N_POINTS)
    with pytest.raises(fp.FosParamError):
        fp.radius_grid(SPHERE, np.zeros((2, 3)), N_POINTS)


# --- diagnostics -----------------------------------------------------------

def test_status_enum_covers_the_c_contract() -> None:
    assert [int(s) for s in fp.Status] == [0, 1, 2, 3, 4, 5,
                                           100, 101, 102, 103, 104, 105]
    assert "unknown" in fp.status_message(6)      # retired, never reused


def test_status_message_covers_every_code() -> None:
    for status in fp.Status:
        message = fp.status_message(status)
        assert message and "unknown" not in message
    assert "unknown" in fp.status_message(999)


def test_limits_are_exported() -> None:
    assert fp.MAX_PARAMS == 50
    assert fp.N_POINTS_FLOOR == 100


# --- removed 1.x and 2.x surface -------------------------------------------

@pytest.mark.parametrize("name", [
    "rho_profile", "RhoProfileResult", "ShapeResult", "MESSAGE_BUFFER_SIZE",
    "CACHE_MAX_PARAMS",
])
def test_removed_names_absent(name: str) -> None:
    assert not hasattr(fp, name)
    assert name not in fp.__all__


# --- inputs ----------------------------------------------------------------

def test_inputs_are_converted_to_contiguous_float64() -> None:
    """Lists of ints, float32 arrays and strided views give the float64 result."""
    thetas = fp.theta_grid(16)
    ref = fp.radius_grid([2.0, 0.0, 0.25], thetas, N_POINTS)   # exact in float32
    assert ref.ok
    variants = [
        fp.radius_grid([2, 0, 0.25], thetas, N_POINTS),
        fp.radius_grid(np.array([2.0, 0.0, 0.25], dtype=np.float32), thetas, N_POINTS),
        fp.radius_grid(np.array([2.0, 9.0, 0.0, 9.0, 0.25])[::2], thetas, N_POINTS),
        fp.radius_grid((2.0, 0.0, 0.25), list(thetas), N_POINTS),
        fp.radius_grid([2.0, 0.0, 0.25], thetas[::-1][::-1], N_POINTS),
    ]
    for res in variants:
        assert res.ok
        assert _same_bits(res.radii, ref.radii)


def test_cache_copies_its_thetas() -> None:
    """Mutating the caller's theta array after construction moves no result."""
    thetas = fp.theta_grid(16)
    with fp.Cache(7, N_POINTS, thetas) as cache:
        before = cache.radius_grid(ASYMMETRIC)
        thetas[:] = 0.3
        after = cache.radius_grid(ASYMMETRIC)
    assert before.ok and after.ok
    assert _same_bits(before.radii, after.radii)


def test_removed_2x_members_absent() -> None:
    assert not hasattr(fp.Status, "tables_not_initialized")
    with fp.Cache(3, N_POINTS, fp.theta_grid(16)) as cache:
        assert not hasattr(cache, "n_params")
