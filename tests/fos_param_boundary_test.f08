!> Contract family 4: boundaries. Every declared limit is asserted on BOTH
!! sides — the last accepted value and the first rejected one — so an off-by-one
!! in either direction fails here.
!!
!!   parameter count, cached tier   8 -> 0,  9 -> SHAPE_ERROR_TOO_MANY_PARAMS
!!   parameter count, standalone   50 -> ok, 51 -> SHAPE_ERROR_TOO_MANY_PARAMS
!!   u-grid resolution            100 -> 0, 99 -> SHAPE_ERROR_INVALID_GRID
!!   theta range                    pi -> 0, 3pi/2 -> SHAPE_ERROR_INVALID_GRID
!!   theta count                     1 -> 0,  0 -> SHAPE_ERROR_INVALID_GRID
!!   vector length                  n_params -> 0, anything else -> 4
!!
!! The two caps are different numbers on purpose. The cached tier's trig tables
!! are built once, at a fixed k_max = (8+2)/2+1, so 8 is its hard ceiling; the
!! standalone forms size their tables per call, which lifts them to
!! FOS_MAX_PARAMS = 50.
!!
!! "ok" for a 50-parameter standalone call means "not rejected on LENGTH": what
!! the shape gates then decide about that particular shape is a different
!! question, and asserting a specific verdict would pin geometry, not the cap.
program fos_param_boundary_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_at_thetas_s, &
            tables_t, tables_init_s, tables_free_s, &
            compute_shape_standalone_s, compute_radius_grid_standalone_s, &
            FOS_MAX_PARAMS, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 8_ik
    integer(kind = ik), parameter :: N_POINTS_FLOOR = 100_ik

    real(kind = rk), parameter :: BASE8(8) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    type(cache_t)      :: cache
    !! Second handle for the theta-count probes: `cache` stays live across them.
    type(cache_t)      :: cache_probe
    type(tables_t)     :: tables
    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: one_theta(1), no_thetas(0)
    real(kind = rk)    :: bad_thetas(N_THETA)
    real(kind = rk)    :: params50(50), params51(51)
    real(kind = rk)    :: radii(N_THETA), dr_dtheta(N_THETA)
    real(kind = rk)    :: radii50(N_THETA)
    real(kind = rk)    :: z_shift, r_north, r_south
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do
    one_theta(1) = 0.5_rk * PI_C

    !---------------------------------------------------------------------------
    ! Parameter count, cached tier: 8 accepted, 9 rejected
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 8_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached tier: 8 params accepted')
    call cache_free_s(cache)

    call cache_init_s(cache, 9_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
            'cached tier: 9 params -> 1')
    call cache_free_s(cache)

    !---------------------------------------------------------------------------
    ! Parameter count, standalone tier: 50 accepted, 51 rejected
    !---------------------------------------------------------------------------
    call assert_int_eq(FOS_MAX_PARAMS, 50_ik, 'FOS_MAX_PARAMS = 50')

    do i = 1_ik, 50_ik
        params50(i) = 1.0e-6_rk
    end do
    params50(1) = 1.5_rk
    do i = 1_ik, 50_ik
        params51(i) = params50(i)
    end do
    params51(51) = 1.0e-6_rk

    call compute_shape_standalone_s(params50, N_POINTS, z_shift, r_north, &
            r_south, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            'standalone: 50 params not rejected on length (shape)')
    call compute_radius_grid_standalone_s(params50, thetas, N_POINTS, radii50, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            'standalone: 50 params not rejected on length (radius grid)')

    call compute_shape_standalone_s(params51, N_POINTS, z_shift, r_north, &
            r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
            'standalone: 51 params -> 1 (shape)')
    call compute_radius_grid_standalone_s(params51, thetas, N_POINTS, radii50, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
            'standalone: 51 params -> 1 (radius grid)')

    !---------------------------------------------------------------------------
    ! u-grid resolution: the floor is 100
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 8_ik, N_POINTS_FLOOR, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'n_points 100 accepted')
    call cache_free_s(cache)

    call cache_init_s(cache, 8_ik, N_POINTS_FLOOR - 1_ik, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'n_points 99 -> 3')
    call cache_free_s(cache)

    call tables_init_s(tables, N_POINTS_FLOOR, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'tables: n_points 100 accepted')
    call tables_free_s(tables)

    call tables_init_s(tables, N_POINTS_FLOOR - 1_ik, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'tables: n_points 99 -> 3')
    call tables_free_s(tables)

    call compute_shape_standalone_s(BASE8, N_POINTS_FLOOR, z_shift, r_north, &
            r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone: n_points 100 accepted')
    call compute_shape_standalone_s(BASE8, N_POINTS_FLOOR - 1_ik, z_shift, &
            r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'standalone: n_points 99 -> 3')

    !---------------------------------------------------------------------------
    ! Theta range: pi is inside, 3pi/2 is not — at tables init AND at at_thetas
    !---------------------------------------------------------------------------
    do i = 1_ik, N_THETA
        bad_thetas(i) = thetas(i)
    end do
    bad_thetas(N_THETA) = PI_C

    call tables_init_s(tables, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'tables: theta = pi accepted')
    call tables_free_s(tables)

    bad_thetas(N_THETA) = 1.5_rk * PI_C
    call tables_init_s(tables, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'tables: theta = 3pi/2 -> 3')
    call tables_free_s(tables)

    call cache_init_s(cache, 8_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'at-thetas: cache init valid')

    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, thetas, radii, &
            dr_dtheta, status)
    call assert_int_eq(status, SHAPE_VALID, 'at_thetas: in-range nodes accepted')

    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, bad_thetas, radii, &
            dr_dtheta, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, &
            'at_thetas: theta = 3pi/2 -> 3')
    call assert_true(all_zero_f(radii) .and. all_zero_f(dr_dtheta), &
            'at_thetas: out-of-range node zero-fills')

    !---------------------------------------------------------------------------
    ! Theta count: one node is enough, zero is not — at every init path that
    ! carries thetas (tables, cache, and the theta-bearing standalone forms).
    ! The theta-LESS standalone forms are exempt by design and not asserted here.
    !---------------------------------------------------------------------------
    call tables_init_s(tables, N_POINTS, one_theta, status)
    call assert_int_eq(status, SHAPE_VALID, 'tables: one theta node accepted')
    call tables_free_s(tables)

    call tables_init_s(tables, N_POINTS, no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'tables: n_theta 0 -> 3')
    call tables_free_s(tables)

    call cache_init_s(cache_probe, 8_ik, N_POINTS, one_theta, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache: one theta node accepted')
    call cache_free_s(cache_probe)

    call cache_init_s(cache_probe, 8_ik, N_POINTS, no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'cache: n_theta 0 -> 3')
    call cache_free_s(cache_probe)

    call compute_radius_grid_standalone_s(BASE8, no_thetas, N_POINTS, &
            radii50(1:0), status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, &
            'standalone radius grid: n_theta 0 -> 3')

    !---------------------------------------------------------------------------
    ! Vector length: exactly n_params, nothing else
    !---------------------------------------------------------------------------
    call cache_radius_grid_s(cache, BASE8, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'exactly 8 params accepted')

    call cache_radius_grid_s(cache, BASE8(1:7), radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, '7 params -> 4')
    call assert_true(all_zero_f(radii), 'short vector zero-fills')

    call cache_radius_grid_s(cache, params50(1:9), radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, '9 params -> 4')

    call cache_radius_grid_s(cache, BASE8(1:0), radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'empty vector -> 4')

    call cache_free_s(cache)

    call test_summary()

contains

    !> .true. iff every element is exactly zero.
    pure function all_zero_f(a) result(zeroed)
        real(kind = rk), intent(in) :: a(:)
        logical :: zeroed
        integer(kind = ik) :: j
        zeroed = .true.
        do j = 1_ik, size(a, kind = ik)
            if (abs(a(j)) > 0.0_rk) zeroed = .false.
        end do
    end function all_zero_f

end program fos_param_boundary_test
