!> Contract family 4: boundaries. Every declared limit is asserted on BOTH
!! sides — the last accepted value and the first rejected one — so an
!! off-by-one in either direction fails here. With L = FOS_MAX_PARAMS = 50:
!!
!!   cache init, max_params        50 -> ok, 51 -> 1, 0 -> 5
!!   cached call, size(params)     max_params -> ok, max_params + 1 -> 4,
!!                                 empty -> 4
!!   one-shot, size(params)        50 -> ok, 51 -> 1, empty -> 4
!!   u-grid resolution            100 -> ok, 99 -> 3 (init and one-shot)
!!   theta count                    1 -> ok,  0 -> 3 (init, one-shot radius
!!                                 forms, at-thetas)
!!   theta range                   pi -> ok, 3pi/2 -> 3
!!
!! At each accepted maximum the LAST parameter must change the output: a
!! library that accepted a long vector and silently ignored its tail would pass
!! every status assertion above.
program fos_param_boundary_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_at_thetas_s, &
            cache_rho_z_grid_s, cache_shape_s, &
            compute_shape_standalone_s, compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            compute_rho_z_grid_standalone_s, &
            FOS_MAX_PARAMS, FOS_N_POINTS_FLOOR, SHAPE_MAX_PARAMS, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT
    use test_utils_mod, only: assert_true, assert_int_eq, arrays_bits_eq_f, &
            all_zero_f, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 8_ik
    integer(kind = ik), parameter :: L = 50_ik

    real(kind = rk), parameter :: BASE8(8) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    type(cache_t)      :: cache
    real(kind = rk)    :: thetas(N_THETA), bad_thetas(N_THETA)
    real(kind = rk)    :: one_theta(1), no_thetas(0), no_params(0)
    real(kind = rk)    :: params50(50), params51(51), nudged50(50)
    real(kind = rk)    :: params9(9), nudged8(8)
    real(kind = rk)    :: radii(N_THETA), dr_dtheta(N_THETA), radii_b(N_THETA)
    real(kind = rk)    :: one_r(1), one_d(1), no_r(0), no_d(0)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho(N_POINTS)
    real(kind = rk)    :: rho_b(N_POINTS), z_floor(100), rho_floor(100), drho_floor(100)
    real(kind = rk)    :: z_shift, r_north, r_south
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do
    one_theta(1) = 0.5_rk * PI_C

    ! Valid long vectors: every slot nonzero, amplitudes small enough that the
    ! high orders do not fold the surface.
    do i = 1_ik, 51_ik
        params51(i) = 1.0e-6_rk
    end do
    params51(1) = 1.5_rk
    do i = 1_ik, 50_ik
        params50(i) = params51(i)
        nudged50(i) = params51(i)
    end do
    nudged50(50) = 2.0e-6_rk
    do i = 1_ik, 8_ik
        params9(i) = BASE8(i)
        nudged8(i) = BASE8(i)
    end do
    params9(9) = 1.0e-3_rk
    nudged8(8) = BASE8(8) + 1.0e-3_rk

    call assert_int_eq(FOS_MAX_PARAMS, L, 'FOS_MAX_PARAMS = 50')
    call assert_int_eq(min(SHAPE_MAX_PARAMS, FOS_MAX_PARAMS), L, &
            'L = min(SHAPE_MAX_PARAMS, FOS_MAX_PARAMS) = 50')

    !---------------------------------------------------------------------------
    ! Cache init: max_params = L accepted, L + 1 -> 1, 0 -> 5
    !---------------------------------------------------------------------------
    call cache_init_s(cache, L + 1_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'init: max_params 51 -> 1')
    call cache_init_s(cache, 0_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init: max_params 0 -> 5')
    call cache_init_s(cache, L, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: max_params 50 accepted')

    !---------------------------------------------------------------------------
    ! Cached call on the L cache: 50 parameters accepted and all of them used;
    ! 51 -> 4; empty -> 4
    !---------------------------------------------------------------------------
    call cache_rho_z_grid_s(cache, params50, z, rho, drho, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: 50 params accepted')
    call cache_rho_z_grid_s(cache, nudged50, z, rho_b, drho, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: nudged 50 params accepted')
    call assert_true(.not. arrays_bits_eq_f(rho, rho_b), &
            'cached: the 50th parameter changes the output')

    call cache_radius_grid_s(cache, params50, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: 50 params give an R(theta) grid')
    call cache_radius_grid_s(cache, nudged50, radii_b, status)
    call assert_true(.not. arrays_bits_eq_f(radii, radii_b), &
            'cached: the 50th parameter changes R(theta)')

    radii = 1.0_rk
    call cache_radius_grid_s(cache, params51, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'cached: 51 params -> 4')
    call assert_true(all_zero_f(radii), 'cached: 4 zero-fills')
    call cache_radius_grid_s(cache, no_params, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'cached: empty params -> 4')
    call cache_free_s(cache)

    !---------------------------------------------------------------------------
    ! Cached call on a max_params = 8 cache: 8 accepted and the 8th used; 9 -> 4
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 8_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: max_params 8 accepted')
    call cache_radius_grid_s(cache, BASE8, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: size = max_params accepted')
    call cache_radius_grid_s(cache, nudged8, radii_b, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: nudged 8 params accepted')
    call assert_true(.not. arrays_bits_eq_f(radii, radii_b), &
            'cached: the last parameter of a full vector changes R(theta)')
    call cache_radius_grid_s(cache, params9, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'cached: max_params + 1 -> 4')
    call cache_shape_s(cache, params9, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            'cached: max_params + 1 -> 4 (shape)')

    !---------------------------------------------------------------------------
    ! Theta range and count on the at-thetas output
    !---------------------------------------------------------------------------
    do i = 1_ik, N_THETA
        bad_thetas(i) = thetas(i)
    end do
    bad_thetas(N_THETA) = PI_C
    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, bad_thetas, radii, &
            dr_dtheta, status)
    call assert_int_eq(status, SHAPE_VALID, 'at-thetas: theta = pi accepted')
    bad_thetas(N_THETA) = 1.5_rk * PI_C
    radii = 1.0_rk
    dr_dtheta = 1.0_rk
    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, bad_thetas, radii, &
            dr_dtheta, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'at-thetas: theta = 3pi/2 -> 3')
    call assert_true(all_zero_f(radii) .and. all_zero_f(dr_dtheta), &
            'at-thetas: 3 zero-fills')

    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, one_theta, one_r, &
            one_d, status)
    call assert_int_eq(status, SHAPE_VALID, 'at-thetas: one theta accepted')
    call cache_radius_and_derivative_at_thetas_s(cache, BASE8, no_thetas, no_r, &
            no_d, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'at-thetas: empty theta set -> 3')
    call cache_free_s(cache)

    !---------------------------------------------------------------------------
    ! One-shot: 50 accepted and all of them used; 51 -> 1; empty -> 4
    !---------------------------------------------------------------------------
    call compute_rho_z_grid_standalone_s(params50, N_POINTS, z, rho, drho, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'one-shot: 50 params accepted')
    call compute_rho_z_grid_standalone_s(nudged50, N_POINTS, z, rho_b, drho, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'one-shot: nudged 50 params accepted')
    call assert_true(.not. arrays_bits_eq_f(rho, rho_b), &
            'one-shot: the 50th parameter changes the output')

    call compute_shape_standalone_s(params51, N_POINTS, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'one-shot: 51 params -> 1 (shape)')
    radii = 1.0_rk
    call compute_radius_grid_standalone_s(params51, thetas, N_POINTS, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, &
            'one-shot: 51 params -> 1 (radius grid)')
    call assert_true(all_zero_f(radii), 'one-shot: 1 zero-fills')

    call compute_shape_standalone_s(no_params, N_POINTS, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'one-shot: empty params -> 4 (shape)')
    call compute_radius_grid_standalone_s(no_params, thetas, N_POINTS, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            'one-shot: empty params -> 4 (radius grid)')

    !---------------------------------------------------------------------------
    ! u-grid resolution: the floor is 100
    !---------------------------------------------------------------------------
    call assert_int_eq(FOS_N_POINTS_FLOOR, 100_ik, 'FOS_N_POINTS_FLOOR = 100')

    call cache_init_s(cache, 8_ik, FOS_N_POINTS_FLOOR, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: n_points 100 accepted')
    call cache_rho_z_grid_s(cache, BASE8, z_floor, rho_floor, drho_floor, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: the floor grid computes')
    call cache_free_s(cache)
    call cache_init_s(cache, 8_ik, FOS_N_POINTS_FLOOR - 1_ik, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: n_points 99 -> 3')

    call compute_shape_standalone_s(BASE8, FOS_N_POINTS_FLOOR, z_shift, r_north, &
            r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'one-shot: n_points 100 accepted')
    call compute_shape_standalone_s(BASE8, FOS_N_POINTS_FLOOR - 1_ik, z_shift, &
            r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'one-shot: n_points 99 -> 3')

    !---------------------------------------------------------------------------
    ! Theta count: one node is enough, zero is not — at init and in the
    ! one-shot radius forms. The theta-less one-shot forms take no thetas.
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 8_ik, N_POINTS, one_theta, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: one theta accepted')
    call cache_radius_grid_s(cache, BASE8, one_r, status)
    call assert_int_eq(status, SHAPE_VALID, 'cached: a one-theta cache computes')
    call cache_free_s(cache)
    call cache_init_s(cache, 8_ik, N_POINTS, no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: empty theta set -> 3')

    call compute_radius_grid_standalone_s(BASE8, one_theta, N_POINTS, one_r, status)
    call assert_int_eq(status, SHAPE_VALID, 'one-shot: one theta accepted')
    call compute_radius_grid_standalone_s(BASE8, no_thetas, N_POINTS, no_r, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, &
            'one-shot: empty theta set -> 3 (radius grid)')
    call compute_radius_and_derivative_standalone_s(BASE8, no_thetas, N_POINTS, &
            no_r, no_d, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, &
            'one-shot: empty theta set -> 3 (radius and derivative)')

    !---------------------------------------------------------------------------
    ! Theta range at init and in the one-shot: pi inside, 3pi/2 not
    !---------------------------------------------------------------------------
    bad_thetas(N_THETA) = PI_C
    call cache_init_s(cache, 8_ik, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: theta = pi accepted')
    call cache_free_s(cache)
    bad_thetas(N_THETA) = 1.5_rk * PI_C
    call cache_init_s(cache, 8_ik, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: theta = 3pi/2 -> 3')
    call compute_radius_grid_standalone_s(BASE8, bad_thetas, N_POINTS, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'one-shot: theta = 3pi/2 -> 3')

    call test_summary()

end program fos_param_boundary_test
