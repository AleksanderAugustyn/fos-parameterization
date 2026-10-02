!> Lifecycle suite of the read-only cache: init and every init rejection in its
!! documented order, the accessors, re-init without a free, copy-assignment,
!! free, and the order of the usage checks on a compute call.
!!
!! A cache holds nothing derived from a shape parameter, so there is no
!! per-shape state to test here: what a compute returns is the business of the
!! `outputs`, `resolve` and `extensions` suites, and what it does NOT keep is
!! the business of `statelessness`.
program fos_param_cache_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_max_params_f, cache_n_points_f, cache_n_thetas_f, &
            cache_is_initialized_f, cache_radius_grid_s, cache_rho_z_grid_s, &
            FOS_MAX_PARAMS, FOS_N_POINTS_FLOOR, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            FOS_ERROR_BUFFER_MISMATCH, FOS_ERROR_INVALID_C
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            all_zero_f, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 301_ik
    integer(kind = ik), parameter :: N_THETA = 16_ik
    integer(kind = ik), parameter :: MAX_PARAMS = 8_ik

    real(kind = rk), parameter :: PARAMS3(3) = [1.5_rk, 0.1_rk, 0.05_rk]
    real(kind = rk), parameter :: BAD_C3(3) = [1.0e-11_rk, 0.1_rk, 0.05_rk]
    real(kind = rk), parameter :: LONG9(MAX_PARAMS + 1_ik) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk, &
             0.001_rk, 0.0005_rk]

    type(cache_t)      :: cache, copy, never
    real(kind = rk)    :: thetas(N_THETA), bad_thetas(N_THETA), few_thetas(5)
    real(kind = rk)    :: no_thetas(0), no_params(0)
    real(kind = rk)    :: radii(N_THETA), radii_copy(N_THETA), radii_few(5)
    real(kind = rk)    :: short_radii(N_THETA - 1_ik)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS), z_shift
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do
    do i = 1_ik, 5_ik
        few_thetas(i) = real(i, rk) * PI_C / 6.0_rk
    end do

    !---------------------------------------------------------------------------
    ! A cache that was never initialized
    !---------------------------------------------------------------------------
    call assert_true(.not. cache_is_initialized_f(never), 'fresh cache: not initialized')
    call assert_int_eq(cache_max_params_f(never), 0_ik, 'fresh cache: max_params 0')
    call assert_int_eq(cache_n_points_f(never), 0_ik, 'fresh cache: n_points 0')
    call assert_int_eq(cache_n_thetas_f(never), 0_ik, 'fresh cache: n_thetas 0')

    radii = 1.0_rk
    call cache_radius_grid_s(never, PARAMS3, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
            'fresh cache: compute -> 2')
    call assert_true(all_zero_f(radii), 'fresh cache: compute zero-fills')

    call cache_free_s(never)
    call assert_true(.not. cache_is_initialized_f(never), &
            'freeing a never-initialized cache is harmless')

    !---------------------------------------------------------------------------
    ! Init rejections, in the documented order: 5, 1, then the grid (3)
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 0_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init: max_params 0 -> 5')
    call cache_init_s(cache, -3_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init: max_params -3 -> 5')
    call cache_init_s(cache, FOS_MAX_PARAMS + 1_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'init: max_params 51 -> 1')
    call assert_true(.not. cache_is_initialized_f(cache), &
            'init: a rejected init leaves the cache uninitialized')

    call cache_init_s(cache, MAX_PARAMS, FOS_N_POINTS_FLOOR - 1_ik, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: n_points 99 -> 3')
    call cache_init_s(cache, MAX_PARAMS, N_POINTS, no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: empty theta set -> 3')

    bad_thetas = thetas
    bad_thetas(N_THETA) = 1.5_rk * PI_C
    call cache_init_s(cache, MAX_PARAMS, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: theta = 3pi/2 -> 3')
    bad_thetas = thetas
    bad_thetas(1) = -1.0e-6_rk
    call cache_init_s(cache, MAX_PARAMS, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init: negative theta -> 3')
    call assert_true(.not. cache_is_initialized_f(cache), &
            'init: a rejected grid leaves the cache uninitialized')

    ! The parameter count is judged before the grid
    call cache_init_s(cache, 0_ik, FOS_N_POINTS_FLOOR - 1_ik, no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init: 5 outranks 3')
    call cache_init_s(cache, FOS_MAX_PARAMS + 1_ik, FOS_N_POINTS_FLOOR - 1_ik, &
            no_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'init: 1 outranks 3')

    !---------------------------------------------------------------------------
    ! A successful init and its accessors
    !---------------------------------------------------------------------------
    call cache_init_s(cache, MAX_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: valid arguments accepted')
    call assert_true(cache_is_initialized_f(cache), 'live cache: initialized')
    call assert_int_eq(cache_max_params_f(cache), MAX_PARAMS, 'live cache: max_params')
    call assert_int_eq(cache_n_points_f(cache), N_POINTS, 'live cache: n_points')
    call assert_int_eq(cache_n_thetas_f(cache), N_THETA, 'live cache: n_thetas')

    ! Both limits of max_params are accepted
    call cache_init_s(copy, 1_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: max_params 1 accepted')
    call cache_init_s(copy, FOS_MAX_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init: max_params 50 accepted')
    call assert_int_eq(cache_max_params_f(copy), FOS_MAX_PARAMS, &
            're-init without a free: max_params follows the new init')
    call cache_free_s(copy)

    !---------------------------------------------------------------------------
    ! Usage checks on a compute, in order: 4, then 105, then the value codes
    !---------------------------------------------------------------------------
    call cache_radius_grid_s(cache, no_params, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'compute: empty params -> 4')
    call cache_radius_grid_s(cache, LONG9, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            'compute: max_params + 1 -> 4')

    short_radii = 1.0_rk
    call cache_radius_grid_s(cache, PARAMS3, short_radii, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'compute: wrong buffer -> 105')
    call assert_true(all_zero_f(short_radii), 'compute: 105 zero-fills')

    call cache_radius_grid_s(cache, LONG9, short_radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'compute: 4 outranks 105')

    call cache_radius_grid_s(cache, BAD_C3, short_radii, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'compute: 105 outranks 102')
    call cache_radius_grid_s(cache, BAD_C3, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'compute: degenerate c -> 102')

    z = 1.0_rk
    rho = 1.0_rk
    drho_dz = 1.0_rk
    call cache_rho_z_grid_s(cache, BAD_C3, z(1:N_POINTS - 1_ik), rho, drho_dz, &
            z_shift, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, &
            'compute: wrong grid buffer -> 105, ahead of 102')
    call assert_true(all_zero_f(rho) .and. all_zero_f(drho_dz), &
            'compute: 105 zero-fills every grid buffer')

    !---------------------------------------------------------------------------
    ! Copy-assignment is a deep copy: the copy survives the original
    !---------------------------------------------------------------------------
    call cache_radius_grid_s(cache, PARAMS3, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'original cache computes')

    ! The cache owns a copy of its thetas: overwriting the caller's array
    ! afterwards does not move a result.
    bad_thetas = thetas
    call cache_init_s(copy, MAX_PARAMS, N_POINTS, bad_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache on a scratch theta array')
    bad_thetas = 0.25_rk * PI_C
    call cache_radius_grid_s(copy, PARAMS3, radii_copy, status)
    call assert_int_eq(status, SHAPE_VALID, 'compute after the caller overwrote its thetas')
    do i = 1_ik, N_THETA
        call assert_bits_eq(radii_copy(i), radii(i), &
                'the cache kept its own thetas, bit for bit')
    end do

    copy = cache
    call cache_free_s(cache)
    call assert_true(.not. cache_is_initialized_f(cache), 'freed cache: not initialized')
    call assert_int_eq(cache_max_params_f(cache), 0_ik, 'freed cache: max_params 0')
    call assert_int_eq(cache_n_points_f(cache), 0_ik, 'freed cache: n_points 0')
    call assert_int_eq(cache_n_thetas_f(cache), 0_ik, 'freed cache: n_thetas 0')

    call cache_radius_grid_s(cache, PARAMS3, radii_copy, status)
    call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, 'freed cache: compute -> 2')
    call cache_free_s(cache)
    call assert_true(.not. cache_is_initialized_f(cache), 'a second free is harmless')

    call assert_true(cache_is_initialized_f(copy), 'the copy outlives the original')
    call cache_radius_grid_s(copy, PARAMS3, radii_copy, status)
    call assert_int_eq(status, SHAPE_VALID, 'the copy computes')
    do i = 1_ik, N_THETA
        call assert_bits_eq(radii_copy(i), radii(i), 'copy == original, bit for bit')
    end do

    !---------------------------------------------------------------------------
    ! Re-init without a free: the old contents are released, the new ones serve
    !---------------------------------------------------------------------------
    call cache_init_s(copy, 3_ik, 201_ik, few_thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 're-init without a free accepted')
    call assert_int_eq(cache_max_params_f(copy), 3_ik, 're-init: max_params')
    call assert_int_eq(cache_n_points_f(copy), 201_ik, 're-init: n_points')
    call assert_int_eq(cache_n_thetas_f(copy), 5_ik, 're-init: n_thetas')
    call cache_radius_grid_s(copy, PARAMS3, radii_few, status)
    call assert_int_eq(status, SHAPE_VALID, 're-init: the new cache computes')
    call cache_radius_grid_s(copy, PARAMS3, radii, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, &
            're-init: the old theta count is now a wrong buffer')

    ! A rejected re-init leaves the cache uninitialized, not half-alive
    call cache_init_s(copy, 3_ik, 50_ik, few_thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'rejected re-init -> 3')
    call assert_true(.not. cache_is_initialized_f(copy), &
            'rejected re-init: cache uninitialized')

    call cache_free_s(copy)
    call test_summary()

end program fos_param_cache_test
