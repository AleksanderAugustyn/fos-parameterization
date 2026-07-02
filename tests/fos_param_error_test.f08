!> Error-path tests: every rejection reports the matching FOS_* code.
program fos_param_error_test

    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: compute_fos_radius_grid_s, &
            FOS_VALID, FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID_SMALL = 101_ik

    real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
    real(kind = rk) :: empty_params(0)
    logical :: is_valid
    character(len = 256) :: message
    integer(kind = ik) :: error_code

    write(*, '(A)') '=== Error-path tests ==='

    ! c = 0
    params = 0.0_rk
    call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
            message, error_code = error_code)
    call assert_true(.not. is_valid, 'c = 0: rejected')
    call assert_int_eq(error_code, FOS_ERROR_INVALID_C, 'c = 0: FOS_ERROR_INVALID_C')

    ! c < 0
    params = 0.0_rk
    params(1) = -1.0_rk
    call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
            message, error_code = error_code)
    call assert_true(.not. is_valid, 'c < 0: rejected')
    call assert_int_eq(error_code, FOS_ERROR_INVALID_C, 'c < 0: FOS_ERROR_INVALID_C')

    ! Empty parameter array
    call compute_fos_radius_grid_s(empty_params, N_GRID_SMALL, radii, z_shift, is_valid, &
            message, error_code = error_code)
    call assert_true(.not. is_valid, 'empty params: rejected')
    call assert_int_eq(error_code, FOS_ERROR_INVALID_C, 'empty params: FOS_ERROR_INVALID_C')

    ! Beak singularity: for a3 = a5 = a6 = 0, a2 = a4/3 and f(0) = 1 - 4 a4 / 3.
    ! a4 = 0.7497 gives f(0) = 4.0e-4: rho(0) > 0 (passes the rho check) but
    ! f_min < F_MIN_THRESHOLD = 1e-3 (fails beak detection).
    params = 0.0_rk
    params(1) = 2.0_rk
    params(3) = 0.7497_rk
    call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
            message, error_code = error_code)
    call assert_true(.not. is_valid, 'near-scission beak: rejected')
    call assert_int_eq(error_code, FOS_ERROR_BEAK_SINGULARITY, &
            'near-scission beak: FOS_ERROR_BEAK_SINGULARITY')

    ! Success path reports FOS_VALID
    params = 0.0_rk
    params(1) = 1.0_rk
    call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
            message, error_code = error_code)
    call assert_true(is_valid, 'sphere: accepted')
    call assert_int_eq(error_code, FOS_VALID, 'sphere: FOS_VALID')

    call test_summary()

end program fos_param_error_test
