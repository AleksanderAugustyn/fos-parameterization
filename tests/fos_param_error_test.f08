!> Error-path tests: every rejection reports the matching FOS_* code.
program fos_param_error_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_radius_grid_standalone_s, &
            status_message, &
            FOS_VALID, FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID_SMALL = 101_ik

    real(kind = rk) :: params(7), radii(N_GRID_SMALL), thetas(N_GRID_SMALL)
    real(kind = rk) :: empty_params(0)
    integer(kind = ik) :: status, i

    write(*, '(A)') '=== Error-path tests ==='

    ! Endpoints pinned: under -ffast-math the last product can land one ulp
    ! above PI_C, which the library rejects as an out-of-domain theta.
    do i = 1_ik, N_GRID_SMALL
        thetas(i) = real(i - 1_ik, rk) * PI_C / real(N_GRID_SMALL - 1_ik, rk)
    end do
    thetas(1) = 0.0_rk
    thetas(N_GRID_SMALL) = PI_C

    ! c = 0
    params = 0.0_rk
    call compute_radius_grid_standalone_s(params, thetas, N_GRID_SMALL, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'c = 0: FOS_ERROR_INVALID_C')
    call assert_true(len_trim(status_message(status)) > 0_ik, 'c = 0: message non-empty')

    ! c < 0
    params = 0.0_rk
    params(1) = -1.0_rk
    call compute_radius_grid_standalone_s(params, thetas, N_GRID_SMALL, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'c < 0: FOS_ERROR_INVALID_C')

    ! Empty parameter array
    call compute_radius_grid_standalone_s(empty_params, thetas, N_GRID_SMALL, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'empty params: FOS_ERROR_INVALID_C')

    ! Beak singularity: for a3 = a5 = a6 = 0, a2 = a4/3 and f(0) = 1 - 4 a4 / 3.
    ! a4 = 0.7497 gives f(0) = 4.0e-4: rho(0) > 0 (passes the rho check) but
    ! f_min < F_MIN_THRESHOLD = 1e-3 (fails beak detection).
    params = 0.0_rk
    params(1) = 2.0_rk
    params(3) = 0.7497_rk
    call compute_radius_grid_standalone_s(params, thetas, N_GRID_SMALL, radii, status)
    call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, &
            'near-scission beak: FOS_ERROR_BEAK_SINGULARITY')
    call assert_true(len_trim(status_message(status)) > 0_ik, 'beak: message non-empty')

    ! Every rejection zero-fills the whole output
    call assert_true(all(abs(radii) <= 0.0_rk), 'rejection zero-fills radii')

    ! Success path reports FOS_VALID
    params = 0.0_rk
    params(1) = 1.0_rk
    call compute_radius_grid_standalone_s(params, thetas, N_GRID_SMALL, radii, status)
    call assert_int_eq(status, FOS_VALID, 'sphere: FOS_VALID')

    call test_summary()

end program fos_param_error_test
