!> Worker-kernel suite: status-reporting a2/z_shift, tabled f-grid and beak
!! scan, rho scaling, and the bounded-bracket Newton radius core.
!!
!! Parity anchors are the 1.x publics (`compute_fos_a2_f`,
!! `compute_fos_z_shift_f`, `compute_fos_f_and_derivatives_s`): the kernels
!! must reproduce them node-for-node, bitwise where the arithmetic is
!! identical.
program fos_param_workers_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_a2_s, compute_z_shift_s, &
            compute_fos_a2_f, compute_fos_z_shift_f, &
            compute_fos_f_and_derivatives_s, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY
    use fos_parameterization_workers_mod, only: tables_t, tables_init_s, &
            tables_free_s, fos_bundle_t, compute_f_grid_s, beak_scan_f_min_s, &
            scale_rho_grid_s, newton_radius_s
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    real(kind = rk), parameter :: PARAMS7(7) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    type(tables_t)     :: tables
    type(fos_bundle_t) :: bundle
    integer(kind = ik) :: status, i, err_code
    real(kind = rk)    :: thetas(4), params51(51), params_beak(3)
    real(kind = rk)    :: a2_new, zs, f_ref, fp_ref, f_min, a4
    real(kind = rk)    :: r, dr, rho_max, c_oblate
    real(kind = rk)    :: f_grid(N_POINTS), fp_grid(N_POINTS)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
    logical            :: beak_ok, rho_positive, converged
    ! NOTE (-Werror hygiene): declare ONLY what this program uses; the
    ! contained probe declares its own scratch outputs.

    do i = 1_ik, 4_ik
        thetas(i) = real(i, rk) * PI_C / 5.0_rk
    end do

    call tables_init_s(tables, N_POINTS, thetas, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'tables built')

    !---------------------------------------------------------------------------
    ! a2: parity with the 1.x pure function, empty ok, oversize rejected
    !---------------------------------------------------------------------------
    call compute_a2_s(PARAMS7, a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'a2 valid')
    call assert_abs_close(a2_new, compute_fos_a2_f(PARAMS7), 0.0_rk, 'a2 bitwise-parity')

    call compute_a2_s(PARAMS7(1:0), a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'a2 empty ok')
    call assert_abs_close(a2_new, 0.0_rk, 0.0_rk, 'a2 empty = 0')

    do i = 1_ik, 51_ik
        params51(i) = 0.01_rk * real(i, rk)
    end do
    call compute_a2_s(params51, a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_TOO_MANY_PARAMS), 'a2 51 params rejected')

    !---------------------------------------------------------------------------
    ! z_shift: parity, empty and degenerate c -> 102
    !---------------------------------------------------------------------------
    call compute_z_shift_s(PARAMS7, zs, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'z_shift valid')
    call assert_abs_close(zs, compute_fos_z_shift_f(PARAMS7), 0.0_rk, 'z_shift parity')

    call compute_z_shift_s(PARAMS7(1:0), zs, status)
    call assert_int_eq(int(status), int(FOS_ERROR_INVALID_C), 'z_shift empty -> 102')
    call compute_z_shift_s([1.0e-11_rk, 0.1_rk], zs, status)
    call assert_int_eq(int(status), int(FOS_ERROR_INVALID_C), 'z_shift degenerate c -> 102')

    !---------------------------------------------------------------------------
    ! f-grid: parity against the live evaluator at every interior node
    !---------------------------------------------------------------------------
    call compute_f_grid_s(tables, PARAMS7, f_grid, fp_grid)
    do i = 2_ik, tables%n_points - 1_ik
        call compute_fos_f_and_derivatives_s(PARAMS7, tables%u(i), f_ref, fp_ref)
        call assert_abs_close(f_grid(i), f_ref, 1.0e-14_rk, 'f_grid node parity')
        call assert_abs_close(fp_grid(i), fp_ref, 1.0e-13_rk, 'fp_grid node parity')
    end do

    !---------------------------------------------------------------------------
    ! Beak scan: sphere passes; a vector the 1.x surface rejects as a beak fails
    !---------------------------------------------------------------------------
    ! Never hardcode an unverified vector: confirm the rejection against the OLD
    ! surface, and if [2.0, 0.0, 0.72] is not rejected, probe a4 in 0.60..0.80
    ! and adopt the first a4 that is.
    params_beak(1) = 2.0_rk
    params_beak(2) = 0.0_rk
    params_beak(3) = 0.72_rk
    call probe_old_surface_s(params_beak, err_code)
    if (err_code /= FOS_ERROR_BEAK_SINGULARITY) then
        do i = 0_ik, 400_ik
            a4 = 0.60_rk + 5.0e-4_rk * real(i, rk)
            params_beak(3) = a4
            call probe_old_surface_s(params_beak, err_code)
            if (err_code == FOS_ERROR_BEAK_SINGULARITY) exit
        end do
    end if
    call assert_int_eq(int(err_code), int(FOS_ERROR_BEAK_SINGULARITY), &
            'probe vector is beak-rejected by the 1.x surface')

    call beak_scan_f_min_s(tables, [1.0_rk], f_min, beak_ok)
    call assert_true(beak_ok, 'sphere passes beak scan')
    call assert_true(f_min > 0.0_rk, 'sphere f_min positive')

    call beak_scan_f_min_s(tables, params_beak, f_min, beak_ok)
    call assert_true(.not. beak_ok, 'beak vector fails scan')

    !---------------------------------------------------------------------------
    ! rho scaling: sphere geometry, tip convention, shift, rho-positivity verdict
    !---------------------------------------------------------------------------
    call compute_f_grid_s(tables, [1.0_rk], f_grid, fp_grid)
    call scale_rho_grid_s(tables, 1.0_rk, 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_true(rho_positive, 'sphere rho grid positive')
    call assert_abs_close(rho_max, 1.0_rk, 1.0e-12_rk, 'sphere rho_max = 1')
    call assert_abs_close(rho(1), 0.0_rk, 0.0_rk, 'south tip rho = 0')
    call assert_abs_close(drho_dz(N_POINTS), 0.0_rk, 0.0_rk, 'north tip drho/dz = 0')
    i = (N_POINTS + 1_ik) / 2_ik
    call assert_abs_close(z(i), 0.0_rk, 1.0e-15_rk, 'equator z = 0')
    call assert_abs_close(rho(i), 1.0_rk, 1.0e-15_rk, 'equator rho = 1')

    ! z(i) = c*u(i) + z_shift_intrinsic, rho scales as 1/sqrt(c)
    call scale_rho_grid_s(tables, 2.0_rk, 0.5_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_abs_close(z(1), -1.5_rk, 1.0e-14_rk, 'z(1) = -c + z_shift')
    call assert_abs_close(z(N_POINTS), 2.5_rk, 1.0e-14_rk, 'z(n) = c + z_shift')
    call assert_abs_close(rho_max, 1.0_rk / sqrt(2.0_rk), 1.0e-12_rk, 'rho_max = 1/sqrt(c)')

    ! f(0) = 0 exactly for a4 = 0.75: the interior node at u = 0 pinches shut
    call compute_f_grid_s(tables, [1.0_rk, 0.0_rk, 0.75_rk], f_grid, fp_grid)
    call scale_rho_grid_s(tables, 1.0_rk, 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_true(.not. rho_positive, 'pinched shape flagged rho <= 0')

    !---------------------------------------------------------------------------
    ! Newton with the analytic bracket: extreme-oblate correctness
    !---------------------------------------------------------------------------
    ! c = 2e-10, single param: rho(z = 0) = sqrt(f(0)/c) = 1/sqrt(c) ~ 7.07e4.
    ! The 1.x doubling bracket silently returns ~1e-7 here.
    c_oblate = 2.0e-10_rk
    bundle%n_params = 1_ik
    bundle%params(1) = c_oblate
    bundle%z_shift = 0.0_rk
    bundle%r_hi_bound = 2.0_rk * sqrt((1.0_rk / sqrt(c_oblate))**2 + c_oblate**2)
    call newton_radius_s(bundle, 0.0_rk, r, dr, converged)
    call assert_true(converged, 'extreme oblate converges')
    call assert_abs_close(r, 1.0_rk / sqrt(c_oblate), 1.0e-6_rk * 7.1e4_rk, &
            'equatorial radius ~ 1/sqrt(c)')

    ! Sphere sanity: r = 1 everywhere, converged
    bundle%n_params = 1_ik
    bundle%params(1) = 1.0_rk
    bundle%z_shift = 0.0_rk
    bundle%r_hi_bound = 2.0_rk * sqrt(2.0_rk)
    call newton_radius_s(bundle, cos(1.1_rk), r, dr, converged)
    call assert_true(converged, 'sphere converges')
    call assert_abs_close(r, 1.0_rk, 1.0e-10_rk, 'sphere r = 1')
    call assert_abs_close(dr, 0.0_rk, 1.0e-8_rk, 'sphere dR/dtheta = 0')

    ! Poles are analytic and always converged
    call newton_radius_s(bundle, 1.0_rk, r, dr, converged)
    call assert_true(converged, 'north pole converged')
    call assert_abs_close(r, 1.0_rk, 0.0_rk, 'north pole r = c + z_shift')
    call newton_radius_s(bundle, -1.0_rk, r, dr, converged)
    call assert_true(converged, 'south pole converged')
    call assert_abs_close(r, 1.0_rk, 0.0_rk, 'south pole r = |-c + z_shift|')

    call tables_free_s(tables)
    call test_summary()

contains

    !> Runs the 1.x radius-grid entry point purely for its rejection code.
    subroutine probe_old_surface_s(p, code)
        use fos_parameterization_mod, only: compute_fos_radius_grid_s
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(out) :: code

        real(kind = rk) :: radii(64), z_shift_out
        logical :: is_valid
        character(len = 256) :: message

        call compute_fos_radius_grid_s(p, 64_ik, radii, z_shift_out, is_valid, &
                message, n_rho_grid = N_POINTS, error_code = code)

    end subroutine probe_old_surface_s

end program fos_param_workers_test
