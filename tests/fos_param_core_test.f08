!> Anchor tests: sphere and spheroids, where every quantity has a closed form.
program fos_param_core_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_fos_radius_grid_s, FOS_VALID
    use fos_test_reference_mod, only: init_quadrature_s, compute_reference_surface_f, &
            evaluate_shape_quality_s, spheroid_surface_area_f
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID_SMALL = 101_ik
    real(kind = rk), parameter :: TOL_VOLUME_A = 1.0e-11_rk
    real(kind = rk), parameter :: TOL_ROUND_TRIP = 1.0e-9_rk

    call init_quadrature_s()
    call run_anchor_tests_s()
    call run_marginal_star_convex_test_s()
    call test_summary()

contains

    subroutine run_anchor_tests_s()
        real(kind = rk), parameter :: SPHEROID_C(4) = [0.8_rk, 1.5_rk, 2.0_rk, 3.0_rk]

        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
        real(kind = rk) :: s_ref, s_exact, dv_rel, ds_rel, rt_max
        logical :: is_valid
        character(len = 256) :: message
        character(len = 64) :: label
        integer(kind = ik) :: i

        write(*, '(A)') '=== Anchor tests (sphere, spheroids) ==='

        ! Sphere: c = 1, all a_k = 0
        params = 0.0_rk
        params(1) = 1.0_rk
        call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, message)
        call assert_true(is_valid, 'sphere: shape valid')
        call assert_abs_close(z_shift, 0.0_rk, 1.0e-14_rk, 'sphere: z_shift = 0')
        call assert_abs_close(maxval(abs(radii - 1.0_rk)), 0.0_rk, 5.0e-12_rk, 'sphere: R(theta) = 1')

        ! For the sphere the S_ref integrand is exactly 1, so S_ref = 4 pi.
        s_ref = compute_reference_surface_f(params)
        call assert_abs_close(s_ref, 4.0_rk * PI_C, 1.0e-12_rk, 'sphere: S_ref = 4 pi')

        call evaluate_shape_quality_s(params, z_shift, s_ref, dv_rel, ds_rel, rt_max)
        call assert_abs_close(dv_rel, 0.0_rk, TOL_VOLUME_A, 'sphere: volume')
        call assert_abs_close(ds_rel, 0.0_rk, 1.0e-11_rk, 'sphere: surface')
        call assert_abs_close(rt_max, 0.0_rk, 1.0e-11_rk, 'sphere: round trip')

        ! Spheroids: c /= 1, all a_k = 0 -> semi-axes (c, 1/sqrt(c)), V = 4 pi / 3
        do i = 1_ik, size(SPHEROID_C, kind = ik)
            params = 0.0_rk
            params(1) = SPHEROID_C(i)
            write(label, '(A,F4.2)') 'spheroid c=', SPHEROID_C(i)

            ! Validate the S_ref machinery itself against the closed form
            s_exact = spheroid_surface_area_f(SPHEROID_C(i))
            s_ref = compute_reference_surface_f(params)
            call assert_abs_close(s_ref, s_exact, 1.0e-10_rk * s_exact, &
                    trim(label) // ': S_ref vs closed form')

            call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, message)
            call assert_true(is_valid, trim(label) // ': shape valid')
            if (.not. is_valid) cycle

            call evaluate_shape_quality_s(params, z_shift, s_exact, dv_rel, ds_rel, rt_max)
            call assert_abs_close(dv_rel, 0.0_rk, TOL_VOLUME_A, trim(label) // ': volume')
            call assert_abs_close(ds_rel, 0.0_rk, 1.0e-10_rk, trim(label) // ': surface vs closed form')
            call assert_abs_close(rt_max, 0.0_rk, TOL_ROUND_TRIP, trim(label) // ': round trip')
        end do
    end subroutine run_anchor_tests_s

    !> The disputed shape (c=2, a3=0.4, a4=0.66) is genuinely star-convex: its
    !! best origin gives max-T = -0.103 R0, past the -0.1 margin. The old coarse
    !! shift search missed the ~0.0085 R0-wide acceptance window and rejected it.
    !! Accepting it must yield an R(theta) that round-trips to quadrature precision.
    subroutine run_marginal_star_convex_test_s()
        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift
        real(kind = rk) :: s_ref, dv_rel, ds_rel, rt_max
        logical :: is_valid
        character(len = 256) :: message
        integer(kind = ik) :: code

        write(*, '(A)') '=== Marginal star-convex shape (c=2, a3=0.4, a4=0.66) ==='
        params = [2.0_rk, 0.4_rk, 0.66_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]

        call compute_fos_radius_grid_s(params, N_GRID_SMALL, radii, z_shift, is_valid, &
                message, error_code = code)
        call assert_true(is_valid, 'marginal: accepted (was falsely rejected)')
        call assert_int_eq(code, FOS_VALID, 'marginal: FOS_VALID')
        call assert_true(abs(z_shift) > 1.0e-6_rk, 'marginal: nonzero star-convex shift')

        if (is_valid) then
            s_ref = compute_reference_surface_f(params)
            call evaluate_shape_quality_s(params, z_shift, s_ref, dv_rel, ds_rel, rt_max)
            call assert_true(rt_max < TOL_ROUND_TRIP, 'marginal: R(theta) round-trips')
        end if
    end subroutine run_marginal_star_convex_test_s

end program fos_param_core_test
