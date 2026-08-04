!> Anchor tests: sphere and spheroids, where every quantity has a closed form.
program fos_param_core_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_radius_grid_standalone_s, FOS_VALID, &
            compute_shape_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY
    use fos_test_reference_mod, only: init_quadrature_s, compute_reference_surface_f, &
            evaluate_shape_quality_s, spheroid_surface_area_f
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID_SMALL = 101_ik
    real(kind = rk), parameter :: TOL_VOLUME_A = 1.0e-11_rk
    real(kind = rk), parameter :: TOL_ROUND_TRIP = 1.0e-9_rk

    !> Uniform theta nodes of the R(theta) grid, built once. The 1.x entry
    !! point owned this grid; the tier-1 form takes the nodes explicitly.
    real(kind = rk) :: grid_thetas(N_GRID_SMALL)
    integer(kind = ik) :: it

    ! Endpoints pinned: under -ffast-math the last product can land one ulp
    ! above PI_C, which the library rejects as an out-of-domain theta.
    do it = 1_ik, N_GRID_SMALL
        grid_thetas(it) = real(it - 1_ik, rk) * PI_C / real(N_GRID_SMALL - 1_ik, rk)
    end do
    grid_thetas(1) = 0.0_rk
    grid_thetas(N_GRID_SMALL) = PI_C

    call init_quadrature_s()
    call run_anchor_tests_s()
    call run_marginal_star_convex_test_s()
    call run_star_convexity_optimum_test_s()
    call test_summary()

contains

    subroutine run_anchor_tests_s()
        real(kind = rk), parameter :: SPHEROID_C(4) = [0.8_rk, 1.5_rk, 2.0_rk, 3.0_rk]

        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift, r_north, r_south
        real(kind = rk) :: s_ref, s_exact, dv_rel, ds_rel, rt_max
        character(len = 64) :: label
        integer(kind = ik) :: i, status

        write(*, '(A)') '=== Anchor tests (sphere, spheroids) ==='

        ! Sphere: c = 1, all a_k = 0
        params = 0.0_rk
        params(1) = 1.0_rk
        call resolve_and_grid_s(params, z_shift, r_north, r_south, radii, status)
        call assert_int_eq(status, FOS_VALID, 'sphere: shape valid')
        call assert_abs_close(z_shift, 0.0_rk, 1.0e-14_rk, 'sphere: z_shift = 0')
        call assert_abs_close(maxval(abs(radii - 1.0_rk)), 0.0_rk, 5.0e-12_rk, 'sphere: R(theta) = 1')

        ! For the sphere the S_ref integrand is exactly 1, so S_ref = 4 pi.
        s_ref = compute_reference_surface_f(params)
        call assert_abs_close(s_ref, 4.0_rk * PI_C, 1.0e-12_rk, 'sphere: S_ref = 4 pi')

        call evaluate_shape_quality_s(params, N_GRID_SMALL, z_shift, s_ref, &
                dv_rel, ds_rel, rt_max, status)
        call assert_int_eq(status, FOS_VALID, 'sphere: quality probe valid')
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

            call resolve_and_grid_s(params, z_shift, r_north, r_south, radii, status)
            call assert_int_eq(status, FOS_VALID, trim(label) // ': shape valid')
            if (status /= FOS_VALID) cycle

            ! The analytic pole radii must be the grid's own endpoints.
            call assert_abs_close(r_north, radii(1), 0.0_rk, trim(label) // ': r_north = R(0)')
            call assert_abs_close(r_south, radii(N_GRID_SMALL), 0.0_rk, &
                    trim(label) // ': r_south = R(pi)')

            call evaluate_shape_quality_s(params, N_GRID_SMALL, z_shift, s_exact, &
                    dv_rel, ds_rel, rt_max, status)
            call assert_int_eq(status, FOS_VALID, trim(label) // ': quality probe valid')
            call assert_abs_close(dv_rel, 0.0_rk, TOL_VOLUME_A, trim(label) // ': volume')
            call assert_abs_close(ds_rel, 0.0_rk, 1.0e-10_rk, trim(label) // ': surface vs closed form')
            call assert_abs_close(rt_max, 0.0_rk, TOL_ROUND_TRIP, trim(label) // ': round trip')
        end do
    end subroutine run_anchor_tests_s

    !> Tier-1 resolve + R(theta) grid at the uniform nodes, the 2.0 replacement
    !! for the 1.x single-call entry point. Both calls run at the same
    !! n_points, so they resolve the same z_shift.
    subroutine resolve_and_grid_s(params, z_shift, r_north, r_south, radii, status)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift, r_north, r_south
        real(kind = rk), intent(out) :: radii(N_GRID_SMALL)
        integer(kind = ik), intent(out) :: status

        radii = 0.0_rk
        call compute_shape_standalone_s(params, N_GRID_SMALL, z_shift, r_north, &
                r_south, status)
        if (status /= FOS_VALID) return
        call compute_radius_grid_standalone_s(params, grid_thetas, N_GRID_SMALL, &
                radii, status)
    end subroutine resolve_and_grid_s

    !> The disputed shape (c=2, a3=0.4, a4=0.66) is genuinely star-convex: its
    !! best origin gives max-T = -0.103 R0, past the -0.1 margin. The old coarse
    !! shift search missed the ~0.0085 R0-wide acceptance window and rejected it.
    !! Accepting it must yield an R(theta) that round-trips to quadrature precision.
    subroutine run_marginal_star_convex_test_s()
        real(kind = rk) :: params(7), radii(N_GRID_SMALL), z_shift, r_north, r_south
        real(kind = rk) :: s_ref, dv_rel, ds_rel, rt_max
        integer(kind = ik) :: code

        write(*, '(A)') '=== Marginal star-convex shape (c=2, a3=0.4, a4=0.66) ==='
        params = [2.0_rk, 0.4_rk, 0.66_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]

        call resolve_and_grid_s(params, z_shift, r_north, r_south, radii, code)
        call assert_int_eq(code, FOS_VALID, 'marginal: accepted (was falsely rejected)')
        call assert_true(abs(z_shift) > 1.0e-6_rk, 'marginal: nonzero star-convex shift')

        if (code == FOS_VALID) then
            s_ref = compute_reference_surface_f(params)
            call evaluate_shape_quality_s(params, N_GRID_SMALL, z_shift, s_ref, &
                    dv_rel, ds_rel, rt_max, code)
            call assert_int_eq(code, FOS_VALID, 'marginal: quality probe valid')
            call assert_true(rt_max < TOL_ROUND_TRIP, 'marginal: R(theta) round-trips')
        end if
    end subroutine run_marginal_star_convex_test_s

    subroutine run_star_convexity_optimum_test_s()
        real(kind = rk) :: params(7), z_shift_total, max_t_opt
        integer(kind = ik) :: code
        integer(kind = ik), parameter :: N = 7201_ik

        write(*, '(A)') '=== Star-convexity optimum diagnostic ==='

        ! Marginal F5: g(s*) = -0.1030 (fixed geometric value, ~just past -0.1),
        ! total shift equals the F5 golden z_shift (-0.16904663...).
        params = [2.0_rk, 0.4_rk, 0.66_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        call compute_star_convexity_optimum_standalone_s(params, N, z_shift_total, &
                max_t_opt, code)
        call assert_int_eq(code, FOS_VALID, 'optimum F5: FOS_VALID')
        call assert_abs_close(max_t_opt, -0.1030_rk, 2.0e-3_rk, 'optimum F5: g(s*) ~ -0.103')
        call assert_abs_close(z_shift_total, -0.16904663095809774_rk, 1.0e-4_rk, &
                'optimum F5: total shift matches golden')

        ! a4=0.67: g(s*) = -0.0873 -- single-valued in principle (g < 0) yet
        ! closer to the grazing-ray limit than F5. Fixed bounds, margin-independent.
        params = [2.0_rk, 0.4_rk, 0.67_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        call compute_star_convexity_optimum_standalone_s(params, N, z_shift_total, &
                max_t_opt, code)
        call assert_int_eq(code, FOS_VALID, 'optimum a4=0.67: FOS_VALID')
        call assert_true(max_t_opt < 0.0_rk, 'optimum a4=0.67: representable (single-valued)')
        call assert_true(max_t_opt > -0.1_rk, 'optimum a4=0.67: closer to limit than F5')

        ! Sphere: strongly star-convex, optimum at s=0, g(s*) = -1 exactly.
        params = 0.0_rk
        params(1) = 1.0_rk
        call compute_star_convexity_optimum_standalone_s(params, N, z_shift_total, &
                max_t_opt, code)
        call assert_int_eq(code, FOS_VALID, 'optimum sphere: FOS_VALID')
        call assert_abs_close(max_t_opt, -1.0_rk, 1.0e-3_rk, 'optimum sphere: g(s*) = -1')
        call assert_abs_close(z_shift_total, 0.0_rk, 1.0e-3_rk, 'optimum sphere: zero shift')

        ! Degenerate: c=0 -> invalid c; beak family -> beak singularity.
        params = 0.0_rk
        call compute_star_convexity_optimum_standalone_s(params, N, z_shift_total, &
                max_t_opt, code)
        call assert_int_eq(code, FOS_ERROR_INVALID_C, 'optimum c=0: FOS_ERROR_INVALID_C')

        params = [2.0_rk, 0.0_rk, 0.7497_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        call compute_star_convexity_optimum_standalone_s(params, N, z_shift_total, &
                max_t_opt, code)
        call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, 'optimum beak: FOS_ERROR_BEAK_SINGULARITY')
    end subroutine run_star_convexity_optimum_test_s

end program fos_param_core_test
