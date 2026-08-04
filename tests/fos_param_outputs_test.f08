!> Outputs suite: the cached R(theta) intermediates (#7, #8) and the
!! at-thetas extension.
!!
!! Parity anchor is the OLD 1.x batch evaluator
!! (the tier-1 `compute_radius_and_derivative_standalone_s` at the tables'
!! thetas, resolving the same shift the tier-1 `compute_shape_standalone_s`
!! reports; the 1.x shift for PARAMS7 is frozen alongside it as
!! SHIFT_PARAMS7_1X). The comparison is a 1e-11 relative
!! tolerance, NOT bitwise: the two Newton solves start from different outer
!! brackets (1.x doubles from 2*max(r_north, r_south); the cached path uses the
!! analytic bound built from rho_max), so the bisection paths differ even though
!! both land on the same root.
!!
!! The at-thetas form must agree with the fixed grid BITWISE when fed the same
!! theta nodes — it is the same bundle and the same cos(theta) values, so any
!! difference would mean the two paths do not share a solver.
program fos_param_outputs_test

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_recompute_count_f, &
            compute_shape_standalone_s, compute_radius_and_derivative_standalone_s, &
            compute_rho_at_z_s, &
            FOS_ERROR_INVALID_C, FOS_ERROR_BUFFER_MISMATCH, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_NOT_STAR_CONVEX
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            SHAPE_ERROR_INVALID_GRID
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_PARAMS = 7_ik
    integer(kind = ik), parameter :: N_THETA = 9_ik

    real(kind = rk), parameter :: PARAMS7(N_PARAMS) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Total z-shift the deleted 1.x `compute_fos_shape_s` resolved for PARAMS7
    !! at N_POINTS, frozen from commit 648428c — the numeric anchor the tier
    !! comparison alone cannot supply.
    real(kind = rk), parameter :: SHIFT_PARAMS7_1X = 6.565141402540683E-002_rk

    !> Odd-heavy asymmetric family; raising c drives g(s*) through the
    !! star-convexity margin (same probe family as the resolve suite).
    real(kind = rk), parameter :: ODD_HEAVY(N_PARAMS) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]

    type(cache_t)      :: cache
    integer(kind = ik) :: status, i, code
    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: radii(N_THETA), radii_d(N_THETA), dr_dtheta(N_THETA)
    real(kind = rk)    :: ref_radii(N_THETA), ref_dr(N_THETA)
    real(kind = rk)    :: big_radii(N_THETA + 3_ik)
    real(kind = rk)    :: ref_shift, ref_north, ref_south
    real(kind = rk)    :: params_beak(N_PARAMS), params_star(N_PARAMS)
    real(kind = rk)    :: at_radii(N_THETA), at_dr(N_THETA)
    real(kind = rk)    :: few_thetas(3), few_radii(3), few_dr(3)
    real(kind = rk)    :: bad_thetas(N_THETA)
    integer(kind = ik) :: before(8)
    logical            :: all_bits_eq

    do i = 1_ik, N_THETA
        thetas(i) = real(i - 1_ik, rk) * PI_C / real(N_THETA - 1_ik, rk)
    end do

    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache init valid')

    !---------------------------------------------------------------------------
    ! Parameter-vector gates, ahead of everything else
    !---------------------------------------------------------------------------
    radii = 1.0_rk
    call cache_radius_grid_s(cache, PARAMS7(1:3), radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'short params -> 4')
    call assert_true(all_zero_f(radii), 'short params zero-fills radii')

    radii = 1.0_rk
    params_star = 0.0_rk
    params_star(1) = 1.0e-11_rk
    call cache_radius_grid_s(cache, params_star, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'degenerate c -> 102')
    call assert_true(all_zero_f(radii), '102 zero-fills radii')

    !---------------------------------------------------------------------------
    ! Buffer mismatch: 105 on both fixed-grid forms and on the at-thetas form
    !---------------------------------------------------------------------------
    big_radii = 1.0_rk
    call cache_radius_grid_s(cache, PARAMS7, big_radii, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'oversized radii -> 105')
    call assert_true(all_zero_f(big_radii), '105 zero-fills the whole extent')

    radii = 1.0_rk
    dr_dtheta = 1.0_rk
    call cache_radius_and_derivative_s(cache, PARAMS7, radii(1:N_THETA - 1_ik), &
            dr_dtheta, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'undersized radii -> 105')

    few_radii = 1.0_rk
    few_dr = 1.0_rk
    few_thetas = [0.0_rk, 0.5_rk * PI_C, PI_C]
    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, few_thetas, &
            few_radii(1:2), few_dr, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'at_thetas size mismatch -> 105')

    !---------------------------------------------------------------------------
    ! Parity with the 1.x batch evaluator
    !---------------------------------------------------------------------------
    call probe_tier1_s(PARAMS7, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, SHAPE_VALID, 'params7 accepted by the tier-1 resolve')
    call assert_close(ref_shift, SHIFT_PARAMS7_1X, 1.0e-12_rk, &
            'params7 tier-1 z_shift matches the frozen 1.x value')
    call compute_radius_and_derivative_standalone_s(PARAMS7, thetas, N_POINTS, &
            ref_radii, ref_dr, code)
    call assert_int_eq(code, SHAPE_VALID, 'tier-1 batch evaluation valid')

    call cache_radius_grid_s(cache, PARAMS7, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius grid valid')
    do i = 1_ik, N_THETA
        call assert_close(radii(i), ref_radii(i), 1.0e-11_rk, 'radius grid parity')
    end do

    call cache_radius_and_derivative_s(cache, PARAMS7, radii_d, dr_dtheta, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius+derivative valid')
    do i = 1_ik, N_THETA
        call assert_close(radii_d(i), ref_radii(i), 1.0e-11_rk, 'radius (deriv form) parity')
        call assert_close(dr_dtheta(i), ref_dr(i), 1.0e-11_rk, 'dR/dtheta parity')
    end do

    ! The pole nodes take the analytic branch in both surfaces
    call assert_abs_close(radii(1), ref_north, 0.0_rk, 'north pole radius exact')
    call assert_abs_close(radii(N_THETA), ref_south, 0.0_rk, 'south pole radius exact')
    call assert_abs_close(dr_dtheta(1), 0.0_rk, 0.0_rk, 'north pole dR/dtheta = 0')
    call assert_abs_close(dr_dtheta(N_THETA), 0.0_rk, 0.0_rk, 'south pole dR/dtheta = 0')

    !---------------------------------------------------------------------------
    ! at_thetas on the tables' own nodes is bitwise identical to the fixed grid
    !---------------------------------------------------------------------------
    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, thetas, &
            at_radii, at_dr, status)
    call assert_int_eq(status, SHAPE_VALID, 'at_thetas on the tables nodes valid')
    all_bits_eq = .true.
    do i = 1_ik, N_THETA
        if (.not. bits_eq_f(at_radii(i), radii_d(i))) all_bits_eq = .false.
        if (.not. bits_eq_f(at_dr(i), dr_dtheta(i))) all_bits_eq = .false.
    end do
    call assert_true(all_bits_eq, 'at_thetas is bitwise equal to the fixed grid')

    ! A shorter, unrelated node set is the whole point of the at-thetas form
    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, few_thetas, &
            few_radii, few_dr, status)
    call assert_int_eq(status, SHAPE_VALID, 'at_thetas on a 3-node set valid')
    call assert_abs_close(few_radii(1), ref_north, 0.0_rk, 'at_thetas north pole')
    call assert_abs_close(few_radii(3), ref_south, 0.0_rk, 'at_thetas south pole')

    !---------------------------------------------------------------------------
    ! Out-of-range theta
    !---------------------------------------------------------------------------
    bad_thetas = thetas
    bad_thetas(5) = 1.5_rk * PI_C
    at_radii = 1.0_rk
    at_dr = 1.0_rk
    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, bad_thetas, &
            at_radii, at_dr, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'theta = 3pi/2 -> invalid grid')
    call assert_true(all_zero_f(at_radii) .and. all_zero_f(at_dr), &
            'invalid grid zero-fills both outputs')

    bad_thetas = thetas
    bad_thetas(2) = -1.0e-6_rk
    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, bad_thetas, &
            at_radii, at_dr, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'negative theta -> invalid grid')

    !---------------------------------------------------------------------------
    ! Minimality: #7 and #8 are stamped, at_thetas stamps nothing
    !---------------------------------------------------------------------------
    call cache_radius_grid_s(cache, PARAMS7, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius grid valid before counting')
    call cache_radius_and_derivative_s(cache, PARAMS7, radii_d, dr_dtheta, status)
    call assert_int_eq(status, SHAPE_VALID, 'radius+derivative valid before counting')
    do i = 1_ik, 8_ik
        before(i) = count_f(i)
    end do

    call cache_radius_grid_s(cache, PARAMS7, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'repeat radius grid valid')
    call cache_radius_and_derivative_s(cache, PARAMS7, radii_d, dr_dtheta, status)
    call assert_int_eq(status, SHAPE_VALID, 'repeat radius+derivative valid')
    do i = 1_ik, 8_ik
        call assert_int_eq(count_f(i), before(i), 'repeat call recomputes nothing')
    end do

    call cache_radius_and_derivative_at_thetas_s(cache, PARAMS7, few_thetas, &
            few_radii, few_dr, status)
    call assert_int_eq(status, SHAPE_VALID, 'at_thetas valid while counting')
    call assert_int_eq(count_f(7_ik), before(7), 'at_thetas does not stamp #7')
    call assert_int_eq(count_f(8_ik), before(8), 'at_thetas does not stamp #8')

    !---------------------------------------------------------------------------
    ! Shape gates propagate to every R(theta) form
    !---------------------------------------------------------------------------
    params_beak = 0.0_rk
    params_beak(1) = 2.0_rk
    params_beak(3) = 0.7495_rk
    ! a4 = 0.7495 is the vector the 1.x beak probe settled on, frozen here.
    call probe_tier1_s(params_beak, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
            'probe vector is beak-rejected by the tier-1 resolve')
    call check_rejection_s(params_beak, FOS_ERROR_BEAK_SINGULARITY, 'beak vector -> 103')

    ! c = 1.3 is the first elongation the 1.x probe rejected with 101, frozen.
    params_star = ODD_HEAVY
    params_star(1) = 1.3_rk
    call probe_tier1_s(params_star, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, FOS_ERROR_NOT_STAR_CONVEX, &
            'probe vector is 101-rejected by the tier-1 resolve')
    call check_rejection_s(params_star, FOS_ERROR_NOT_STAR_CONVEX, &
            'non-star-convex vector -> 101')

    call cache_free_s(cache)

    call extreme_oblate_s()

    call test_summary()

contains

    !> Recompute counter as an `ik` integer, for the assertion helpers.
    function count_f(intermediate) result(n)
        integer(kind = ik), intent(in) :: intermediate
        integer(kind = ik) :: n
        n = int(cache_recompute_count_f(cache, intermediate), ik)
    end function count_f

    !> Bit-level equality without a real `==` (banned by -Werror).
    pure function bits_eq_f(a, b) result(same)
        real(kind = rk), intent(in) :: a, b
        logical :: same
        same = transfer(a, 0_ikl) == transfer(b, 0_ikl)
    end function bits_eq_f

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

    !> The tier-1 resolve at the cache's own resolution, for its verdict and
    !! its z_shift — the cross-tier oracle that replaced the 1.x resolve.
    subroutine probe_tier1_s(p, tier1_code, tier1_shift, tier1_north, tier1_south)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(out) :: tier1_code
        real(kind = rk), intent(out) :: tier1_shift, tier1_north, tier1_south

        call compute_shape_standalone_s(p, N_POINTS, tier1_shift, tier1_north, &
                tier1_south, tier1_code)

    end subroutine probe_tier1_s

    !> A rejected shape yields the same code and a zero-filled output on all
    !! three R(theta) forms.
    subroutine check_rejection_s(p, want_code, label)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(in) :: want_code
        character(len = *), intent(in) :: label

        real(kind = rk) :: out_r(N_THETA), out_d(N_THETA), out_t(N_THETA)
        integer(kind = ik) :: code_out

        out_r = 1.0_rk
        call cache_radius_grid_s(cache, p, out_r, code_out)
        call assert_int_eq(code_out, want_code, label // ' (grid)')
        call assert_true(all_zero_f(out_r), label // ': grid zeroed')

        out_r = 1.0_rk
        out_d = 1.0_rk
        call cache_radius_and_derivative_s(cache, p, out_r, out_d, code_out)
        call assert_int_eq(code_out, want_code, label // ' (grid+deriv)')
        call assert_true(all_zero_f(out_r) .and. all_zero_f(out_d), &
                label // ': grid+deriv zeroed')

        out_r = 1.0_rk
        out_t = 1.0_rk
        call cache_radius_and_derivative_at_thetas_s(cache, p, thetas, out_r, &
                out_t, code_out)
        call assert_int_eq(code_out, want_code, label // ' (at_thetas)')
        call assert_true(all_zero_f(out_r) .and. all_zero_f(out_t), &
                label // ': at_thetas zeroed')

    end subroutine check_rejection_s

    !> Extreme-oblate regression: c = 2e-10 is a pancake whose equatorial radius
    !! (sqrt(1/c) ~ 7.07e4) exceeds the polar extent by 14 orders of magnitude.
    !! 1.x doubled its bracket only 8 times and silently returned the initial
    !! guess; the analytic bracket must reach the true root.
    !!
    !! Probed, not assumed: if the star-convexity margin rejects the pancake the
    !! test walks c up until a vector is accepted and checks the same invariants
    !! on the most extreme shape that passes.
    subroutine extreme_oblate_s()

        type(cache_t)      :: flat
        real(kind = rk)    :: p1(1), t3(3), r3(3), d3(3)
        real(kind = rk)    :: z_shift, r_north, r_south, rho_surface, residual
        real(kind = rk)    :: c_used
        integer(kind = ik) :: st, j
        logical            :: accepted

        t3 = [0.0_rk, 0.5_rk * PI_C, PI_C]

        accepted = .false.
        c_used = 0.0_rk
        do j = 0_ik, 12_ik
            p1(1) = 2.0e-10_rk * 10.0_rk ** real(j, rk)
            call cache_init_s(flat, 1_ik, N_POINTS, t3, st)
            call assert_int_eq(st, SHAPE_VALID, 'flat cache init valid')
            call cache_radius_and_derivative_at_thetas_s(flat, p1, t3, r3, d3, st)
            if (st == SHAPE_VALID) then
                accepted = .true.
                c_used = p1(1)
                exit
            end if
            call cache_free_s(flat)
        end do

        call assert_true(accepted, 'an extreme-oblate vector is accepted')
        if (.not. accepted) return

        call assert_abs_close(c_used, 2.0e-10_rk, 1.0e-18_rk, &
                'c = 2e-10 itself passes validation')

        ! Equatorial radius: sqrt(1/c) up to the u-offset of cos(pi/2) ~ 6e-17
        call assert_close(r3(2), sqrt(1.0_rk / c_used), 1.0e-2_rk, &
                'equatorial radius reaches sqrt(1/c)')
        call assert_true(r3(2) > 1.0e4_rk, 'equatorial radius is not the initial guess')

        ! Residual check on the returned root, scaled like the solver's own
        ! tolerance (NR_TOLERANCE * max(1, r)).
        call cache_shape_s(flat, p1, z_shift, r_north, r_south, st)
        call assert_int_eq(st, SHAPE_VALID, 'flat shape valid')
        call compute_rho_at_z_s(p1, r3(2) * cos(t3(2)), z_shift, rho_surface)
        residual = abs(r3(2) * sin(t3(2)) - rho_surface)
        call assert_true(residual < 1.0e-9_rk * max(1.0_rk, r3(2)), &
                'equatorial root satisfies the surface equation')

        call cache_free_s(flat)

    end subroutine extreme_oblate_s

end program fos_param_outputs_test
