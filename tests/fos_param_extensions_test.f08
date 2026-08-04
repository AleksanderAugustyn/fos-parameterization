!> Extensions suite: the two cached tier-1 diagnostics that deliberately do NOT
!! share the R(theta) gate set — `cache_neck_s` and
!! `cache_star_convexity_optimum_s`.
!!
!! Parity anchors are the tier-1 forms (`compute_neck_standalone_s`,
!! `compute_star_convexity_optimum_standalone_s`), which run the same kernels
!! outside the cache, to 1e-11 relative, plus a g(s*) literal frozen from the
!! deleted 1.x diagnostic (commit 648428c): the cached path
!! derives u from the trig tables while 1.x computes z*(1/c), so the comparison
!! is close, not bitwise.
!!
!! The suite's real subject is the GATING ASYMMETRY. The neck lives on the
!! cylindrical profile, so a beak-singular shape still has one; the optimum
!! exists precisely to report g(s*) for shapes the star-convexity margin
!! rejects, so it must not gate on that margin.
!!
!! Gate vectors are never hardcoded blind: each one is confirmed against the
!! 1.x surface in-test, with a scan if the first guess is accepted.
program fos_param_extensions_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_neck_s, cache_star_convexity_optimum_s, cache_shape_s, &
            cache_rho_z_grid_s, cache_radius_grid_s, cache_recompute_count_f, &
            compute_neck_standalone_s, compute_star_convexity_optimum_standalone_s, &
            compute_shape_standalone_s, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_NOT_STAR_CONVEX, &
            STAR_CONVEXITY_MARGIN
    use shape_core_mod, only: SHAPE_VALID
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    !> The 1.x neck scanner defaulted to a 1001-point grid; the cache and the
    !! tier-1 reference are both built at that resolution so every parity
    !! comparison here compares like with like.
    integer(kind = ik), parameter :: N_POINTS = 1001_ik
    integer(kind = ik), parameter :: N_PARAMS = 7_ik
    integer(kind = ik), parameter :: N_THETA = 4_ik

    !> Parity tolerance, mixed absolute/relative (see Task 6's rationale).
    real(kind = rk), parameter :: TOL_PARITY = 1.0e-11_rk

    !> Elongations of the symmetric (c, a4) neck family, from the neck suite.
    real(kind = rk), parameter :: NECK_C(5) = &
            [1.0_rk, 1.5_rk, 2.0_rk, 2.5_rk, 3.0_rk]

    real(kind = rk), parameter :: PARAMS7(N_PARAMS) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> The ungated star-convexity optimum the deleted 1.x diagnostic reported
    !! for PARAMS7 at N_POINTS, frozen from commit 648428c. Both the shift and
    !! g(s*) are geometric properties of the shape, so they anchor the 2.0
    !! optimum to a value no 2.0 code path produced.
    real(kind = rk), parameter :: OPT_SHIFT_PARAMS7_1X = 3.7350130220951994E-001_rk
    real(kind = rk), parameter :: OPT_GOPT_PARAMS7_1X = -8.183218439306816E-001_rk

    !> Odd-heavy asymmetric family: raising c drives g(s*) up through
    !! -STAR_CONVEXITY_MARGIN, which is how the 101 vector is probed.
    real(kind = rk), parameter :: ODD_HEAVY(N_PARAMS) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]

    type(cache_t)      :: cache
    integer(kind = ik) :: status, code, i, j
    real(kind = rk)    :: thetas(N_THETA), radii(N_THETA)
    real(kind = rk)    :: a4
    real(kind = rk)    :: params_neck(N_PARAMS), params_beak(N_PARAMS)
    real(kind = rk)    :: params_star(N_PARAMS)
    real(kind = rk)    :: z_neck, rho_neck, ref_z, ref_rho
    real(kind = rk)    :: z_shift_total, g_opt, ref_shift_total, ref_g_opt
    real(kind = rk)    :: z_shift, r_north, r_south, grid_shift
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
    integer(kind = ik) :: before(8)
    logical            :: found, ref_found
    character(len = 64) :: label

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / 5.0_rk
    end do

    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache init valid')

    !---------------------------------------------------------------------------
    ! Neck parity with the 1.x scanner on the symmetric (c, a4) family
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Neck parity with the tier-1 scanner ==='

    do i = 1_ik, 5_ik
        do j = 0_ik, 12_ik
            a4 = 0.10_rk + 0.05_rk * real(j, rk)
            params_neck = 0.0_rk
            params_neck(1) = NECK_C(i)
            params_neck(3) = a4
            write(label, '(A,F4.2,A,F4.2)') 'neck c=', NECK_C(i), ' a4=', a4

            call compute_neck_standalone_s(params_neck, N_POINTS, ref_z, ref_rho, &
                    ref_found, code)
            call assert_int_eq(code, SHAPE_VALID, trim(label) // ': tier-1 status valid')
            call cache_neck_s(cache, params_neck, z_neck, rho_neck, found, status)

            call assert_int_eq(status, SHAPE_VALID, trim(label) // ': status valid')
            call assert_true(found .eqv. ref_found, trim(label) // ': found matches tier-1')
            if (.not. found) cycle
            call assert_close(z_neck, ref_z, TOL_PARITY, trim(label) // ': z_neck parity')
            call assert_close(rho_neck, ref_rho, TOL_PARITY, &
                    trim(label) // ': rho_neck parity')
        end do
    end do

    !---------------------------------------------------------------------------
    ! Star-convexity optimum parity with the 1.x diagnostic
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Star-convexity optimum parity (tier-1 + frozen 1.x) ==='

    call compute_star_convexity_optimum_standalone_s(PARAMS7, N_POINTS, &
            ref_shift_total, ref_g_opt, code)
    call assert_int_eq(code, SHAPE_VALID, 'params7 tier-1 optimum reports valid')
    call assert_close(ref_shift_total, OPT_SHIFT_PARAMS7_1X, TOL_PARITY, &
            'params7 optimum shift matches the frozen 1.x value')
    call assert_close(ref_g_opt, OPT_GOPT_PARAMS7_1X, TOL_PARITY, &
            'params7 optimum g(s*) matches the frozen 1.x value')

    call cache_star_convexity_optimum_s(cache, PARAMS7, z_shift_total, g_opt, status)
    call assert_int_eq(status, SHAPE_VALID, 'params7 optimum valid')
    call assert_close(z_shift_total, ref_shift_total, TOL_PARITY, &
            'params7 optimum z_shift_total parity')
    call assert_close(g_opt, ref_g_opt, TOL_PARITY, 'params7 optimum g(s*) parity')

    ! params7 is well-conditioned at its COM, so the R(theta) origin stays there
    ! while s* does not: the two shifts must NOT agree. This is what separates
    ! "intrinsic + s*" from the resolve's stored, branch-selected z_shift_total.
    call cache_shape_s(cache, PARAMS7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'params7 shape valid')
    call assert_true(abs(z_shift_total - z_shift) > 1.0e-6_rk, &
            'params7 optimum shift differs from the branch-selected R(theta) origin')

    !---------------------------------------------------------------------------
    ! Gating asymmetry, case 1: a beak-singular shape still has a neck
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Gating asymmetry: beak vector ==='

    params_beak = 0.0_rk
    params_beak(1) = 2.0_rk
    ! a4 = 0.7495 is the vector the 1.x beak probe settled on, frozen here.
    params_beak(3) = 0.7495_rk
    call probe_tier1_s(params_beak, code)
    call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
            'probe vector is beak-rejected by the tier-1 resolve')

    call compute_neck_standalone_s(params_beak, N_POINTS, ref_z, ref_rho, ref_found, code)
    call assert_int_eq(code, SHAPE_VALID, 'beak vector: tier-1 neck succeeds')
    call assert_true(ref_found, 'beak vector: the tier-1 scanner finds a neck')

    call cache_neck_s(cache, params_beak, z_neck, rho_neck, found, status)
    call assert_int_eq(status, SHAPE_VALID, 'beak vector: neck succeeds (no 103 gate)')
    call assert_true(found, 'beak vector: neck found on the cylindrical grid')
    call assert_close(z_neck, ref_z, TOL_PARITY, 'beak vector: z_neck parity')
    call assert_close(rho_neck, ref_rho, TOL_PARITY, 'beak vector: rho_neck parity')

    z_shift_total = 1.0_rk
    g_opt = 1.0_rk
    call cache_star_convexity_optimum_s(cache, params_beak, z_shift_total, g_opt, status)
    call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, 'beak vector: optimum -> 103')
    call assert_abs_close(z_shift_total, 0.0_rk, 0.0_rk, &
            'beak vector: 103 zero-fills z_shift_total')
    call assert_abs_close(g_opt, 0.0_rk, 0.0_rk, 'beak vector: 103 zero-fills g(s*)')

    radii = 1.0_rk
    call cache_radius_grid_s(cache, params_beak, radii, status)
    call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, 'beak vector: R(theta) -> 103')

    !---------------------------------------------------------------------------
    ! Gating asymmetry, case 2: the optimum reports on a 101-rejected shape
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Gating asymmetry: non-star-convex vector ==='

    ! c = 1.3 is the first elongation the 1.x probe rejected with 101, frozen.
    params_star = ODD_HEAVY
    params_star(1) = 1.3_rk
    call probe_tier1_s(params_star, code)
    call assert_int_eq(code, FOS_ERROR_NOT_STAR_CONVEX, &
            'probe vector is 101-rejected by the tier-1 resolve')

    call cache_neck_s(cache, params_star, z_neck, rho_neck, found, status)
    call assert_int_eq(status, SHAPE_VALID, 'non-star-convex vector: neck succeeds')

    call cache_rho_z_grid_s(cache, params_star, z, rho, drho_dz, grid_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'non-star-convex vector: rho(z) grid succeeds')

    call cache_shape_s(cache, params_star, z_shift, r_north, r_south, status)
    call assert_int_eq(status, FOS_ERROR_NOT_STAR_CONVEX, &
            'non-star-convex vector: shape -> 101')

    call cache_radius_grid_s(cache, params_star, radii, status)
    call assert_int_eq(status, FOS_ERROR_NOT_STAR_CONVEX, &
            'non-star-convex vector: R(theta) -> 101')

    call cache_star_convexity_optimum_s(cache, params_star, z_shift_total, g_opt, status)
    call assert_int_eq(status, SHAPE_VALID, 'non-star-convex vector: optimum still reports')
    call assert_true(g_opt > -STAR_CONVEXITY_MARGIN, &
            'non-star-convex vector: g(s*) is above the acceptance margin')

    call compute_star_convexity_optimum_standalone_s(params_star, N_POINTS, &
            ref_shift_total, ref_g_opt, code)
    call assert_int_eq(code, SHAPE_VALID, 'non-star-convex vector: tier-1 optimum also reports')
    call assert_close(z_shift_total, ref_shift_total, TOL_PARITY, &
            'non-star-convex vector: z_shift_total parity')
    call assert_close(g_opt, ref_g_opt, TOL_PARITY, &
            'non-star-convex vector: g(s*) parity')

    !---------------------------------------------------------------------------
    ! Minimality: the optimum reads the STORED resolve, and nothing else
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Minimality of the optimum ==='

    call cache_shape_s(cache, PARAMS7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'params7 shape valid')
    do i = 1_ik, 8_ik
        before(i) = count_f(i)
    end do
    call cache_star_convexity_optimum_s(cache, PARAMS7, z_shift_total, g_opt, status)
    call assert_int_eq(status, SHAPE_VALID, 'optimum after shape valid')
    do i = 1_ik, 8_ik
        call assert_int_eq(count_f(i), before(i), &
                'optimum after a successful shape recomputes nothing')
    end do

    ! A 101 rejection returns the engine to cold, so the optimum on the SAME
    ! vector has to rebuild #1-#6 from scratch before it can report on it.
    call cache_shape_s(cache, params_star, z_shift, r_north, r_south, status)
    call assert_int_eq(status, FOS_ERROR_NOT_STAR_CONVEX, 'rejected shape -> 101')
    do i = 1_ik, 6_ik
        before(i) = count_f(i)
    end do
    call cache_star_convexity_optimum_s(cache, params_star, z_shift_total, g_opt, status)
    call assert_int_eq(status, SHAPE_VALID, 'optimum after a 101 rejection valid')
    do i = 1_ik, 6_ik
        call assert_int_eq(count_f(i), before(i) + 1_ik, &
                'optimum after a 101 rejection recomputes cold')
    end do

    call cache_free_s(cache)
    call test_summary()

contains

    !> Recompute counter as an `ik` integer, for the assertion helpers.
    function count_f(intermediate) result(n)
        integer(kind = ik), intent(in) :: intermediate
        integer(kind = ik) :: n
        n = int(cache_recompute_count_f(cache, intermediate), ik)
    end function count_f

    !> The tier-1 resolve at the cache's own resolution, for its verdict alone.
    subroutine probe_tier1_s(p, tier1_code)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(out) :: tier1_code

        real(kind = rk) :: tier1_shift, tier1_north, tier1_south

        call compute_shape_standalone_s(p, N_POINTS, tier1_shift, tier1_north, &
                tier1_south, tier1_code)

    end subroutine probe_tier1_s

end program fos_param_extensions_test
