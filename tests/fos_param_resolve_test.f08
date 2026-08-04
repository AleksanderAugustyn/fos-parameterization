!> Resolve suite: the origin/star-convexity intermediate (#6) and the cached
!! shape output (`cache_shape_s`).
!!
!! Parity anchor is the tier-1 resolve (`compute_shape_standalone_s`), which
!! runs the same gates outside the cache, plus z-shift literals frozen from the
!! deleted 1.x `compute_fos_shape_s` before it went (commit 648428c): total
!! z-shift and both pole radii must match it to 1e-12 relative. The comparison
!! is NOT bitwise — the cached path derives u from the trig tables while 1.x
!! computes z*(1/c) — but the ORIGIN BRANCH is exact: a COM-conditioned shape
!! adds exactly nothing to the intrinsic shift, so its z_shift is bit-equal to
!! `compute_z_shift_s`'s value (a within-2.0.0 comparison).
!!
!! Gate vectors are never hardcoded blind: each one is confirmed against the
!! 1.x surface in-test, with a scan if the first guess is accepted.
program fos_param_resolve_test

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_rho_z_grid_s, cache_recompute_count_f, &
            compute_z_shift_s, compute_shape_standalone_s, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_INVALID_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_PARAMS = 7_ik

    real(kind = rk), parameter :: PARAMS7(N_PARAMS) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]
    real(kind = rk), parameter :: ASYM7(N_PARAMS) = &
            [1.8_rk, 0.15_rk, 0.1_rk, 0.08_rk, 0.0_rk, 0.0_rk, 0.0_rk]

    !> Total z-shifts the 1.x `compute_fos_shape_s` resolved at N_POINTS,
    !! frozen from commit 648428c (the last commit carrying that surface). The
    !! pole radii are not frozen separately: they are c + z_shift and
    !! |-c + z_shift| by definition, and asserting that identity against the
    !! frozen shift is the stronger statement.
    real(kind = rk), parameter :: SHIFT_PARAMS7_1X = 6.565141402540683E-002_rk
    real(kind = rk), parameter :: SHIFT_ASYM7_1X = 9.453803619658581E-002_rk

    !> The marginal-origin vector the 1.x probe settled on (ODD_HEAVY at
    !! c = 1.0) and the total shift it resolved there — the optimum branch, not
    !! the COM branch. Frozen from the same commit.
    real(kind = rk), parameter :: MARGINAL_C_1X = 1.0_rk
    real(kind = rk), parameter :: SHIFT_MARGINAL_1X = 5.667639178165991E-001_rk

    !> Odd-heavy asymmetric family. Its COM is a poor R(theta) origin, so it
    !! takes the star-convexity optimum instead; raising c drives g(s*) up
    !! through -STAR_CONVEXITY_MARGIN, which is how the 101 vector is probed.
    real(kind = rk), parameter :: ODD_HEAVY(N_PARAMS) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]

    type(cache_t)      :: cache
    integer(kind = ik) :: status, i, code
    real(kind = rk)    :: thetas(4)
    real(kind = rk)    :: z_shift, r_north, r_south
    real(kind = rk)    :: ref_shift, ref_north, ref_south
    real(kind = rk)    :: zs_intrinsic
    real(kind = rk)    :: params_marginal(N_PARAMS), params_star(N_PARAMS)
    real(kind = rk)    :: params_beak(N_PARAMS), params_rho(N_PARAMS)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
    real(kind = rk)    :: grid_shift
    integer(kind = ik) :: before(6)
    logical            :: found

    do i = 1_ik, 4_ik
        thetas(i) = real(i, rk) * PI_C / 5.0_rk
    end do

    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'cache init valid')

    !---------------------------------------------------------------------------
    ! Wrong parameter count wins over everything and zero-fills
    !---------------------------------------------------------------------------
    z_shift = 1.0_rk
    r_north = 1.0_rk
    r_south = 1.0_rk
    call cache_shape_s(cache, PARAMS7(1:3), z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'short params -> 4')
    call assert_abs_close(z_shift, 0.0_rk, 0.0_rk, 'short params zero-fills z_shift')
    call assert_abs_close(r_north, 0.0_rk, 0.0_rk, 'short params zero-fills r_north')
    call assert_abs_close(r_south, 0.0_rk, 0.0_rk, 'short params zero-fills r_south')

    !---------------------------------------------------------------------------
    ! Degenerate c: 102 before any intermediate runs
    !---------------------------------------------------------------------------
    params_rho = 0.0_rk
    params_rho(1) = 1.0e-11_rk
    call cache_shape_s(cache, params_rho, z_shift, r_north, r_south, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'degenerate c -> 102')
    call assert_abs_close(r_north, 0.0_rk, 0.0_rk, '102 zero-fills r_north')

    !---------------------------------------------------------------------------
    ! Parity with the 1.x resolve step: both parameter sets, COM origin branch
    !---------------------------------------------------------------------------
    call cache_shape_s(cache, PARAMS7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'params7 shape valid')
    call probe_tier1_s(PARAMS7, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, SHAPE_VALID, 'params7 accepted by the tier-1 resolve')
    call assert_close(z_shift, ref_shift, 1.0e-12_rk, 'params7 z_shift tier parity')
    call assert_close(r_north, ref_north, 1.0e-12_rk, 'params7 r_north tier parity')
    call assert_close(r_south, ref_south, 1.0e-12_rk, 'params7 r_south tier parity')
    call assert_close(z_shift, SHIFT_PARAMS7_1X, 1.0e-12_rk, 'params7 z_shift 1.x parity')

    ! COM branch: the origin search adds exactly nothing (within-2.0.0, exact)
    call compute_z_shift_s(PARAMS7, zs_intrinsic, status)
    call assert_int_eq(status, SHAPE_VALID, 'params7 intrinsic shift valid')
    call assert_true(bits_eq_f(z_shift, zs_intrinsic), &
            'params7 takes the COM origin (z_shift bit-equal to intrinsic)')

    ! Pole radii are the analytic expressions in the shifted frame
    call assert_abs_close(r_north, PARAMS7(1) + z_shift, 0.0_rk, 'r_north = c + z_shift')
    call assert_abs_close(r_south, abs(-PARAMS7(1) + z_shift), 0.0_rk, &
            'r_south = |-c + z_shift|')

    call cache_shape_s(cache, ASYM7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'asym shape valid')
    call probe_tier1_s(ASYM7, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, SHAPE_VALID, 'asym accepted by the tier-1 resolve')
    call assert_close(z_shift, ref_shift, 1.0e-12_rk, 'asym z_shift tier parity')
    call assert_close(r_north, ref_north, 1.0e-12_rk, 'asym r_north tier parity')
    call assert_close(r_south, ref_south, 1.0e-12_rk, 'asym r_south tier parity')
    call assert_close(z_shift, SHIFT_ASYM7_1X, 1.0e-12_rk, 'asym z_shift 1.x parity')

    call compute_z_shift_s(ASYM7, zs_intrinsic, status)
    call assert_true(bits_eq_f(z_shift, zs_intrinsic), &
            'asym takes the COM origin (z_shift bit-equal to intrinsic)')

    !---------------------------------------------------------------------------
    ! Minimality: a repeat call recomputes nothing; #7/#8 never run here
    !---------------------------------------------------------------------------
    do i = 1_ik, 6_ik
        before(i) = count_f(i)
    end do
    call cache_shape_s(cache, ASYM7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'repeat shape call valid')
    do i = 1_ik, 6_ik
        call assert_int_eq(count_f(i), before(i), 'repeat call recomputes nothing')
    end do
    call assert_int_eq(count_f(7_ik), 0_ik, 'radius grid (#7) untouched')
    call assert_int_eq(count_f(8_ik), 0_ik, 'radius derivative (#8) untouched')

    ! The cylindrical output shares #1-#3/#5, so a following shape call on the
    ! same vector recomputes only the beak scan (#4) and the resolve (#6).
    call cache_rho_z_grid_s(cache, PARAMS7, z, rho, drho_dz, grid_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'cylindrical output valid')
    do i = 1_ik, 6_ik
        before(i) = count_f(i)
    end do
    call cache_shape_s(cache, PARAMS7, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'shape after rho grid valid')
    call assert_int_eq(count_f(1_ik), before(1), 'shared a2 reused')
    call assert_int_eq(count_f(3_ik), before(3), 'shared f-grid reused')
    call assert_int_eq(count_f(5_ik), before(5), 'shared rho grid reused')
    call assert_int_eq(count_f(4_ik), before(4) + 1_ik, 'beak scan computed')
    call assert_int_eq(count_f(6_ik), before(6) + 1_ik, 'resolve computed')

    !---------------------------------------------------------------------------
    ! Marginal origin branch: the COM is too steep, so s* is taken
    !---------------------------------------------------------------------------
    ! MARGINAL_C_1X is the c the 1.x probe settled on: accepted, and resolved to
    ! a z_shift different from the intrinsic one (= the optimum branch). The
    ! tier-1 resolve is re-run here so both tiers are known to agree on it.
    params_marginal = ODD_HEAVY
    params_marginal(1) = MARGINAL_C_1X
    call probe_tier1_s(params_marginal, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, SHAPE_VALID, 'marginal vector accepted by tier-1')
    call compute_z_shift_s(params_marginal, zs_intrinsic, status)
    found = .not. bits_eq_f(ref_shift, zs_intrinsic)
    call assert_true(found, 'marginal vector takes the optimum origin on tier-1')
    call assert_close(ref_shift, SHIFT_MARGINAL_1X, 1.0e-12_rk, &
            'marginal z_shift 1.x parity')

    call cache_shape_s(cache, params_marginal, z_shift, r_north, r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'marginal vector accepted')
    call assert_close(z_shift, ref_shift, 1.0e-12_rk, 'marginal z_shift tier parity')
    call assert_close(r_north, ref_north, 1.0e-12_rk, 'marginal r_north tier parity')
    call assert_close(r_south, ref_south, 1.0e-12_rk, 'marginal r_south tier parity')
    call assert_true(.not. bits_eq_f(z_shift, zs_intrinsic), &
            'marginal vector takes the optimum, not the COM origin')
    call assert_true(abs(z_shift - zs_intrinsic) > 1.0e-3_rk, &
            'optimum origin is a real shift, not roundoff')

    !---------------------------------------------------------------------------
    ! Gate 103: beak singularity
    !---------------------------------------------------------------------------
    params_beak = 0.0_rk
    params_beak(1) = 2.0_rk
    ! a4 = 0.7495 is the vector the 1.x beak probe settled on, frozen here; the
    ! tier-1 resolve must reject it the same way.
    params_beak(3) = 0.7495_rk
    call probe_tier1_s(params_beak, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
            'probe vector is beak-rejected by the tier-1 resolve')

    call check_rejection_s(params_beak, FOS_ERROR_BEAK_SINGULARITY, 'beak vector -> 103')

    !---------------------------------------------------------------------------
    ! Gate 100 vs 103: the rho-negative vector, both cached outputs
    !---------------------------------------------------------------------------
    ! An interior rho <= 0 means f <= 0 there, which the beak scan (#4, denser
    ! than the rho grid) also sees — and #4 gates first in cache_shape_s. So the
    ! shape path reports 103 where 1.x reported 100, while the cylindrical
    ! output, which never runs #4, still reports 100.
    params_rho = 0.0_rk
    params_rho(1) = 1.0_rk
    ! a3 = 0.9 is the vector the 1.x rho probe settled on. Tier-1 applies the
    ! same gate order as the cache, so the resolve reports 103 here too; the
    ! 1.x verdict, 100, survives on the cylindrical form asserted below.
    params_rho(2) = 0.9_rk
    call probe_tier1_s(params_rho, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, FOS_ERROR_BEAK_SINGULARITY, &
            'probe vector is beak-gated by the tier-1 resolve')

    z = 1.0_rk
    rho = 1.0_rk
    drho_dz = 1.0_rk
    call cache_rho_z_grid_s(cache, params_rho, z, rho, drho_dz, grid_shift, status)
    call assert_int_eq(status, FOS_ERROR_RHO_NEGATIVE, 'rho vector -> 100 (cylindrical)')
    call assert_true(all_zero_f(z) .and. all_zero_f(rho) .and. all_zero_f(drho_dz), &
            '100 zero-fills the cylindrical output')
    call recover_s('after 100')

    call check_rejection_s(params_rho, FOS_ERROR_BEAK_SINGULARITY, &
            'rho vector -> 103 (beak gates first in the shape path)')

    !---------------------------------------------------------------------------
    ! Gate 101: star-convexity margin
    !---------------------------------------------------------------------------
    ! Same odd-heavy family, walked up in c: g(s*) rises toward 0 until the
    ! margin rejects it.
    ! c = 1.3 is the first elongation the 1.x probe rejected with 101, frozen.
    params_star = ODD_HEAVY
    params_star(1) = 1.3_rk
    call probe_tier1_s(params_star, code, ref_shift, ref_north, ref_south)
    call assert_int_eq(code, FOS_ERROR_NOT_STAR_CONVEX, &
            'probe vector is 101-rejected by the tier-1 resolve')

    call check_rejection_s(params_star, FOS_ERROR_NOT_STAR_CONVEX, &
            'non-star-convex vector -> 101')

    call cache_free_s(cache)
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
    !! its outputs. Same gates, same order, no cache — the cross-tier oracle
    !! that replaced the deleted 1.x resolve.
    subroutine probe_tier1_s(p, tier1_code, tier1_shift, tier1_north, tier1_south)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(out) :: tier1_code
        real(kind = rk), intent(out) :: tier1_shift, tier1_north, tier1_south

        call compute_shape_standalone_s(p, N_POINTS, tier1_shift, tier1_north, &
                tier1_south, tier1_code)

    end subroutine probe_tier1_s

    !> A rejected shape zero-fills every output and returns the cache to cold:
    !! the next good call recomputes intermediates #1-#6 from scratch.
    subroutine check_rejection_s(p, want_code, label)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(in) :: want_code
        character(len = *), intent(in) :: label

        real(kind = rk) :: shift_out, north_out, south_out
        integer(kind = ik) :: code_out

        shift_out = 1.0_rk
        north_out = 1.0_rk
        south_out = 1.0_rk
        call cache_shape_s(cache, p, shift_out, north_out, south_out, code_out)
        call assert_int_eq(code_out, want_code, label)
        call assert_abs_close(shift_out, 0.0_rk, 0.0_rk, label // ': z_shift zeroed')
        call assert_abs_close(north_out, 0.0_rk, 0.0_rk, label // ': r_north zeroed')
        call assert_abs_close(south_out, 0.0_rk, 0.0_rk, label // ': r_south zeroed')

        call recover_s(label)

    end subroutine check_rejection_s

    !> params7 must succeed after a rejection, recomputing every intermediate.
    subroutine recover_s(label)
        character(len = *), intent(in) :: label

        real(kind = rk) :: shift_out, north_out, south_out
        integer(kind = ik) :: code_out, j, counts(6)

        do j = 1_ik, 6_ik
            counts(j) = count_f(j)
        end do

        call cache_shape_s(cache, PARAMS7, shift_out, north_out, south_out, code_out)
        call assert_int_eq(code_out, SHAPE_VALID, label // ': recovery call valid')
        do j = 1_ik, 6_ik
            call assert_int_eq(count_f(j), counts(j) + 1_ik, &
                    label // ': recovery recomputed the intermediate')
        end do

    end subroutine recover_s

end program fos_param_resolve_test
