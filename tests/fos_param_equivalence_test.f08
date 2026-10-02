!> Contract family 1: equivalence.
!!
!! One-shot and cached calls must return the same status and bitwise-identical
!! outputs, for any cache with `max_params >= size(params)` at the same
!! resolution; and a short vector must equal its zero-padded form. No tolerance
!! appears anywhere in this suite: every comparison is on the IEEE754 bit
!! patterns, scalars included.
!!
!!   E1  one-shot == cached     seven outputs, vector lengths 1/5/8/50, caches
!!                              with max_params = n, n + 3 and 50, both origin
!!                              branches and a necked shape
!!   E2  rejections             100, 101, 102, 103, 105 and the grid code 3:
!!                              same code and zero-filled outputs in both tiers
!!   E3  short == zero-padded   a length-3 vector against its forms padded to 5
!!                              and 8 with +0.0 and with -0.0, every output of
!!                              both tiers plus a2 and z_shift; and a vector
!!                              whose whole tail is zero against c alone
!!   E4  unchecked profile      a separated shape is 100 checked and valid
!!                              unchecked; a connected one gives the same bits;
!!                              the checked form is a sampled check
!!
!! Every output of one tier is collected into a `result_t` and two results are
!! compared field by field. The at-thetas output has no one-shot form (the
!! one-shot radius forms already take caller thetas), so its fields are compared
!! against the fixed-grid output it must reproduce.
!!
!! Code 104 has no regime here: no known vector reaches it through a gated
!! output. The `standalone` suite covers it through the diagnostic conversion.
program fos_param_equivalence_test

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_shape_s, &
            cache_rho_z_grid_s, cache_rho_z_grid_unchecked_s, cache_neck_s, &
            cache_star_convexity_optimum_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, &
            compute_rho_z_grid_unchecked_standalone_s, compute_neck_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            compute_a2_s, compute_z_shift_s, FOS_MAX_PARAMS, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_BUFFER_MISMATCH
    use test_utils_mod, only: assert_true, assert_int_eq, bits_eq_f, &
            arrays_bits_eq_f, is_zero_f, all_zero_f, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 32_ik

    !> IEEE754 -0.0; built from bits so that -fno-signed-zeros cannot fold it
    !! into +0.0 the way a `-0.0_rk` literal could.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    !> Valid vectors of length 1, 5 and 8. BASE8 keeps the COM origin.
    real(kind = rk), parameter :: P1(1) = [1.3_rk]
    real(kind = rk), parameter :: P5(5) = [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk]
    real(kind = rk), parameter :: BASE8(8) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Odd-heavy family at c = 1.0: accepted, but its COM is too steep an
    !! origin, so the resolve moves to the star-convexity optimum. At c = 1.3
    !! the margin rejects it (101).
    real(kind = rk), parameter :: MARGINAL7(7) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]
    real(kind = rk), parameter :: STAR7(7) = &
            [1.3_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]

    !> Symmetric (c, a4) family: a neck at a4 = 0.4, a beak at a4 = 0.74985
    !! (f(0) = 2.0e-4 < F_MIN_THRESHOLD), two fragments at a4 = 0.9
    !! (f(0) = 1 - 4 a4 / 3 = -0.2).
    real(kind = rk), parameter :: NECK3(3) = [2.0_rk, 0.0_rk, 0.4_rk]
    real(kind = rk), parameter :: BEAK3(3) = [2.0_rk, 0.0_rk, 0.74985_rk]
    real(kind = rk), parameter :: SPLIT3(3) = [2.0_rk, 0.0_rk, 0.9_rk]

    !> Degenerate c (102), and a3 = 0.9: rho <= 0 near a pole (100 on the
    !! cylindrical output, 103 on the R(theta) outputs).
    real(kind = rk), parameter :: BAD_C2(2) = [1.0e-11_rk, 0.1_rk]
    real(kind = rk), parameter :: RHO_NEG2(2) = [1.0_rk, 0.9_rk]

    !> The short vector of E3.
    real(kind = rk), parameter :: SHORT3(3) = [1.6_rk, 0.12_rk, 0.08_rk]

    !> a4 a hair above the touching point 0.75: f(0) = 1 - 4 a4 / 3 = -1.3e-6,
    !! a gap narrower than the node spacing of a 100-point grid.
    real(kind = rk), parameter :: GAP3(3) = [1.0_rk, 0.0_rk, 0.750001_rk]

    !> Every output of one tier for one parameter vector.
    type :: result_t
        integer(kind = ik) :: st_grid = -1_ik, st_deriv = -1_ik, st_at = -1_ik
        integer(kind = ik) :: st_shape = -1_ik, st_rho = -1_ik, st_raw = -1_ik
        integer(kind = ik) :: st_neck = -1_ik, st_opt = -1_ik
        real(kind = rk) :: radii(N_THETA), radii_d(N_THETA), dr(N_THETA)
        real(kind = rk) :: at_radii(N_THETA), at_dr(N_THETA)
        real(kind = rk) :: z_shift, r_north, r_south
        real(kind = rk) :: z(N_POINTS), rho(N_POINTS), drho(N_POINTS), grid_shift
        real(kind = rk) :: raw_z(N_POINTS), raw_rho(N_POINTS), raw_drho(N_POINTS)
        real(kind = rk) :: raw_shift
        real(kind = rk) :: z_neck, rho_neck
        logical :: found = .false.
        real(kind = rk) :: opt_shift, g_opt
    end type result_t

    real(kind = rk) :: thetas(N_THETA)
    real(kind = rk) :: p50(50)
    integer(kind = ik) :: i

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do

    ! A valid 50-parameter vector: every slot nonzero, amplitudes small enough
    ! that the high orders do not fold the surface.
    do i = 2_ik, 50_ik
        p50(i) = 1.0e-6_rk
    end do
    p50(1) = 1.5_rk

    call assert_int_eq(FOS_MAX_PARAMS, 50_ik, 'L = 50')

    call run_e1_s()
    call run_e2_s()
    call run_e3_s()
    call run_e4_s()

    call test_summary()

contains

    !===========================================================================
    ! E1 — one-shot == cached
    !===========================================================================
    subroutine run_e1_s()

        type(result_t) :: res
        real(kind = rk) :: zs
        integer(kind = ik) :: st

        write(*, '(A)') '=== E1: one-shot == cached ==='

        call compare_tiers_s(P1, .true., 'len 1')
        call compare_tiers_s(P5, .true., 'len 5')
        call compare_tiers_s(BASE8, .true., 'len 8')
        call compare_tiers_s(p50, .true., 'len 50')
        call compare_tiers_s(MARGINAL7, .true., 'marginal origin')
        call compare_tiers_s(NECK3, .true., 'necked')

        ! Both origin branches are exercised: the COM origin adds exactly
        ! nothing to the intrinsic shift, the optimum origin a real amount.
        call collect_one_shot_s(BASE8, res)
        call compute_z_shift_s(BASE8, zs, st)
        call assert_int_eq(st, SHAPE_VALID, 'E1: BASE8 intrinsic shift valid')
        call assert_true(bits_eq_f(res%z_shift, zs), 'E1: BASE8 keeps the COM origin')

        call collect_one_shot_s(MARGINAL7, res)
        call compute_z_shift_s(MARGINAL7, zs, st)
        call assert_int_eq(st, SHAPE_VALID, 'E1: MARGINAL7 intrinsic shift valid')
        call assert_true(abs(res%z_shift - zs) > 1.0e-3_rk, &
                'E1: MARGINAL7 moves to the optimum origin')

        ! The necked shape really has a neck, so that output is not compared on
        ! "not found" alone.
        call collect_one_shot_s(NECK3, res)
        call assert_true(res%found, 'E1: NECK3 has a neck')

    end subroutine run_e1_s

    !> One vector, three caches (max_params = n, n + 3 where it fits, and 50):
    !! every cached output equals the one-shot output.
    subroutine compare_tiers_s(params, expect_valid, label)

        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: expect_valid
        character(len = *), intent(in) :: label

        type(result_t) :: one_shot, cached
        type(cache_t) :: cache
        integer(kind = ik) :: sizes(3), k, n, st

        n = size(params, kind = ik)
        sizes(1) = n
        sizes(2) = min(n + 3_ik, FOS_MAX_PARAMS)
        sizes(3) = FOS_MAX_PARAMS

        call collect_one_shot_s(params, one_shot)
        if (expect_valid) then
            call assert_int_eq(one_shot%st_deriv, SHAPE_VALID, &
                    label // ': the vector is valid (comparison is not vacuous)')
        end if

        do k = 1_ik, 3_ik
            call cache_init_s(cache, sizes(k), N_POINTS, thetas, st)
            call assert_int_eq(st, SHAPE_VALID, label // ': cache init')
            call collect_cached_s(cache, params, cached)
            call assert_same_s(cached, one_shot, label // ': cached == one-shot')
            call cache_free_s(cache)
        end do

    end subroutine compare_tiers_s

    !===========================================================================
    ! E2 — rejections
    !===========================================================================
    subroutine run_e2_s()

        type(result_t) :: res
        type(cache_t) :: cache
        real(kind = rk) :: bad_thetas(N_THETA), no_thetas(0)
        real(kind = rk) :: short_r(N_THETA - 1_ik), short_d(N_THETA - 1_ik)
        real(kind = rk) :: none_r(0)
        integer(kind = ik) :: st_cached, st_one_shot

        write(*, '(A)') '=== E2: rejections ==='

        ! Each vector is compared across the tiers on every output, then its
        ! code is pinned on the output that defines it.
        call compare_tiers_s(RHO_NEG2, .false., 'rho-negative vector')
        call collect_one_shot_s(RHO_NEG2, res)
        call assert_int_eq(res%st_rho, FOS_ERROR_RHO_NEGATIVE, 'E2: 100 on the cylindrical output')
        call assert_int_eq(res%st_neck, FOS_ERROR_RHO_NEGATIVE, 'E2: 100 on the neck')
        call assert_int_eq(res%st_shape, FOS_ERROR_BEAK_SINGULARITY, &
                'E2: the same vector is 103 on the shape (beak gates first)')
        call assert_true(all_zero_f(res%rho) .and. all_zero_f(res%z), 'E2: 100 zero-fills')

        call compare_tiers_s(STAR7, .false., 'non-star-convex vector')
        call collect_one_shot_s(STAR7, res)
        call assert_int_eq(res%st_shape, FOS_ERROR_NOT_STAR_CONVEX, 'E2: 101 on the shape')
        call assert_int_eq(res%st_deriv, FOS_ERROR_NOT_STAR_CONVEX, 'E2: 101 on the radii')
        call assert_int_eq(res%st_opt, SHAPE_VALID, 'E2: the optimum still reports on a 101 shape')
        call assert_true(all_zero_f(res%radii_d) .and. all_zero_f(res%dr), 'E2: 101 zero-fills')

        call compare_tiers_s(BAD_C2, .false., 'degenerate c')
        call collect_one_shot_s(BAD_C2, res)
        call assert_int_eq(res%st_grid, FOS_ERROR_INVALID_C, 'E2: 102 on the radius grid')
        call assert_int_eq(res%st_raw, FOS_ERROR_INVALID_C, 'E2: 102 on the unchecked profile')
        call assert_int_eq(res%st_opt, FOS_ERROR_INVALID_C, 'E2: 102 on the optimum')

        call compare_tiers_s(BEAK3, .false., 'beak vector')
        call collect_one_shot_s(BEAK3, res)
        call assert_int_eq(res%st_shape, FOS_ERROR_BEAK_SINGULARITY, 'E2: 103 on the shape')
        call assert_int_eq(res%st_opt, FOS_ERROR_BEAK_SINGULARITY, 'E2: 103 on the optimum')
        call assert_int_eq(res%st_rho, SHAPE_VALID, 'E2: a beak shape still has a profile')

        ! 105: a wrong buffer, in both tiers, ahead of the shape
        call cache_init_s(cache, 8_ik, N_POINTS, thetas, st_cached)
        call assert_int_eq(st_cached, SHAPE_VALID, 'E2: cache init')

        short_r = 1.0_rk
        short_d = 1.0_rk
        call cache_radius_and_derivative_s(cache, BEAK3, short_r, short_d, st_cached)
        call assert_int_eq(st_cached, FOS_ERROR_BUFFER_MISMATCH, 'E2: 105 cached, ahead of 103')
        call assert_true(all_zero_f(short_r) .and. all_zero_f(short_d), 'E2: 105 cached zero-fills')
        short_r = 1.0_rk
        short_d = 1.0_rk
        call compute_radius_and_derivative_standalone_s(BEAK3, thetas, N_POINTS, &
                short_r, short_d, st_one_shot)
        call assert_int_eq(st_one_shot, st_cached, 'E2: 105 one-shot == cached')
        call assert_true(all_zero_f(short_r) .and. all_zero_f(short_d), 'E2: 105 one-shot zero-fills')
        call cache_free_s(cache)

        ! 3: the grid. The cached tier rejects at init, the one-shot at the call.
        call cache_init_s(cache, 8_ik, 99_ik, thetas, st_cached)
        call compute_radius_grid_standalone_s(BASE8, thetas, 99_ik, res%radii, st_one_shot)
        call assert_int_eq(st_cached, SHAPE_ERROR_INVALID_GRID, 'E2: n_points 99, init -> 3')
        call assert_int_eq(st_one_shot, st_cached, 'E2: n_points 99, one-shot == init')
        call assert_true(all_zero_f(res%radii), 'E2: n_points 99 zero-fills')

        call cache_init_s(cache, 8_ik, N_POINTS, no_thetas, st_cached)
        call compute_radius_grid_standalone_s(BASE8, no_thetas, N_POINTS, none_r, st_one_shot)
        call assert_int_eq(st_cached, SHAPE_ERROR_INVALID_GRID, 'E2: empty thetas, init -> 3')
        call assert_int_eq(st_one_shot, st_cached, 'E2: empty thetas, one-shot == init')

        bad_thetas = thetas
        bad_thetas(7) = 1.5_rk * PI_C
        call cache_init_s(cache, 8_ik, N_POINTS, bad_thetas, st_cached)
        call compute_radius_grid_standalone_s(BASE8, bad_thetas, N_POINTS, res%radii, st_one_shot)
        call assert_int_eq(st_cached, SHAPE_ERROR_INVALID_GRID, 'E2: theta = 3pi/2, init -> 3')
        call assert_int_eq(st_one_shot, st_cached, 'E2: theta = 3pi/2, one-shot == init')

    end subroutine run_e2_s

    !===========================================================================
    ! E3 — short == zero-padded
    !===========================================================================
    subroutine run_e3_s()

        type(result_t) :: short_cached, short_one_shot, padded_res
        type(result_t) :: c_only_cached, c_only_one_shot
        type(cache_t) :: cache
        real(kind = rk) :: padded5(5), padded8(8)
        real(kind = rk) :: a2_short, a2_pad, zs_short, zs_pad
        real(kind = rk), volatile :: pad
        integer(kind = ik) :: sign_case, st, j
        character(len = 8) :: tag

        write(*, '(A)') '=== E3: short == zero-padded ==='

        call cache_init_s(cache, 8_ik, N_POINTS, thetas, st)
        call assert_int_eq(st, SHAPE_VALID, 'E3: cache init')

        call collect_cached_s(cache, SHORT3, short_cached)
        call collect_one_shot_s(SHORT3, short_one_shot)
        call assert_int_eq(short_cached%st_deriv, SHAPE_VALID, 'E3: the short vector is valid')
        call compute_a2_s(SHORT3, a2_short, st)
        call compute_z_shift_s(SHORT3, zs_short, st)

        call collect_cached_s(cache, P1, c_only_cached)
        call collect_one_shot_s(P1, c_only_one_shot)

        do sign_case = 1_ik, 2_ik
            if (sign_case == 1_ik) then
                pad = 0.0_rk
                tag = '+0.0 pad'
            else
                pad = transfer(NEG_ZERO_BITS, 1.0_rk)
                tag = '-0.0 pad'
            end if

            do j = 1_ik, 5_ik
                padded5(j) = pad
            end do
            do j = 1_ik, 8_ik
                padded8(j) = pad
            end do
            do j = 1_ik, 3_ik
                padded5(j) = SHORT3(j)
                padded8(j) = SHORT3(j)
            end do
            if (sign_case == 2_ik) then
                call assert_true(transfer(padded8(8), 0_ikl) == NEG_ZERO_BITS, &
                        'E3: the -0.0 padding reached memory')
            end if

            call collect_cached_s(cache, padded5, padded_res)
            call assert_same_s(padded_res, short_cached, 'E3 cached, to 5, ' // tag)
            call collect_cached_s(cache, padded8, padded_res)
            call assert_same_s(padded_res, short_cached, 'E3 cached, to 8, ' // tag)

            call collect_one_shot_s(padded5, padded_res)
            call assert_same_s(padded_res, short_one_shot, 'E3 one-shot, to 5, ' // tag)
            call collect_one_shot_s(padded8, padded_res)
            call assert_same_s(padded_res, short_one_shot, 'E3 one-shot, to 8, ' // tag)

            call compute_a2_s(padded8, a2_pad, st)
            call assert_int_eq(st, SHAPE_VALID, 'E3: a2 of the padded vector valid')
            call assert_true(bits_eq_f(a2_pad, a2_short), 'E3: a2, ' // tag)
            call compute_z_shift_s(padded8, zs_pad, st)
            call assert_int_eq(st, SHAPE_VALID, 'E3: z_shift of the padded vector valid')
            call assert_true(bits_eq_f(zs_pad, zs_short), 'E3: z_shift, ' // tag)

            ! Every coefficient zero: the vector trims down to c alone, and no
            ! further — the trim stops at one parameter.
            do j = 2_ik, 5_ik
                padded5(j) = pad
            end do
            padded5(1) = P1(1)
            call collect_cached_s(cache, padded5, padded_res)
            call assert_same_s(padded_res, c_only_cached, 'E3 cached, all-zero tail, ' // tag)
            call collect_one_shot_s(padded5, padded_res)
            call assert_same_s(padded_res, c_only_one_shot, &
                    'E3 one-shot, all-zero tail, ' // tag)
        end do

        call cache_free_s(cache)

    end subroutine run_e3_s

    !===========================================================================
    ! E4 — unchecked cylindrical profile
    !===========================================================================
    subroutine run_e4_s()

        type(result_t) :: res
        real(kind = rk) :: z100(100), rho100(100), drho100(100), shift, north, south
        integer(kind = ik) :: mid, j, st
        logical :: left_body, right_body

        write(*, '(A)') '=== E4: unchecked profile ==='

        ! Both tiers agree on every output of the separated shape
        call compare_tiers_s(SPLIT3, .false., 'separated shape')

        call collect_one_shot_s(SPLIT3, res)
        call assert_int_eq(res%st_rho, FOS_ERROR_RHO_NEGATIVE, 'E4: checked profile -> 100')
        call assert_true(all_zero_f(res%rho), 'E4: checked profile zero-filled')
        call assert_int_eq(res%st_raw, SHAPE_VALID, 'E4: unchecked profile valid')

        ! The void: rho = 0 and drho/dz = 0 at the middle node, a body on each side
        mid = (N_POINTS + 1_ik) / 2_ik
        call assert_true(is_zero_f(res%raw_rho(mid)), 'E4: rho = 0 in the void')
        call assert_true(is_zero_f(res%raw_drho(mid)), 'E4: drho/dz = 0 in the void')
        left_body = .false.
        right_body = .false.
        do j = 2_ik, mid - 1_ik
            if (res%raw_rho(j) > 0.1_rk) left_body = .true.
        end do
        do j = mid + 1_ik, N_POINTS - 1_ik
            if (res%raw_rho(j) > 0.1_rk) right_body = .true.
        end do
        call assert_true(left_body .and. right_body, 'E4: a fragment on each side of the void')

        ! A connected shape: checked and unchecked are the same bits
        call collect_one_shot_s(BASE8, res)
        call assert_int_eq(res%st_rho, SHAPE_VALID, 'E4: connected shape, checked valid')
        call assert_int_eq(res%st_raw, SHAPE_VALID, 'E4: connected shape, unchecked valid')
        call assert_true(arrays_bits_eq_f(res%raw_z, res%z), 'E4: unchecked z == checked z')
        call assert_true(arrays_bits_eq_f(res%raw_rho, res%rho), 'E4: unchecked rho == checked rho')
        call assert_true(arrays_bits_eq_f(res%raw_drho, res%drho), &
                'E4: unchecked drho/dz == checked drho/dz')
        call assert_true(bits_eq_f(res%raw_shift, res%grid_shift), &
                'E4: unchecked z_shift == checked z_shift')

        ! The checked profile is a SAMPLED check, not a connectivity verdict: a
        ! 100-point grid has no node at u = 0, where GAP3 has f < 0, so the gap
        ! passes. The R(theta) outputs see it through the 1001-point beak scan.
        call compute_rho_z_grid_standalone_s(GAP3, 100_ik, z100, rho100, drho100, &
                shift, st)
        call assert_int_eq(st, SHAPE_VALID, &
                'E4: a gap between two nodes passes the checked profile')
        call compute_shape_standalone_s(GAP3, 100_ik, shift, north, south, st)
        call assert_int_eq(st, FOS_ERROR_BEAK_SINGULARITY, &
                'E4: the same shape is 103 on the R(theta) side')

    end subroutine run_e4_s

    !===========================================================================
    ! Collectors and comparison
    !===========================================================================

    !> Every cached output of `params` on `cache`.
    subroutine collect_cached_s(cache, params, res)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res

        call cache_radius_grid_s(cache, params, res%radii, res%st_grid)
        call cache_radius_and_derivative_s(cache, params, res%radii_d, res%dr, res%st_deriv)
        call cache_radius_and_derivative_at_thetas_s(cache, params, thetas, &
                res%at_radii, res%at_dr, res%st_at)
        call cache_shape_s(cache, params, res%z_shift, res%r_north, res%r_south, res%st_shape)
        call cache_rho_z_grid_s(cache, params, res%z, res%rho, res%drho, &
                res%grid_shift, res%st_rho)
        call cache_rho_z_grid_unchecked_s(cache, params, res%raw_z, res%raw_rho, &
                res%raw_drho, res%raw_shift, res%st_raw)
        call cache_neck_s(cache, params, res%z_neck, res%rho_neck, res%found, res%st_neck)
        call cache_star_convexity_optimum_s(cache, params, res%opt_shift, res%g_opt, &
                res%st_opt)

    end subroutine collect_cached_s

    !> Every one-shot output of `params`. The at-thetas fields are filled by a
    !! second radius-and-derivative call: there is no one-shot at-thetas form.
    subroutine collect_one_shot_s(params, res)

        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res

        call compute_radius_grid_standalone_s(params, thetas, N_POINTS, res%radii, &
                res%st_grid)
        call compute_radius_and_derivative_standalone_s(params, thetas, N_POINTS, &
                res%radii_d, res%dr, res%st_deriv)
        call compute_radius_and_derivative_standalone_s(params, thetas, N_POINTS, &
                res%at_radii, res%at_dr, res%st_at)
        call compute_shape_standalone_s(params, N_POINTS, res%z_shift, res%r_north, &
                res%r_south, res%st_shape)
        call compute_rho_z_grid_standalone_s(params, N_POINTS, res%z, res%rho, &
                res%drho, res%grid_shift, res%st_rho)
        call compute_rho_z_grid_unchecked_standalone_s(params, N_POINTS, res%raw_z, &
                res%raw_rho, res%raw_drho, res%raw_shift, res%st_raw)
        call compute_neck_standalone_s(params, N_POINTS, res%z_neck, res%rho_neck, &
                res%found, res%st_neck)
        call compute_star_convexity_optimum_standalone_s(params, N_POINTS, &
                res%opt_shift, res%g_opt, res%st_opt)

    end subroutine collect_one_shot_s

    !> Two results are the same: every status equal, every output bit-identical.
    !! Also asserts the two internal identities of one result: the radii of the
    !! two fixed-grid outputs, and the at-thetas output against the fixed grid.
    subroutine assert_same_s(a, b, label)

        type(result_t), intent(in) :: a, b
        character(len = *), intent(in) :: label

        call assert_int_eq(a%st_grid, b%st_grid, label // ': radius_grid status')
        call assert_int_eq(a%st_deriv, b%st_deriv, label // ': radius_and_derivative status')
        call assert_int_eq(a%st_at, b%st_at, label // ': at-thetas status')
        call assert_int_eq(a%st_shape, b%st_shape, label // ': shape status')
        call assert_int_eq(a%st_rho, b%st_rho, label // ': rho_z_grid status')
        call assert_int_eq(a%st_raw, b%st_raw, label // ': rho_z_grid_unchecked status')
        call assert_int_eq(a%st_neck, b%st_neck, label // ': neck status')
        call assert_int_eq(a%st_opt, b%st_opt, label // ': optimum status')

        call assert_true(arrays_bits_eq_f(a%radii, b%radii), label // ': radius_grid radii')
        call assert_true(arrays_bits_eq_f(a%radii_d, b%radii_d), &
                label // ': radius_and_derivative radii')
        call assert_true(arrays_bits_eq_f(a%dr, b%dr), label // ': dR/dtheta')
        call assert_true(arrays_bits_eq_f(a%at_radii, b%at_radii), label // ': at-thetas radii')
        call assert_true(arrays_bits_eq_f(a%at_dr, b%at_dr), label // ': at-thetas dR/dtheta')

        call assert_true(bits_eq_f(a%z_shift, b%z_shift) .and. bits_eq_f(a%r_north, b%r_north) &
                .and. bits_eq_f(a%r_south, b%r_south), label // ': shape scalars')

        call assert_true(arrays_bits_eq_f(a%z, b%z) .and. arrays_bits_eq_f(a%rho, b%rho) &
                .and. arrays_bits_eq_f(a%drho, b%drho) &
                .and. bits_eq_f(a%grid_shift, b%grid_shift), label // ': rho_z_grid')
        call assert_true(arrays_bits_eq_f(a%raw_z, b%raw_z) &
                .and. arrays_bits_eq_f(a%raw_rho, b%raw_rho) &
                .and. arrays_bits_eq_f(a%raw_drho, b%raw_drho) &
                .and. bits_eq_f(a%raw_shift, b%raw_shift), label // ': rho_z_grid_unchecked')

        call assert_true((a%found .eqv. b%found) .and. bits_eq_f(a%z_neck, b%z_neck) &
                .and. bits_eq_f(a%rho_neck, b%rho_neck), label // ': neck')
        call assert_true(bits_eq_f(a%opt_shift, b%opt_shift) .and. bits_eq_f(a%g_opt, b%g_opt), &
                label // ': optimum')

        ! Within one result
        call assert_true(arrays_bits_eq_f(a%radii, a%radii_d), &
                label // ': radius_grid == radius_and_derivative radii')
        call assert_true(arrays_bits_eq_f(a%at_radii, a%radii_d) &
                .and. arrays_bits_eq_f(a%at_dr, a%dr), &
                label // ': at-thetas on the primary thetas == fixed grid')

    end subroutine assert_same_s

end program fos_param_equivalence_test
