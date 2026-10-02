!> Contract family 2: statelessness.
!!
!! A cache stores nothing derived from a shape parameter, so no call can
!! influence a later one. Three properties, all bitwise:
!!
!!   S1  a repeated call returns the same outputs
!!   S2  a call after a rejected call (4, 105, 3, 100, 101, 102, 103) returns
!!       the same outputs as before the rejection
!!   S3  any sequence of calls on one shared cache — different outputs,
!!       different vector lengths, different theta sets, rejections in between
!!       — returns for each call the outputs that call gives alone on a fresh
!!       cache
!!
!! S3 is also the library's cross-call-site check. The shared-cache calls are
!! made in `shared_forward_s` / `shared_reverse_s`; the fresh-cache reference
!! is built and called in `reference_s`. Under `-flto=auto -ffast-math` a
!! kernel inlined into each of those procedures is optimized per call site and
!! can differ in the last bits. The library compiles its own objects without
!! LTO and keeps the arithmetic in a different file from the tiers, so every
!! path here executes the same machine code; this suite is what fails if that
!! ever stops being true.
program fos_param_statelessness_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, cache_radius_and_derivative_s, &
            cache_radius_and_derivative_at_thetas_s, cache_shape_s, &
            cache_rho_z_grid_s, cache_rho_z_grid_unchecked_s, cache_neck_s, &
            cache_star_convexity_optimum_s, &
            SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, SHAPE_ERROR_WRONG_PARAM_COUNT, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_BUFFER_MISMATCH
    use test_utils_mod, only: assert_true, assert_int_eq, bits_eq_f, &
            arrays_bits_eq_f, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 301_ik
    integer(kind = ik), parameter :: N_THETA = 24_ik
    integer(kind = ik), parameter :: N_OTHER = 7_ik
    integer(kind = ik), parameter :: MAX_PARAMS = 8_ik

    real(kind = rk), parameter :: P1(1) = [1.3_rk]
    real(kind = rk), parameter :: SHORT3(3) = [1.6_rk, 0.12_rk, 0.08_rk]
    real(kind = rk), parameter :: NECK3(3) = [2.0_rk, 0.0_rk, 0.4_rk]
    real(kind = rk), parameter :: BASE8(8) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]
    real(kind = rk), parameter :: MARGINAL7(7) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]

    ! Rejected vectors: 101, 103, 102, and rho <= 0 (100 cylindrical, 103 R(theta))
    real(kind = rk), parameter :: STAR7(7) = &
            [1.3_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk]
    real(kind = rk), parameter :: BEAK3(3) = [2.0_rk, 0.0_rk, 0.74985_rk]
    real(kind = rk), parameter :: BAD_C2(2) = [1.0e-11_rk, 0.1_rk]
    real(kind = rk), parameter :: RHO_NEG2(2) = [1.0_rk, 0.9_rk]
    real(kind = rk), parameter :: LONG9(MAX_PARAMS + 1_ik) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk, &
             0.001_rk, 0.0005_rk]

    !> Every cached output for one parameter vector, the at-thetas output on a
    !! second theta set.
    type :: result_t
        integer(kind = ik) :: st(8) = -1_ik
        real(kind = rk) :: radii(N_THETA), radii_d(N_THETA), dr(N_THETA)
        real(kind = rk) :: at_radii(N_OTHER), at_dr(N_OTHER)
        real(kind = rk) :: z_shift, r_north, r_south
        real(kind = rk) :: z(N_POINTS), rho(N_POINTS), drho(N_POINTS), grid_shift
        real(kind = rk) :: raw_z(N_POINTS), raw_rho(N_POINTS), raw_drho(N_POINTS)
        real(kind = rk) :: raw_shift
        real(kind = rk) :: z_neck, rho_neck
        logical :: found = .false.
        real(kind = rk) :: opt_shift, g_opt
    end type result_t

    type(cache_t) :: shared
    type(result_t) :: first, again, ref
    real(kind = rk) :: thetas(N_THETA), other_thetas(N_OTHER)
    real(kind = rk) :: bad_thetas(N_OTHER), no_thetas(0), no_r(0), no_d(0)
    real(kind = rk) :: radii(N_THETA), short_radii(N_THETA - 1_ik)
    real(kind = rk) :: at_r(N_OTHER), at_d(N_OTHER)
    real(kind = rk) :: z(N_POINTS), rho(N_POINTS), drho(N_POINTS), shift, north, south
    integer(kind = ik) :: i, status, pass

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do
    ! The second theta set includes both poles
    do i = 1_ik, N_OTHER
        other_thetas(i) = real(i - 1_ik, rk) * PI_C / real(N_OTHER - 1_ik, rk)
    end do
    other_thetas(1) = 0.0_rk
    other_thetas(N_OTHER) = PI_C

    call cache_init_s(shared, MAX_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared cache init')

    !---------------------------------------------------------------------------
    ! S1 — a repeated call
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== S1: repeated calls ==='

    call shared_forward_s(shared, BASE8, first)
    call assert_true(all(first%st == SHAPE_VALID), 'S1: BASE8 valid on every output')
    call shared_forward_s(shared, BASE8, again)
    call assert_same_s(again, first, 'S1: second call')
    call shared_reverse_s(shared, BASE8, again)
    call assert_same_s(again, first, 'S1: third call, outputs in reverse order')

    !---------------------------------------------------------------------------
    ! S2 — a call after each rejection
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== S2: after a rejected call ==='

    call cache_radius_grid_s(shared, LONG9, radii, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'S2: rejection 4')
    call after_rejection_s('S2 after 4')

    call cache_radius_grid_s(shared, BASE8, short_radii, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'S2: rejection 105')
    call after_rejection_s('S2 after 105')

    bad_thetas = other_thetas
    bad_thetas(3) = 1.5_rk * PI_C
    call cache_radius_and_derivative_at_thetas_s(shared, BASE8, bad_thetas, at_r, &
            at_d, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'S2: rejection 3 (theta range)')
    call after_rejection_s('S2 after 3 (theta range)')

    call cache_radius_and_derivative_at_thetas_s(shared, BASE8, no_thetas, no_r, &
            no_d, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'S2: rejection 3 (empty thetas)')
    call after_rejection_s('S2 after 3 (empty thetas)')

    call cache_rho_z_grid_s(shared, RHO_NEG2, z, rho, drho, shift, status)
    call assert_int_eq(status, FOS_ERROR_RHO_NEGATIVE, 'S2: rejection 100')
    call after_rejection_s('S2 after 100')

    call cache_shape_s(shared, STAR7, shift, north, south, status)
    call assert_int_eq(status, FOS_ERROR_NOT_STAR_CONVEX, 'S2: rejection 101')
    call after_rejection_s('S2 after 101')

    call cache_radius_grid_s(shared, BAD_C2, radii, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'S2: rejection 102')
    call after_rejection_s('S2 after 102')

    call cache_radius_grid_s(shared, BEAK3, radii, status)
    call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, 'S2: rejection 103')
    call after_rejection_s('S2 after 103')

    !---------------------------------------------------------------------------
    ! S3 — mixed sequences on the shared cache == each call alone on a fresh
    ! cache, the reference computed in another procedure
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== S3: mixed sequences vs fresh caches ==='

    do pass = 1_ik, 2_ik
        call mixed_s(BASE8, pass, 'BASE8')
        call mixed_s(STAR7, pass, 'STAR7 (101)')
        call mixed_s(P1, pass, 'P1')
        call mixed_s(BEAK3, pass, 'BEAK3 (103)')
        call mixed_s(MARGINAL7, pass, 'MARGINAL7')
        call mixed_s(BAD_C2, pass, 'BAD_C2 (102)')
        call mixed_s(NECK3, pass, 'NECK3')
        call mixed_s(RHO_NEG2, pass, 'RHO_NEG2 (100)')
        call mixed_s(SHORT3, pass, 'SHORT3')
        call mixed_s(BASE8, pass, 'BASE8 again')
    end do

    call cache_free_s(shared)
    call test_summary()

contains

    !> BASE8 on the shared cache must still equal its first result.
    subroutine after_rejection_s(label)

        character(len = *), intent(in) :: label

        type(result_t) :: now

        call shared_forward_s(shared, BASE8, now)
        call assert_same_s(now, first, label)

    end subroutine after_rejection_s

    !> One step of the mixed sequence: the shared-cache result (outputs in
    !! forward order on pass 1, reverse order on pass 2) against the fresh-cache
    !! reference.
    subroutine mixed_s(params, which_pass, label)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: which_pass
        character(len = *), intent(in) :: label

        type(result_t) :: got

        call reference_s(params, ref)
        if (which_pass == 1_ik) then
            call shared_forward_s(shared, params, got)
            call assert_same_s(got, ref, 'S3 pass 1, ' // label)
        else
            call shared_reverse_s(shared, params, got)
            call assert_same_s(got, ref, 'S3 pass 2, ' // label)
        end if

    end subroutine mixed_s

    !> All eight outputs on a caller's cache, in declaration order.
    subroutine shared_forward_s(cache, params, res)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res

        call cache_radius_grid_s(cache, params, res%radii, res%st(1))
        call cache_radius_and_derivative_s(cache, params, res%radii_d, res%dr, res%st(2))
        call cache_radius_and_derivative_at_thetas_s(cache, params, other_thetas, &
                res%at_radii, res%at_dr, res%st(3))
        call cache_shape_s(cache, params, res%z_shift, res%r_north, res%r_south, res%st(4))
        call cache_rho_z_grid_s(cache, params, res%z, res%rho, res%drho, &
                res%grid_shift, res%st(5))
        call cache_rho_z_grid_unchecked_s(cache, params, res%raw_z, res%raw_rho, &
                res%raw_drho, res%raw_shift, res%st(6))
        call cache_neck_s(cache, params, res%z_neck, res%rho_neck, res%found, res%st(7))
        call cache_star_convexity_optimum_s(cache, params, res%opt_shift, res%g_opt, &
                res%st(8))

    end subroutine shared_forward_s

    !> The same eight outputs in the opposite order.
    subroutine shared_reverse_s(cache, params, res)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res

        call cache_star_convexity_optimum_s(cache, params, res%opt_shift, res%g_opt, &
                res%st(8))
        call cache_neck_s(cache, params, res%z_neck, res%rho_neck, res%found, res%st(7))
        call cache_rho_z_grid_unchecked_s(cache, params, res%raw_z, res%raw_rho, &
                res%raw_drho, res%raw_shift, res%st(6))
        call cache_rho_z_grid_s(cache, params, res%z, res%rho, res%drho, &
                res%grid_shift, res%st(5))
        call cache_shape_s(cache, params, res%z_shift, res%r_north, res%r_south, res%st(4))
        call cache_radius_and_derivative_at_thetas_s(cache, params, other_thetas, &
                res%at_radii, res%at_dr, res%st(3))
        call cache_radius_and_derivative_s(cache, params, res%radii_d, res%dr, res%st(2))
        call cache_radius_grid_s(cache, params, res%radii, res%st(1))

    end subroutine shared_reverse_s

    !> Each output alone on its own fresh cache, built and freed here. A
    !! different procedure from the shared-cache callers on purpose: the
    !! comparison must cross call sites.
    subroutine reference_s(params, res)

        real(kind = rk), intent(in) :: params(:)
        type(result_t), intent(out) :: res

        type(cache_t) :: fresh
        integer(kind = ik) :: k, st

        do k = 1_ik, 8_ik
            call cache_init_s(fresh, MAX_PARAMS, N_POINTS, thetas, st)
            call assert_int_eq(st, SHAPE_VALID, 'reference cache init')
            select case (k)
            case (1_ik)
                call cache_radius_grid_s(fresh, params, res%radii, res%st(1))
            case (2_ik)
                call cache_radius_and_derivative_s(fresh, params, res%radii_d, &
                        res%dr, res%st(2))
            case (3_ik)
                call cache_radius_and_derivative_at_thetas_s(fresh, params, &
                        other_thetas, res%at_radii, res%at_dr, res%st(3))
            case (4_ik)
                call cache_shape_s(fresh, params, res%z_shift, res%r_north, &
                        res%r_south, res%st(4))
            case (5_ik)
                call cache_rho_z_grid_s(fresh, params, res%z, res%rho, res%drho, &
                        res%grid_shift, res%st(5))
            case (6_ik)
                call cache_rho_z_grid_unchecked_s(fresh, params, res%raw_z, &
                        res%raw_rho, res%raw_drho, res%raw_shift, res%st(6))
            case (7_ik)
                call cache_neck_s(fresh, params, res%z_neck, res%rho_neck, &
                        res%found, res%st(7))
            case (8_ik)
                call cache_star_convexity_optimum_s(fresh, params, res%opt_shift, &
                        res%g_opt, res%st(8))
            end select
            call cache_free_s(fresh)
        end do

    end subroutine reference_s

    !> Two results are the same: every status equal, every output bit-identical.
    subroutine assert_same_s(a, b, label)

        type(result_t), intent(in) :: a, b
        character(len = *), intent(in) :: label

        integer(kind = ik) :: k

        do k = 1_ik, 8_ik
            call assert_int_eq(a%st(k), b%st(k), label // ': status')
        end do

        call assert_true(arrays_bits_eq_f(a%radii, b%radii), label // ': radius_grid')
        call assert_true(arrays_bits_eq_f(a%radii_d, b%radii_d) &
                .and. arrays_bits_eq_f(a%dr, b%dr), label // ': radius_and_derivative')
        call assert_true(arrays_bits_eq_f(a%at_radii, b%at_radii) &
                .and. arrays_bits_eq_f(a%at_dr, b%at_dr), label // ': at-thetas')
        call assert_true(bits_eq_f(a%z_shift, b%z_shift) .and. bits_eq_f(a%r_north, b%r_north) &
                .and. bits_eq_f(a%r_south, b%r_south), label // ': shape')
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

    end subroutine assert_same_s

end program fos_param_statelessness_test
