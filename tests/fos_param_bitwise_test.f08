!> Contract family 1: bitwise warm == cold.
!!
!! An incremental compute — a cache walked through a parameter sweep, reusing
!! every intermediate its diff did not invalidate — must reproduce a COLD
!! compute (fresh cache, same init arguments, same final parameter vector) bit
!! for bit. No tolerance appears anywhere in this suite; every comparison is on
!! the IEEE754 bit patterns.
!!
!! The sweep runs over FOUR base regimes, so the property is pinned on both
!! origin branches and on both sides of the shape gates:
!!
!!   (a) symmetric      — every odd coefficient zero, intrinsic shift 0
!!   (b) asymmetric     — the full 8-parameter base
!!   (c) near-beak      — a4 probed upward to the last value before 103
!!   (d) marginal-origin— the odd-heavy family, whose COM is too steep an origin
!!
!! Branch detection is the z_shift itself: a base takes the COM origin iff
!! `cache_shape_s`'s z_shift is bit-equal to `compute_z_shift_s`'s intrinsic
!! value (origin s = 0 adds exactly nothing), and the optimum origin otherwise.
!! `cache_star_convexity_optimum_s` cannot answer this — it ALWAYS reports the
!! optimum, by contract, so it cannot distinguish the branch actually taken.
!!
!! MEASURED, not assumed: (a) and (b) are COM-branch, (c) and (d) are
!! optimum-branch. A near-beak shape is a deeply necked one, and a deep neck is
!! exactly what makes the COM a badly conditioned origin (g(0) above
!! -ORIGIN_CONDITION_MARGIN), so the near-beak regime routes to the optimum. The
!! contract's requirement — both branches demonstrably exercised — is met by
!! (a)/(b) versus (c)/(d), and this suite asserts the branch each regime is
!! measured to take rather than a guess about it.
!!
!! Rejections are part of the property: warm and cold must return the SAME
!! status, and a rejected call zero-fills, so the bit comparison covers the
!! rejected steps too. One deliberate mid-sweep beak rejection then a resumed
!! valid step pins the recovery path.
program fos_param_bitwise_test

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_radius_grid_s, compute_z_shift_s, &
            SHAPE_VALID, FOS_ERROR_BEAK_SINGULARITY
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            test_summary

    implicit none

    integer(kind = ik), parameter :: N_PARAMS = 8_ik
    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 64_ik
    integer(kind = ik), parameter :: N_REGIMES = 4_ik

    !> Sweep steps: five values of the swept parameter, base + k*delta,
    !! k = -2..2. c moves on a coarser step than the coefficients.
    real(kind = rk), parameter :: DELTA_COEFF = 0.02_rk
    real(kind = rk), parameter :: DELTA_C = 0.1_rk

    real(kind = rk), parameter :: BASE8(N_PARAMS) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Odd-heavy family whose COM is a poor R(theta) origin — the marginal-origin
    !! regime's probe base (same family as the resolve suite's).
    real(kind = rk), parameter :: ODD_HEAVY(N_PARAMS) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk, 0.0_rk]

    character(len = 16), parameter :: REGIME_NAMES(N_REGIMES) = &
            ['symmetric       ', 'asymmetric      ', 'near-beak       ', &
             'marginal-origin ']

    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: bases(N_PARAMS, N_REGIMES)
    real(kind = rk)    :: beak_params(N_PARAMS)
    integer(kind = ik) :: i, r
    logical            :: com_branch(N_REGIMES)

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do

    !---------------------------------------------------------------------------
    ! Regime bases. (c) and (d) are PROBED against the surface that defines
    ! them; nothing here is a hardcoded guess.
    !---------------------------------------------------------------------------
    do i = 1_ik, N_PARAMS
        bases(i, 1) = BASE8(i)
        bases(i, 2) = BASE8(i)
    end do
    ! (a) symmetric: every odd coefficient (a3, a5, a7, a9 = params 2, 4, 6, 8) zero
    bases(2, 1) = 0.0_rk
    bases(4, 1) = 0.0_rk
    bases(6, 1) = 0.0_rk
    bases(8, 1) = 0.0_rk

    call probe_near_beak_s(bases(:, 1), bases(:, 3), beak_params)
    call probe_marginal_s(bases(:, 4))

    !---------------------------------------------------------------------------
    ! Origin branch of every regime, and the coverage requirement
    !---------------------------------------------------------------------------
    do r = 1_ik, N_REGIMES
        call classify_branch_s(bases(:, r), com_branch(r), trim(REGIME_NAMES(r)))
    end do
    call assert_true(com_branch(1), 'symmetric base takes the COM origin')
    call assert_true(com_branch(2), 'asymmetric base takes the COM origin')
    call assert_true(.not. com_branch(3), 'near-beak base takes the optimum origin')
    call assert_true(.not. com_branch(4), 'marginal base takes the optimum origin')
    call assert_true(any(com_branch) .and. .not. all(com_branch), &
            'the sweep exercises both origin branches')

    !---------------------------------------------------------------------------
    ! The sweep itself
    !---------------------------------------------------------------------------
    do r = 1_ik, N_REGIMES
        call run_sweep_s(bases(:, r), trim(REGIME_NAMES(r)))
    end do

    call run_invalid_resume_s(bases(:, 2), beak_params)

    call test_summary()

contains

    !> Last a4 (params(3)) the shape gates still accept before the beak
    !! rejection fires, walked up from the symmetric base. Also hands back the
    !! first REJECTED vector, which the mid-sweep recovery case needs.
    subroutine probe_near_beak_s(sym_base, near_beak, rejected)

        real(kind = rk), intent(in) :: sym_base(N_PARAMS)
        real(kind = rk), intent(out) :: near_beak(N_PARAMS)
        real(kind = rk), intent(out) :: rejected(N_PARAMS)

        type(cache_t) :: probe
        real(kind = rk) :: p(N_PARAMS), z_shift, r_north, r_south
        integer(kind = ik) :: j, k, status
        logical :: hit

        call cache_init_s(probe, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'near-beak probe cache init')

        do j = 1_ik, N_PARAMS
            near_beak(j) = sym_base(j)
            rejected(j) = sym_base(j)
        end do

        hit = .false.
        do k = 0_ik, 300_ik
            do j = 1_ik, N_PARAMS
                p(j) = sym_base(j)
            end do
            p(3) = sym_base(3) + 5.0e-3_rk * real(k, rk)
            call cache_shape_s(probe, p, z_shift, r_north, r_south, status)
            if (status == SHAPE_VALID) then
                do j = 1_ik, N_PARAMS
                    near_beak(j) = p(j)
                end do
            else
                call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, &
                        'the a4 probe leaves the valid region through the beak gate')
                do j = 1_ik, N_PARAMS
                    rejected(j) = p(j)
                end do
                hit = .true.
                exit
            end if
        end do

        call assert_true(hit, 'a4 probe reached the beak rejection')
        call assert_true(near_beak(3) > sym_base(3), &
                'near-beak base is a genuinely necked shape (a4 raised)')

        call cache_free_s(probe)

    end subroutine probe_near_beak_s

    !> First vector of the odd-heavy family that is accepted AND whose resolved
    !! z_shift differs bitwise from the intrinsic one — i.e. the origin search
    !! moved it off the COM, which is the definition of the optimum branch.
    subroutine probe_marginal_s(marginal)

        real(kind = rk), intent(out) :: marginal(N_PARAMS)

        type(cache_t) :: probe
        real(kind = rk) :: p(N_PARAMS), z_shift, r_north, r_south, zs_intrinsic
        integer(kind = ik) :: j, k, status, zstatus
        logical :: hit

        call cache_init_s(probe, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'marginal probe cache init')

        do j = 1_ik, N_PARAMS
            marginal(j) = ODD_HEAVY(j)
        end do

        hit = .false.
        do k = 0_ik, 20_ik
            do j = 1_ik, N_PARAMS
                p(j) = ODD_HEAVY(j)
            end do
            p(1) = ODD_HEAVY(1) - 0.025_rk * real(k, rk)
            call cache_shape_s(probe, p, z_shift, r_north, r_south, status)
            if (status /= SHAPE_VALID) cycle
            call compute_z_shift_s(p, zs_intrinsic, zstatus)
            if (zstatus /= SHAPE_VALID) cycle
            if (.not. bits_eq_f(z_shift, zs_intrinsic)) then
                do j = 1_ik, N_PARAMS
                    marginal(j) = p(j)
                end do
                call assert_true(abs(z_shift - zs_intrinsic) > 1.0e-3_rk, &
                        'marginal probe found a real shift, not roundoff')
                hit = .true.
                exit
            end if
        end do

        call assert_true(hit, 'probe found a marginal-origin vector')

        call cache_free_s(probe)

    end subroutine probe_marginal_s

    !> Origin branch of one base: COM iff the resolved z_shift is bit-equal to
    !! the intrinsic COM shift.
    subroutine classify_branch_s(base, is_com, label)

        real(kind = rk), intent(in) :: base(N_PARAMS)
        logical, intent(out) :: is_com
        character(len = *), intent(in) :: label

        type(cache_t) :: probe
        real(kind = rk) :: z_shift, r_north, r_south, zs_intrinsic
        integer(kind = ik) :: status, zstatus

        call cache_init_s(probe, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, label // ': branch probe init')

        call cache_shape_s(probe, base, z_shift, r_north, r_south, status)
        call assert_int_eq(status, SHAPE_VALID, label // ': base vector is accepted')
        call compute_z_shift_s(base, zs_intrinsic, zstatus)
        call assert_int_eq(zstatus, SHAPE_VALID, label // ': intrinsic shift valid')

        is_com = bits_eq_f(z_shift, zs_intrinsic)

        call cache_free_s(probe)

    end subroutine classify_branch_s

    !> One regime: ONE warm cache walked through all eight parameters, five
    !! values each, against a fresh cold cache at every step.
    subroutine run_sweep_s(base, label)

        real(kind = rk), intent(in) :: base(N_PARAMS)
        character(len = *), intent(in) :: label

        type(cache_t) :: warm
        real(kind = rk) :: params(N_PARAMS), delta
        integer(kind = ik) :: j, k, m, status
        character(len = 64) :: step_label

        call cache_init_s(warm, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, label // ': warm cache init')

        do j = 1_ik, N_PARAMS
            delta = DELTA_COEFF
            if (j == 1_ik) delta = DELTA_C

            do k = -2_ik, 2_ik
                do m = 1_ik, N_PARAMS
                    params(m) = base(m)
                end do
                params(j) = base(j) + real(k, rk) * delta

                write(step_label, '(A,A,I0,A,I0)') label, ', param ', j, ', step ', k
                call compare_warm_cold_s(warm, params, trim(step_label))
            end do
        end do

        call cache_free_s(warm)

    end subroutine run_sweep_s

    !> The property itself: the warm cache's radius grid and resolved shape
    !! must match a cold cache's, status and bits.
    subroutine compare_warm_cold_s(warm, params, label)

        type(cache_t), intent(inout) :: warm
        real(kind = rk), intent(in) :: params(N_PARAMS)
        character(len = *), intent(in) :: label

        type(cache_t) :: cold
        real(kind = rk) :: radii_w(N_THETA), radii_c(N_THETA)
        real(kind = rk) :: zs_w, rn_w, rs_w, zs_c, rn_c, rs_c
        integer(kind = ik) :: grid_w, grid_c, shape_w, shape_c, init_c, n

        call cache_radius_grid_s(warm, params, radii_w, grid_w)
        call cache_shape_s(warm, params, zs_w, rn_w, rs_w, shape_w)

        call cache_init_s(cold, N_PARAMS, N_POINTS, thetas, init_c)
        call assert_int_eq(init_c, SHAPE_VALID, label // ': cold cache init')
        call cache_radius_grid_s(cold, params, radii_c, grid_c)
        call cache_shape_s(cold, params, zs_c, rn_c, rs_c, shape_c)
        call cache_free_s(cold)

        call assert_int_eq(grid_w, grid_c, label // ': radius grid status warm == cold')
        call assert_int_eq(shape_w, shape_c, label // ': shape status warm == cold')

        do n = 1_ik, N_THETA
            call assert_bits_eq(radii_w(n), radii_c(n), &
                    label // ': radius node warm == cold, bitwise')
        end do
        call assert_bits_eq(zs_w, zs_c, label // ': z_shift warm == cold, bitwise')
        call assert_bits_eq(rn_w, rn_c, label // ': r_north warm == cold, bitwise')
        call assert_bits_eq(rs_w, rs_c, label // ': r_south warm == cold, bitwise')

    end subroutine compare_warm_cold_s

    !> A rejection in the middle of a sweep must not contaminate what follows:
    !! the engine goes cold, and the next accepted vector still matches a cold
    !! cache bit for bit.
    subroutine run_invalid_resume_s(base, rejected)

        real(kind = rk), intent(in) :: base(N_PARAMS)
        real(kind = rk), intent(in) :: rejected(N_PARAMS)

        type(cache_t) :: warm
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii(N_THETA), z_shift, r_north, r_south
        integer(kind = ik) :: j, status

        call cache_init_s(warm, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'resume: warm cache init')

        do j = 1_ik, N_PARAMS
            params(j) = base(j)
        end do

        ! Warm the cache on a valid vector.
        call compare_warm_cold_s(warm, params, 'resume, before rejection')

        ! Mid-sweep beak rejection: zero-filled, engine back to cold.
        call cache_radius_grid_s(warm, rejected, radii, status)
        call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, &
                'resume: mid-sweep vector is beak-rejected')
        call cache_shape_s(warm, rejected, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, &
                'resume: mid-sweep shape call is beak-rejected too')

        ! Resume: still bit-equal to cold.
        call compare_warm_cold_s(warm, params, 'resume, after rejection')

        params(1) = base(1) + DELTA_C
        call compare_warm_cold_s(warm, params, 'resume, stepped after rejection')

        call cache_free_s(warm)

    end subroutine run_invalid_resume_s

    !> Bit-level equality without a real `==` (banned by -Wcompare-reals).
    pure function bits_eq_f(a, b) result(same)
        real(kind = rk), intent(in) :: a, b
        logical :: same
        same = transfer(a, 0_ikl) == transfer(b, 0_ikl)
    end function bits_eq_f

end program fos_param_bitwise_test
