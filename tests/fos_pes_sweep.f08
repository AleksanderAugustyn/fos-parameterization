!> PES-range sweep gate: incremental == cold, bit for bit, over a PES box.
!!
!! The contract's bitwise rule is asserted by `fos_param_bitwise_test` on four
!! hand-chosen regimes. This gate asserts the SAME property statistically, over
!! the parameter box a PES scan actually walks, with every one of the eight
!! cached slots moving:
!!
!!   pass 1 (coarse) — c x a3 x a4 on a 15 x 8 x 8 grid, higher coefficients 0.
!!   pass 2 (fine)   — a reduced (c, a3, a4) base times a5, a6 (5 values each),
!!                     a7, a8 and a9 (3 values each). a9 is p8, the last cached
!!                     slot; without it one mask column would never be dirtied.
!!
!! Every point is computed twice: once in ONE cache walked through the whole
!! grid in odometer order (so each step reuses whatever its parameter diff did
!! not invalidate), and once in a cache created fresh for that point alone. The
!! two radius grids must be bit-identical and the two statuses must agree —
!! rejections included, since a rejected call zero-fills and the zero pattern is
!! part of the property.
!!
!! Rejections are a normal outcome of a PES box and are counted, not failed. The
!! histogram and the throughput are PRINTED for the validation record; neither
!! is asserted, so no CI threshold can go flaky. Asserted: zero bitwise
!! mismatches, zero status mismatches, no misuse status (1-6) anywhere, and a
!! nonzero valid count (a sweep that rejected everything would pass the bitwise
!! property vacuously).
program fos_pes_sweep

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, SHAPE_VALID, FOS_ERROR_RHO_NEGATIVE, &
            FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE, &
            FOS_ERROR_BUFFER_MISMATCH
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            test_summary

    implicit none

    integer(kind = ik), parameter :: N_DIMS = 8_ik
    integer(kind = ik), parameter :: N_POINTS = 201_ik
    integer(kind = ik), parameter :: N_THETA = 41_ik

    !> Histogram buckets: the seven reachable codes plus one catch-all. A count
    !! in the catch-all means a misuse status (1-6) or an unknown code, and the
    !! gate fails on it.
    integer(kind = ik), parameter :: N_CODES = 7_ik
    integer(kind = ik), parameter :: CODES(N_CODES) = &
            [SHAPE_VALID, FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
             FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY, &
             FOS_ERROR_CONVERGENCE, FOS_ERROR_BUFFER_MISMATCH]
    character(len = 26), parameter :: CODE_NAMES(N_CODES) = &
            ['  0 valid                 ', '100 rho <= 0              ', &
             '101 not star-convex       ', '102 invalid c             ', &
             '103 beak singularity      ', '104 newton not converged  ', &
             '105 buffer mismatch       ']

    !> Coarse pass: the PES box proper. Slots 4-8 pinned to 0.
    real(kind = rk), parameter :: COARSE_LO(N_DIMS) = &
            [1.00_rk, 0.00_rk, -0.09_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    real(kind = rk), parameter :: COARSE_HI(N_DIMS) = &
            [2.40_rk, 0.21_rk,  0.21_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    integer(kind = ik), parameter :: COARSE_N(N_DIMS) = &
            [15_ik, 8_ik, 8_ik, 1_ik, 1_ik, 1_ik, 1_ik, 1_ik]

    !> Fine pass: a reduced (c, a3, a4) base — the coarse pass already covers
    !! that plane — times the five higher coefficients. The half-widths are set
    !! wide enough that both shape gates fire inside the pass (measured: ~47 %
    !! beak, a handful of star-convexity rejections) while a majority stays
    !! valid: the bitwise property has to hold across gate boundaries, not only
    !! on accepted shapes.
    real(kind = rk), parameter :: FINE_LO(N_DIMS) = &
            [1.20_rk, 0.00_rk, 0.00_rk, -0.18_rk, -0.18_rk, -0.10_rk, -0.10_rk, -0.08_rk]
    real(kind = rk), parameter :: FINE_HI(N_DIMS) = &
            [2.00_rk, 0.12_rk, 0.12_rk,  0.18_rk,  0.18_rk,  0.10_rk,  0.10_rk,  0.08_rk]
    integer(kind = ik), parameter :: FINE_N(N_DIMS) = &
            [3_ik, 2_ik, 2_ik, 5_ik, 5_ik, 3_ik, 3_ik, 3_ik]

    real(kind = rk)    :: thetas(N_THETA)
    integer(kind = ik) :: i

    ! Open-uniform nodes: strictly inside (0, pi), so no endpoint can overshoot
    ! the domain bound by an ulp under -ffast-math.
    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do

    call run_pass_s('coarse (c, a3, a4)      ', COARSE_LO, COARSE_HI, COARSE_N, thetas)
    call run_pass_s('fine   (all eight slots)', FINE_LO, FINE_HI, FINE_N, thetas)

    call test_summary()

contains

    !> Sweeps one rectangular grid, comparing incremental against cold.
    !!
    !! @param[in] label   Pass name for the printed block and assertion labels
    !! @param[in] lo      Per-slot lower bound
    !! @param[in] hi      Per-slot upper bound (ignored where n(d) == 1)
    !! @param[in] n       Per-slot value count, >= 1
    !! @param[in] thetas  Polar nodes both caches are initialized with
    subroutine run_pass_s(label, lo, hi, n, thetas)

        character(len = *), intent(in) :: label
        real(kind = rk),    intent(in) :: lo(N_DIMS), hi(N_DIMS)
        integer(kind = ik), intent(in) :: n(N_DIMS)
        real(kind = rk),    intent(in) :: thetas(N_THETA)

        type(cache_t)       :: warm, cold
        real(kind = rk)     :: params(N_DIMS), step(N_DIMS)
        real(kind = rk)     :: r_warm(N_THETA), r_cold(N_THETA)
        real(kind = rk)     :: wall_seconds, rate_hz
        integer(kind = ik)  :: idx(N_DIMS), d, b, t, status_warm, status_cold
        integer(kind = ik)  :: status, n_cold_init_fail
        integer(kind = ikl) :: n_total, hist(N_CODES + 1_ik)
        integer(kind = ikl) :: tick_start, tick_end, tick_rate

        do d = 1_ik, N_DIMS
            if (n(d) > 1_ik) then
                step(d) = (hi(d) - lo(d)) / real(n(d) - 1_ik, rk)
            else
                step(d) = 0.0_rk
            end if
        end do

        n_total = 1_ikl
        do d = 1_ik, N_DIMS
            n_total = n_total * int(n(d), ikl)
        end do

        call cache_init_s(warm, N_DIMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, label // ': warm cache init')
        if (status /= SHAPE_VALID) return

        hist(:) = 0_ikl
        n_cold_init_fail = 0_ik
        idx(:) = 1_ik
        call system_clock(count = tick_start, count_rate = tick_rate)

        do
            do d = 1_ik, N_DIMS
                params(d) = lo(d) + real(idx(d) - 1_ik, rk) * step(d)
            end do

            call cache_radius_grid_s(warm, params, r_warm, status_warm)

            call cache_init_s(cold, N_DIMS, N_POINTS, thetas, status)
            if (status /= SHAPE_VALID) then
                n_cold_init_fail = n_cold_init_fail + 1_ik
            else
                call cache_radius_grid_s(cold, params, r_cold, status_cold)
                call cache_free_s(cold)

                call assert_int_eq(status_warm, status_cold, &
                        label // ': incremental status == cold status')
                do t = 1_ik, N_THETA
                    call assert_bits_eq(r_warm(t), r_cold(t), &
                            label // ': incremental radius == cold radius')
                end do
            end if

            b = code_bucket_f(status_warm)
            hist(b) = hist(b) + 1_ikl

            ! Odometer step, last slot fastest: consecutive points differ in as
            ! few parameters as possible, which is what makes the incremental
            ! cache do incremental work.
            d = N_DIMS
            do
                idx(d) = idx(d) + 1_ik
                if (idx(d) <= n(d)) exit
                idx(d) = 1_ik
                d = d - 1_ik
                if (d < 1_ik) exit
            end do
            if (d < 1_ik) exit
        end do

        call system_clock(count = tick_end)
        call cache_free_s(warm)

        wall_seconds = real(tick_end - tick_start, rk) / real(tick_rate, rk)
        if (wall_seconds > 0.0_rk) then
            rate_hz = real(n_total, rk) / wall_seconds
        else
            rate_hz = 0.0_rk
        end if

        write(*, '(A)') repeat('-', 68)
        write(*, '(A,A)')  'PES sweep pass: ', label
        write(*, '(A)') repeat('-', 68)
        write(*, '(A,8(1X,I0))') '  per-slot counts    :', n(:)
        write(*, '(A,I0)')       '  grid points        : ', n_total
        do b = 1_ik, N_CODES
            write(*, '(A,A,A,I0,A,F7.3,A)') '  ', CODE_NAMES(b), ': ', hist(b), &
                    '  (', 100.0_rk * real(hist(b), rk) / real(n_total, rk), ' %)'
        end do
        write(*, '(A,I0)')    '  unexpected statuses : ', hist(N_CODES + 1_ik)
        write(*, '(A,F12.3)') '  wall seconds        : ', wall_seconds
        write(*, '(A,F14.1)') '  shapes / second     : ', rate_hz

        ! Asserted outcomes. The histogram above is a record, not a gate; these
        ! four are the gate.
        call assert_int_eq(n_cold_init_fail, 0_ik, label // ': cold cache inits all succeeded')
        call assert_true(hist(N_CODES + 1_ik) == 0_ikl, &
                label // ': no misuse or unknown status')
        call assert_true(hist(1) > 0_ikl, label // ': at least one valid shape')

    end subroutine run_pass_s

    !> Histogram bucket for a status code; N_CODES + 1 for anything unlisted.
    !!
    !! @param[in] code  Status returned by a compute call
    !! @return          Bucket index in 1 .. N_CODES + 1
    pure function code_bucket_f(code) result(bucket)

        integer(kind = ik), intent(in) :: code
        integer(kind = ik) :: bucket

        integer(kind = ik) :: j

        bucket = N_CODES + 1_ik
        do j = 1_ik, N_CODES
            if (code == CODES(j)) then
                bucket = j
                exit
            end if
        end do

    end function code_bucket_f

end program fos_pes_sweep
