!> PES-range sweep gate: shared cache == fresh cache, bit for bit, over a PES box.
!!
!! The contract's equivalence and statelessness rules are asserted by
!! `fos_param_equivalence_test` and `fos_param_statelessness_test` on
!! hand-chosen regimes. This gate asserts the SAME property statistically, over
!! the parameter box a PES scan actually walks, with every one of the eight
!! parameter slots moving:
!!
!!   pass 1 (coarse) — c x a3 x a4 on a 15 x 8 x 8 grid, higher coefficients 0.
!!   pass 2 (fine)   — a reduced (c, a3, a4) base times a5, a6 (5 values each),
!!                     a7, a8 and a9 (3 values each).
!!
!! Every point is computed twice: once on ONE shared cache walked through the
!! whole grid in odometer order — the way a PES scan uses the library — and
!! once on a cache built for that point alone, in a different procedure
!! (`fresh_point_s`), so the comparison also crosses call sites. The two radius
!! grids must be bit-identical and the two statuses must agree — rejections
!! included, since a rejected call zero-fills and the zero pattern is part of
!! the property.
!!
!! Rejections are a normal outcome of a PES box and are counted, not failed. The
!! histogram and the two throughputs are PRINTED for the validation record; none
!! of them is asserted, so no CI threshold can go flaky. "shared" times the
!! compute alone on the shared cache — the rate a PES scan sees, and the number
!! the no-LTO cost is measured on; "fresh" times a cache's whole life (table
!! build + compute + free, dominated by the table build). Each counts one shape
!! per grid point. Asserted: zero bitwise mismatches, zero status mismatches, no
!! usage status (1-5, 105) anywhere, and a nonzero valid count (a sweep that
!! rejected everything would pass the bitwise property vacuously).
program fos_pes_sweep

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_grid_s, SHAPE_VALID, FOS_ERROR_RHO_NEGATIVE, &
            FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            test_summary

    implicit none

    integer(kind = ik), parameter :: N_DIMS = 8_ik
    integer(kind = ik), parameter :: N_POINTS = 201_ik
    integer(kind = ik), parameter :: N_THETA = 41_ik

    !> Histogram buckets: success and the five value codes, plus one catch-all.
    !! A count in the catch-all means a usage status (1-5, or 105 — the sweep
    !! passes correctly sized buffers) or an unknown code, and the gate fails
    !! on it.
    integer(kind = ik), parameter :: N_CODES = 6_ik
    integer(kind = ik), parameter :: CODES(N_CODES) = &
            [SHAPE_VALID, FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
             FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY, &
             FOS_ERROR_CONVERGENCE]
    character(len = 26), parameter :: CODE_NAMES(N_CODES) = &
            ['  0 valid                 ', '100 rho <= 0              ', &
             '101 not star-convex       ', '102 invalid c             ', &
             '103 beak singularity      ', '104 newton not converged  ']

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

    !> Sweeps one rectangular grid, comparing the shared cache against a fresh
    !! cache per point.
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

        type(cache_t)       :: shared
        real(kind = rk)     :: params(N_DIMS), step(N_DIMS)
        real(kind = rk)     :: r_shared(N_THETA), r_fresh(N_THETA)
        real(kind = rk)     :: shared_seconds, fresh_seconds
        integer(kind = ik)  :: idx(N_DIMS), d, b, t, status_shared, status_fresh
        integer(kind = ik)  :: status, n_fresh_init_fail
        integer(kind = ikl) :: n_total, hist(N_CODES + 1_ik)
        integer(kind = ikl) :: shared_ticks, fresh_ticks, t0, t1, tick_rate
        logical             :: fresh_ok

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

        call cache_init_s(shared, N_DIMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, label // ': shared cache init')
        if (status /= SHAPE_VALID) return

        hist(:) = 0_ikl
        n_fresh_init_fail = 0_ik
        shared_ticks = 0_ikl
        fresh_ticks = 0_ikl
        idx(:) = 1_ik

        call system_clock(count_rate = tick_rate)

        do
            do d = 1_ik, N_DIMS
                params(d) = lo(d) + real(idx(d) - 1_ik, rk) * step(d)
            end do

            ! Timed segment 1: one compute on the shared cache. This is the
            ! rate a PES scan actually sees.
            call system_clock(count = t0)
            call cache_radius_grid_s(shared, params, r_shared, status_shared)
            call system_clock(count = t1)
            shared_ticks = shared_ticks + (t1 - t0)

            ! Timed segment 2: the fresh-cache reference — table build, one
            ! compute and the free, in another procedure.
            call system_clock(count = t0)
            call fresh_point_s(params, thetas, r_fresh, status_fresh, fresh_ok)
            call system_clock(count = t1)
            fresh_ticks = fresh_ticks + (t1 - t0)

            if (.not. fresh_ok) then
                n_fresh_init_fail = n_fresh_init_fail + 1_ik
            else
                ! Comparison is outside both timed segments: the assertions are
                ! the gate's work, not the library's.
                call assert_int_eq(status_shared, status_fresh, &
                        label // ': shared status == fresh status')
                do t = 1_ik, N_THETA
                    call assert_bits_eq(r_shared(t), r_fresh(t), &
                            label // ': shared radius == fresh radius')
                end do
            end if

            b = code_bucket_f(status_shared)
            hist(b) = hist(b) + 1_ikl

            ! Odometer step, last slot fastest.
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

        call cache_free_s(shared)

        shared_seconds = real(shared_ticks, rk) / real(tick_rate, rk)
        fresh_seconds = real(fresh_ticks, rk) / real(tick_rate, rk)

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
        ! Both rates count ONE shape per grid point, over the segment named.
        write(*, '(A,F12.3)') '  shared seconds      : ', shared_seconds
        write(*, '(A,F14.1)') '  shared shapes/s     : ', rate_f(n_total, shared_seconds)
        write(*, '(A,F12.3)') '  fresh seconds       : ', fresh_seconds
        write(*, '(A,F14.1)') '  fresh shapes/s (incl. cache init + free): ', &
                rate_f(n_total, fresh_seconds)

        ! Asserted outcomes. The histogram above is a record, not a gate; these
        ! are the gate.
        call assert_int_eq(n_fresh_init_fail, 0_ik, label // ': fresh cache inits all succeeded')
        call assert_true(hist(N_CODES + 1_ik) == 0_ikl, &
                label // ': no usage or unknown status')
        call assert_true(hist(1) > 0_ikl, label // ': at least one valid shape')

    end subroutine run_pass_s

    !> One point on a cache built for it alone. A separate procedure from the
    !! shared-cache call on purpose: the comparison must cross call sites.
    !!
    !! @param[in]  params  Parameter vector of the grid point
    !! @param[in]  thetas  Polar nodes
    !! @param[out] radii   R(theta); zero-filled on any rejection
    !! @param[out] status  Status of the compute (undefined when init failed)
    !! @param[out] ok      .false. iff the fresh cache could not be built
    subroutine fresh_point_s(params, thetas, radii, status, ok)

        real(kind = rk),    intent(in)  :: params(N_DIMS)
        real(kind = rk),    intent(in)  :: thetas(N_THETA)
        real(kind = rk),    intent(out) :: radii(N_THETA)
        integer(kind = ik), intent(out) :: status
        logical,            intent(out) :: ok

        type(cache_t) :: fresh
        integer(kind = ik) :: init_status

        radii = 0.0_rk
        status = SHAPE_VALID

        call cache_init_s(fresh, N_DIMS, N_POINTS, thetas, init_status)
        ok = init_status == SHAPE_VALID
        if (.not. ok) return

        call cache_radius_grid_s(fresh, params, radii, status)
        call cache_free_s(fresh)

    end subroutine fresh_point_s

    !> Shapes per second, 0 when the segment was too short for the clock.
    !!
    !! @param[in] n_shapes  Shapes resolved in the segment
    !! @param[in] seconds   Segment duration
    !! @return              Rate in Hz
    pure function rate_f(n_shapes, seconds) result(rate_hz)

        integer(kind = ikl), intent(in) :: n_shapes
        real(kind = rk),     intent(in) :: seconds
        real(kind = rk) :: rate_hz

        if (seconds > 0.0_rk) then
            rate_hz = real(n_shapes, rk) / seconds
        else
            rate_hz = 0.0_rk
        end if

    end function rate_f

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
