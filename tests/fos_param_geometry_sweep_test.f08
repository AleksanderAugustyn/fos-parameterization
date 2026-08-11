!> Geometry sweep over the FoS parameter box, ported from WMMM
!! tests/fos_geometry_validation_test.f08 and adapted to the 2.0 tier-1 API.
!!
!! Tier A classifies every grid point through the production gates
!! (`compute_shape_standalone_s`) and pins the verdict counts against golden
!! tallies captured on the default coarse grid. Tier B (accuracy) and the
!! representability probe are added by later tasks.
!!
!! Grid (the VALIDATED box): c in [0.75, 3.50], a3 in [0.00, 0.60],
!! a4 in [-0.20, 0.75], a5 in [-0.15, 0.15], a6 in [-0.10, 0.10]; default
!! step 0.05 on every axis (509,600 points), `--full` refines a5/a6 to 0.01
!! (~9.5M points). The a6 bound is part of the 2026-08-11 F_MIN_THRESHOLD
!! statement: at |a6| = 0.15 an a4 = 0.75, a3 = 0.55, c <= 0.85 family with
!! healthy f_min degrades GL-4096 V/S (see the constant's doc comment);
!! --a5-abs / --a6-abs re-run wider boxes as experiments.
module fos_sweep_support_mod

    use precision_utilities_mod, only: ik, rk

    implicit none

    private

    public :: sweep_config_t, parse_args_s, print_config_s, grid_params_f
    public :: C_LO, C_HI, A3_LO, A3_HI, A4_LO, A4_HI, A5_LO, A5_HI, A6_LO, A6_HI
    public :: N_POINTS_SWEEP

    real(kind = rk), parameter :: C_LO = 0.75_rk, C_HI = 3.50_rk
    real(kind = rk), parameter :: A3_LO = 0.00_rk, A3_HI = 0.60_rk
    real(kind = rk), parameter :: A4_LO = -0.20_rk, A4_HI = 0.75_rk
    real(kind = rk), parameter :: A5_LO = -0.15_rk, A5_HI = 0.15_rk
    real(kind = rk), parameter :: A6_LO = -0.10_rk, A6_HI = 0.10_rk

    !> u-grid resolution of every sweep-tier conversion (wmmm N_RHO_INTERNAL
    !! parity).
    integer(kind = ik), parameter :: N_POINTS_SWEEP = 1001_ik

    type :: sweep_config_t
        real(kind = rk) :: dc = 0.05_rk, da3 = 0.05_rk, da4 = 0.05_rk
        real(kind = rk) :: da5 = 0.05_rk, da6 = 0.05_rk
        !> Axis bounds; runtime-overridable (--a5-abs / --a6-abs) for
        !! box-trimming experiments.
        real(kind = rk) :: a5_lo = A5_LO, a5_hi = A5_HI
        real(kind = rk) :: a6_lo = A6_LO, a6_hi = A6_HI
        integer(kind = ik) :: n_c = 0_ik, n_a3 = 0_ik, n_a4 = 0_ik
        integer(kind = ik) :: n_a5 = 0_ik, n_a6 = 0_ik
        logical :: probe = .false., skip_tier_b = .false.
        !> .true. iff no step/range override and no --full: golden tallies apply.
        logical :: golden_grid = .true.
    end type sweep_config_t

contains

    pure function count_steps_f(lo, hi, step) result(n)

        real(kind = rk), intent(in) :: lo, hi, step
        integer(kind = ik) :: n

        n = floor((hi - lo) / step + 1.0e-6_rk, kind = ik) + 1_ik

    end function count_steps_f

    subroutine parse_args_s(cfg)

        type(sweep_config_t), intent(inout) :: cfg

        integer(kind = ik) :: i, n_args
        character(len = 512) :: arg

        n_args = int(command_argument_count(), ik)
        i = 1_ik
        do while (i <= n_args)
            call get_command_argument(int(i), arg)
            select case (trim(arg))
            case ('--full')
                cfg%dc = 0.05_rk
                cfg%da3 = 0.05_rk
                cfg%da4 = 0.05_rk
                cfg%da5 = 0.01_rk
                cfg%da6 = 0.01_rk
                cfg%golden_grid = .false.
            case ('--probe')
                cfg%probe = .true.
            case ('--skip-tier-b')
                cfg%skip_tier_b = .true.
            case ('--dc')
                call read_real_arg_s(i, cfg%dc)
                cfg%golden_grid = .false.
            case ('--da3')
                call read_real_arg_s(i, cfg%da3)
                cfg%golden_grid = .false.
            case ('--da4')
                call read_real_arg_s(i, cfg%da4)
                cfg%golden_grid = .false.
            case ('--da5')
                call read_real_arg_s(i, cfg%da5)
                cfg%golden_grid = .false.
            case ('--da6')
                call read_real_arg_s(i, cfg%da6)
                cfg%golden_grid = .false.
            case ('--a5-abs')
                call read_real_arg_s(i, cfg%a5_hi)
                cfg%a5_lo = -cfg%a5_hi
                cfg%golden_grid = .false.
            case ('--a6-abs')
                call read_real_arg_s(i, cfg%a6_hi)
                cfg%a6_lo = -cfg%a6_hi
                cfg%golden_grid = .false.
            case default
                write(*, '(A,A)') 'Unknown argument: ', trim(arg)
                write(*, '(A)') 'Usage: fos_param_geometry_sweep_test [--full]' // &
                        ' [--probe] [--skip-tier-b]' // &
                        ' [--dc X] [--da3 X] [--da4 X] [--da5 X] [--da6 X]' // &
                        ' [--a5-abs X] [--a6-abs X]'
                stop 2
            end select
            i = i + 1_ik
        end do

        cfg%n_c = count_steps_f(C_LO, C_HI, cfg%dc)
        cfg%n_a3 = count_steps_f(A3_LO, A3_HI, cfg%da3)
        cfg%n_a4 = count_steps_f(A4_LO, A4_HI, cfg%da4)
        cfg%n_a5 = count_steps_f(cfg%a5_lo, cfg%a5_hi, cfg%da5)
        cfg%n_a6 = count_steps_f(cfg%a6_lo, cfg%a6_hi, cfg%da6)

    contains

        subroutine read_real_arg_s(pos, val)

            integer(kind = ik), intent(inout) :: pos
            real(kind = rk), intent(out) :: val

            character(len = 512) :: buf
            integer :: ios

            pos = pos + 1_ik
            if (pos > n_args) then
                write(*, '(A)') 'Missing value after step-size flag'
                stop 2
            end if
            call get_command_argument(int(pos), buf)
            read(buf, *, iostat = ios) val
            if (ios /= 0 .or. val <= 0.0_rk) then
                write(*, '(A,A)') 'Invalid step value: ', trim(buf)
                stop 2
            end if

        end subroutine read_real_arg_s

    end subroutine parse_args_s

    subroutine print_config_s(cfg)

        type(sweep_config_t), intent(in) :: cfg

        write(*, '(A)') 'Sweep configuration:'
        write(*, '(A,F6.3,A,I0,A)') '  dc  = ', cfg%dc, '  (', cfg%n_c, ' values)'
        write(*, '(A,F6.3,A,I0,A)') '  da3 = ', cfg%da3, '  (', cfg%n_a3, ' values)'
        write(*, '(A,F6.3,A,I0,A)') '  da4 = ', cfg%da4, '  (', cfg%n_a4, ' values)'
        write(*, '(A,F6.3,A,I0,A)') '  da5 = ', cfg%da5, '  (', cfg%n_a5, ' values)'
        write(*, '(A,F6.3,A,I0,A)') '  da6 = ', cfg%da6, '  (', cfg%n_a6, ' values)'
        write(*, '(A,F6.3,A,F6.3,A)') '  a5 in [', cfg%a5_lo, ', ', cfg%a5_hi, ']'
        write(*, '(A,F6.3,A,F6.3,A)') '  a6 in [', cfg%a6_lo, ', ', cfg%a6_hi, ']'
        write(*, '(A,I0)') '  grid points: ', &
                cfg%n_c * cfg%n_a3 * cfg%n_a4 * cfg%n_a5 * cfg%n_a6
        write(*, '(A,L1)') '  golden grid: ', cfg%golden_grid

    end subroutine print_config_s

    !> Grid point (i_c, i3, i4, i5, i6) -> 7-slot parameter vector
    !! [c, a3, a4, a5, a6, 0, 0] (a2 implicit via the volume constraint).
    pure function grid_params_f(cfg, i_c, i3, i4, i5, i6) result(params)

        type(sweep_config_t), intent(in) :: cfg
        integer(kind = ik), intent(in) :: i_c, i3, i4, i5, i6
        real(kind = rk) :: params(7)

        params = [C_LO + real(i_c - 1_ik, rk) * cfg%dc, &
                A3_LO + real(i3 - 1_ik, rk) * cfg%da3, &
                A4_LO + real(i4 - 1_ik, rk) * cfg%da4, &
                cfg%a5_lo + real(i5 - 1_ik, rk) * cfg%da5, &
                cfg%a6_lo + real(i6 - 1_ik, rk) * cfg%da6, &
                0.0_rk, 0.0_rk]

    end function grid_params_f

end module fos_sweep_support_mod

program fos_param_geometry_sweep_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_shape_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            compute_f_min_standalone_s, compute_conversion_diagnostic_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            F_MIN_THRESHOLD, STAR_CONVEXITY_MARGIN, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
            FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE
    use shape_core_mod, only: SHAPE_VALID
    use fos_sweep_support_mod, only: sweep_config_t, parse_args_s, &
            print_config_s, grid_params_f, N_POINTS_SWEEP
    use fos_test_reference_mod, only: init_quadrature_s, &
            compute_reference_surface_f, evaluate_shape_quality_s, &
            evaluate_shape_quality_diag_s, find_neck_in_radii_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            test_summary

    implicit none

    !> Tier B tolerances (absolute), re-baselined 2026-08-11 on the FINE grid
    !! of the validated box (a6 in [-0.10, 0.10], a5/a6 step 0.01, 9.5M pts)
    !! with F_MIN_THRESHOLD = 5e-4, against the GL-4096 measurand.
    !!
    !! Fine-probe worst over the accepted f_min bins above 5e-4: dV/V
    !! 8.5e-8, dS/S 4.3e-5, round-trip 6.3e-12 — asserted roughly one decade
    !! above. The V/S envelope is set by high-f_min family tails
    !! (f_min ~ 2e-3..4e-3, high c with a3 at the box edge, resolved origin
    !! near the south tip; conversion exact, GL-4096 quadrature-limited —
    !! dV/V converges ~N^-7 while round-trip stays ~5e-12) that NO legal
    !! threshold can evict; the threshold-adjacent quadrature band
    !! ([3.162e-4, 4.642e-4), fine dS 2.4e-4) sits below 5e-4 and is
    !! rejected. On the pre-trim box (|a6| = 0.15) the a4 = 0.75 corner
    !! family broke dS to 6.2e-3 at healthy f_min — outside the validated
    !! box. Round-trip is the conversion-correctness pin.
    real(kind = rk), parameter :: TOL_VOLUME_B = 1.0e-6_rk
    real(kind = rk), parameter :: TOL_SURFACE_REL = 5.0e-4_rk
    real(kind = rk), parameter :: TOL_ROUND_TRIP = 1.0e-10_rk

    !> Neck-context thresholds, mirroring wmmm's Level 2c physics filter
    !! (NECK_MIN_DEPTH, NECK_MIN_ELONGATION). Context for the printed report
    !! ONLY — deliberately NOT a gate on the worst-case statistics: fos-param's
    !! accuracy guarantee covers every shape it accepts.
    real(kind = rk), parameter :: NECK_DEPTH_MIRROR = 0.25_rk
    real(kind = rk), parameter :: NECK_ELONGATION_MIRROR = 1.2_rk

    !> Tier A verdict counts on the default coarse grid of the validated box
    !! (a6 in [-0.10, 0.10], F_MIN_THRESHOLD = 5e-4), captured 2026-08-11.
    !! Index: 0 = valid, 1..4 = codes 100..103, 5 = other.
    integer(kind = ik), parameter :: GOLDEN_TALLY(0:5) = &
            [390409_ik, 0_ik, 1087_ik, 0_ik, 118104_ik, 0_ik]

    type(sweep_config_t) :: cfg
    integer(kind = ik) :: tally(0:5), n_total

    call parse_args_s(cfg)
    call print_config_s(cfg)
    n_total = cfg%n_c * cfg%n_a3 * cfg%n_a4 * cfg%n_a5 * cfg%n_a6

    call run_tier_a_s(cfg, tally)

    call assert_int_eq(sum(tally), n_total, 'tier A: every point classified')
    call assert_int_eq(tally(5), 0_ik, 'tier A: no unexpected codes')
    if (cfg%golden_grid) then
        call assert_int_eq(tally(0), GOLDEN_TALLY(0), 'tier A golden: valid')
        call assert_int_eq(tally(1), GOLDEN_TALLY(1), 'tier A golden: rho (100)')
        call assert_int_eq(tally(2), GOLDEN_TALLY(2), 'tier A golden: star (101)')
        call assert_int_eq(tally(3), GOLDEN_TALLY(3), 'tier A golden: c (102)')
        call assert_int_eq(tally(4), GOLDEN_TALLY(4), 'tier A golden: beak (103)')
        call assert_int_eq(tally(5), GOLDEN_TALLY(5), 'tier A golden: other')
    else
        write(*, '(A)') 'Non-default grid: golden tally asserts skipped.'
    end if

    if (.not. cfg%skip_tier_b .or. cfg%probe) call init_quadrature_s()
    if (.not. cfg%skip_tier_b) call run_tier_b_s(cfg)
    if (cfg%probe) call run_probe_s(cfg)

    call test_summary()

contains

    !> Tier A: classify every grid point through the production standalone
    !! gates and tally the verdicts.
    subroutine run_tier_a_s(cfg, tally_out)

        type(sweep_config_t), intent(in) :: cfg
        integer(kind = ik), intent(out) :: tally_out(0:5)

        integer(kind = ik) :: i_c, i3, i4, i5, i6, status
        real(kind = rk) :: params(7), z_shift, r_north, r_south
        integer(kind = ik) :: tally(0:5)
        integer(kind = ik) :: t_start, t_end, t_rate

        write(*, '(A)') '=== Tier A: gate classification sweep ==='

        tally = 0_ik
        call system_clock(t_start, t_rate)

        !$omp parallel do collapse(3) schedule(dynamic) default(shared) &
        !$omp& private(i_c, i3, i4, i5, i6, params, z_shift, r_north, r_south, status) &
        !$omp& reduction(+:tally)
        do i_c = 1_ik, cfg%n_c
            do i3 = 1_ik, cfg%n_a3
                do i4 = 1_ik, cfg%n_a4
                    do i5 = 1_ik, cfg%n_a5
                        do i6 = 1_ik, cfg%n_a6
                            params = grid_params_f(cfg, i_c, i3, i4, i5, i6)
                            call compute_shape_standalone_s(params, N_POINTS_SWEEP, &
                                    z_shift, r_north, r_south, status)
                            select case (status)
                            case (SHAPE_VALID);                tally(0) = tally(0) + 1_ik
                            case (FOS_ERROR_RHO_NEGATIVE);     tally(1) = tally(1) + 1_ik
                            case (FOS_ERROR_NOT_STAR_CONVEX);  tally(2) = tally(2) + 1_ik
                            case (FOS_ERROR_INVALID_C);        tally(3) = tally(3) + 1_ik
                            case (FOS_ERROR_BEAK_SINGULARITY); tally(4) = tally(4) + 1_ik
                            case default;                      tally(5) = tally(5) + 1_ik
                            end select
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        call system_clock(t_end)

        write(*, '(A,I0)') '  valid (0):            ', tally(0)
        write(*, '(A,I0)') '  rho negative (100):   ', tally(1)
        write(*, '(A,I0)') '  not star-convex (101): ', tally(2)
        write(*, '(A,I0)') '  invalid c (102):      ', tally(3)
        write(*, '(A,I0)') '  beak f_min (103):     ', tally(4)
        write(*, '(A,I0)') '  other:                ', tally(5)
        write(*, '(A,F10.2,A)') '  Tier A wall time: ', &
                real(t_end - t_start, rk) / real(t_rate, rk), ' s'

        tally_out = tally

    end subroutine run_tier_a_s

    !> Tier B: conversion accuracy on every production-accepted shape.
    !!
    !! Per accepted point: cylindrical reference surface, then the GL-4096
    !! V/S/round-trip metrics through the production conversion. Worst cases
    !! are tracked over ALL accepted shapes — no physics filter (deliberate
    !! departure from wmmm): fos-param's guarantee covers every shape it
    !! accepts, and filtering would hide a conversion defect exactly where
    !! wmmm's replacement energy test won't look.
    subroutine run_tier_b_s(cfg)

        type(sweep_config_t), intent(in) :: cfg

        integer(kind = ik), parameter :: N_NECK_GRID = 101_ik

        integer(kind = ik) :: i_c, i3, i4, i5, i6, status, q_status, n_status, i
        real(kind = rk) :: params(7), z_shift, r_north, r_south
        real(kind = rk) :: s_ref, dv_rel, ds_rel, rt_max
        real(kind = rk) :: worst_dv, worst_ds, worst_rt
        real(kind = rk) :: worst_dv_params(5), worst_ds_params(5), worst_rt_params(5)
        integer(kind = ik) :: n_valid, n_neck, n_eval_fail
        real(kind = rk) :: thetas_neck(N_NECK_GRID)
        real(kind = rk) :: radii_neck(N_NECK_GRID), dr_neck(N_NECK_GRID)
        logical :: has_neck
        real(kind = rk) :: neck_radius, neck_depth
        integer(kind = ik) :: t_start, t_end, t_rate

        write(*, '(A)') '=== Tier B: conversion accuracy sweep ==='

        do i = 1_ik, N_NECK_GRID
            thetas_neck(i) = real(i - 1_ik, rk) * PI_C / real(N_NECK_GRID - 1_ik, rk)
        end do
        thetas_neck(N_NECK_GRID) = PI_C

        worst_dv = 0.0_rk
        worst_ds = 0.0_rk
        worst_rt = 0.0_rk
        worst_dv_params = 0.0_rk
        worst_ds_params = 0.0_rk
        worst_rt_params = 0.0_rk
        n_valid = 0_ik
        n_neck = 0_ik
        n_eval_fail = 0_ik

        call system_clock(t_start, t_rate)

        ! The unsynchronized reads of worst_* in the trigger condition are a
        ! benign race: worst_* grows monotonically and a stale (smaller) value
        ! only causes a redundant critical-section entry, never a missed update.
        !$omp parallel do collapse(3) schedule(dynamic) default(shared) &
        !$omp& private(i_c, i3, i4, i5, i6, params, z_shift, r_north, r_south, &
        !$omp&         status, q_status, n_status, s_ref, dv_rel, ds_rel, rt_max, &
        !$omp&         radii_neck, dr_neck, has_neck, neck_radius, neck_depth) &
        !$omp& reduction(+:n_valid, n_neck, n_eval_fail)
        do i_c = 1_ik, cfg%n_c
            do i3 = 1_ik, cfg%n_a3
                do i4 = 1_ik, cfg%n_a4
                    do i5 = 1_ik, cfg%n_a5
                        do i6 = 1_ik, cfg%n_a6
                            params = grid_params_f(cfg, i_c, i3, i4, i5, i6)
                            call compute_shape_standalone_s(params, N_POINTS_SWEEP, &
                                    z_shift, r_north, r_south, status)
                            if (status /= SHAPE_VALID) cycle
                            n_valid = n_valid + 1_ik

                            s_ref = compute_reference_surface_f(params)
                            call evaluate_shape_quality_s(params, N_POINTS_SWEEP, &
                                    z_shift, s_ref, dv_rel, ds_rel, rt_max, q_status)
                            ! A production-accepted shape that fails to evaluate
                            ! returns huge() metrics: counted, and the worst
                            ! trackers trip the tolerance asserts loudly.
                            if (q_status /= SHAPE_VALID) n_eval_fail = n_eval_fail + 1_ik

                            call compute_radius_and_derivative_standalone_s(params, &
                                    thetas_neck, N_POINTS_SWEEP, radii_neck, dr_neck, &
                                    n_status)
                            if (n_status == SHAPE_VALID) then
                                call find_neck_in_radii_s(radii_neck, has_neck, &
                                        neck_radius, neck_depth)
                                if (has_neck .and. neck_depth > NECK_DEPTH_MIRROR &
                                        .and. 0.5_rk * (radii_neck(1) &
                                        + radii_neck(N_NECK_GRID)) &
                                        < NECK_ELONGATION_MIRROR) then
                                    n_neck = n_neck + 1_ik
                                end if
                            end if

                            if (abs(dv_rel) > worst_dv .or. abs(ds_rel) > worst_ds &
                                    .or. rt_max > worst_rt) then
                                !$omp critical (tier_b_worst)
                                if (abs(dv_rel) > worst_dv) then
                                    worst_dv = abs(dv_rel)
                                    worst_dv_params = params(1:5)
                                end if
                                if (abs(ds_rel) > worst_ds) then
                                    worst_ds = abs(ds_rel)
                                    worst_ds_params = params(1:5)
                                end if
                                if (rt_max > worst_rt) then
                                    worst_rt = rt_max
                                    worst_rt_params = params(1:5)
                                end if
                                !$omp end critical (tier_b_worst)
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        call system_clock(t_end)

        write(*, '(A,I0)') '  accepted shapes evaluated: ', n_valid
        write(*, '(A,I0)') '  pronounced neck at low elongation (context): ', n_neck
        write(*, '(A,I0)') '  evaluation failures on accepted shapes: ', n_eval_fail
        write(*, '(A,ES10.3,A,5F8.3)') '  max |dV/V|:     ', worst_dv, &
                '  at c,a3,a4,a5,a6 =', worst_dv_params
        write(*, '(A,ES10.3,A,5F8.3)') '  max |dS/S|:     ', worst_ds, &
                '  at c,a3,a4,a5,a6 =', worst_ds_params
        write(*, '(A,ES10.3,A,5F8.3)') '  max round-trip: ', worst_rt, &
                '  at c,a3,a4,a5,a6 =', worst_rt_params
        write(*, '(A,F10.2,A)') '  Tier B wall time: ', &
                real(t_end - t_start, rk) / real(t_rate, rk), ' s'

        call assert_true(n_valid > 0_ik, 'tier B: at least one accepted shape')
        call assert_int_eq(n_eval_fail, 0_ik, 'tier B: every accepted shape evaluates')
        call assert_true(worst_dv <= TOL_VOLUME_B, 'tier B: |dV/V| within tolerance')
        call assert_true(worst_ds <= TOL_SURFACE_REL, 'tier B: |dS/S| within tolerance')
        call assert_true(worst_rt <= TOL_ROUND_TRIP, 'tier B: round-trip within tolerance')

    end subroutine run_tier_b_s

    !> Representability probe (--probe): July-parity g(s*) bins plus f_min
    !! bins split by branch (polar = boundary-clamped scan minimum, neck =
    !! interior), all measured through the beak- and star-ungated diagnostic
    !! conversion. REPORTS only — the cliff readout asserts nothing; the
    !! threshold decision is the Task 7 user checkpoint.
    subroutine run_probe_s(cfg)

        type(sweep_config_t), intent(in) :: cfg

        ! July-parity g bins: [-0.15, 0.02) in 0.01 steps; shapes with
        ! g(s*) < 0.01 are included (thin over-limit band), g < -0.15 clamps
        ! into bin 1.
        integer(kind = ik), parameter :: N_G_BINS = 17_ik
        real(kind = rk), parameter :: G_BIN_LO = -0.15_rk
        real(kind = rk), parameter :: G_BIN_W = 0.01_rk
        !> July-parity resolve resolution for the g values.
        integer(kind = ik), parameter :: N_RHO_JULY = 7201_ik
        ! f_min bins: 6/decade over [1e-6, 1e-1] (bins 1..30) + underflow (0).
        integer(kind = ik), parameter :: N_F_BINS = 30_ik

        integer(kind = ik) :: g_count(N_G_BINS), g_fail(N_G_BINS)
        real(kind = rk) :: g_dv(N_G_BINS), g_ds(N_G_BINS), g_rt(N_G_BINS)
        real(kind = rk) :: g_rp(N_G_BINS), g_rpp(N_G_BINS)
        integer(kind = ik) :: f_count(0:N_F_BINS, 2), f_fail(0:N_F_BINS, 2)
        real(kind = rk) :: f_dv(0:N_F_BINS, 2), f_ds(0:N_F_BINS, 2)
        real(kind = rk) :: f_rt(0:N_F_BINS, 2), f_rp(0:N_F_BINS, 2)
        real(kind = rk) :: f_rpp(0:N_F_BINS, 2)

        integer(kind = ik) :: i_c, i3, i4, i5, i6, stf, stg, stq, gb, fb, br, b
        real(kind = rk) :: params(7), f_min, u_min, g7, zt7
        real(kind = rk) :: s_ref, dv_rel, ds_rel, rt_max, rp_max, rpp_max
        real(kind = rk) :: zs_d, g_d, g_spot, zt_spot
        real(kind = rk) :: th1(1), r1(1), dr1(1)
        logical :: interior_min, beak_pass
        integer(kind = ik) :: t_start, t_end, t_rate

        write(*, '(A)') '=== Representability probe ==='

        ! Spot check: for a beak-passing shape the July-parity star optimum
        ! and the diagnostic resolve must report the SAME g(s*) bitwise (both
        ! run the identical resolve at N_RHO_JULY).
        params = [1.5_rk, 0.1_rk, 0.1_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
        call compute_star_convexity_optimum_standalone_s(params, N_RHO_JULY, &
                zt_spot, g_spot, stg)
        call assert_int_eq(stg, SHAPE_VALID, 'probe: spot-check star optimum valid')
        th1(1) = 0.5_rk * PI_C
        call compute_conversion_diagnostic_standalone_s(params, th1, N_RHO_JULY, &
                r1, dr1, zs_d, g_d, stq)
        call assert_int_eq(stq, SHAPE_VALID, 'probe: spot-check diagnostic valid')
        call assert_bits_eq(g_spot, g_d, 'probe: g(s*) bitwise across both paths')

        g_count = 0_ik
        g_fail = 0_ik
        g_dv = 0.0_rk
        g_ds = 0.0_rk
        g_rt = 0.0_rk
        g_rp = 0.0_rk
        g_rpp = 0.0_rk
        f_count = 0_ik
        f_fail = 0_ik
        f_dv = 0.0_rk
        f_ds = 0.0_rk
        f_rt = 0.0_rk
        f_rp = 0.0_rk
        f_rpp = 0.0_rk

        call system_clock(t_start, t_rate)

        !$omp parallel do collapse(3) schedule(dynamic) default(shared) &
        !$omp& private(i_c, i3, i4, i5, i6, params, f_min, u_min, interior_min, &
        !$omp&         stf, stg, stq, gb, fb, br, g7, zt7, s_ref, dv_rel, ds_rel, &
        !$omp&         rt_max, rp_max, rpp_max, zs_d, g_d, th1, r1, dr1, beak_pass) &
        !$omp& reduction(+:g_count, g_fail, f_count, f_fail) &
        !$omp& reduction(max:g_dv, g_ds, g_rt, g_rp, g_rpp, f_dv, f_ds, f_rt, f_rp, f_rpp)
        do i_c = 1_ik, cfg%n_c
            do i3 = 1_ik, cfg%n_a3
                do i4 = 1_ik, cfg%n_a4
                    do i5 = 1_ik, cfg%n_a5
                        do i6 = 1_ik, cfg%n_a6
                            params = grid_params_f(cfg, i_c, i3, i4, i5, i6)

                            call compute_f_min_standalone_s(params, f_min, u_min, &
                                    interior_min, stf)
                            if (stf /= SHAPE_VALID) cycle
                            beak_pass = f_min > F_MIN_THRESHOLD

                            ! g(s*) at July parity. The star-optimum form is
                            ! beak-gated, so it only serves beak-passing
                            ! shapes; the rest resolve through the diagnostic
                            ! path (g is set before the theta solve, so a 104
                            ! there still reports a valid g).
                            !
                            ! Known sliver: bin membership uses g at 7201
                            ! while the measured conversion resolves at
                            ! N_POINTS_SWEEP, so shapes within the two
                            ! resolutions' g difference (~1e-5) of the -0.01
                            ! margin can be classified across it. Tier B
                            ! (production gate set, 1001) cross-checks the
                            ! envelope: its worst cases match the probe's
                            ! bin envelope on both grids.
                            if (beak_pass) then
                                call compute_star_convexity_optimum_standalone_s( &
                                        params, N_RHO_JULY, zt7, g7, stg)
                                if (stg /= SHAPE_VALID) cycle  ! rho-rejected
                            else
                                th1(1) = 0.5_rk * PI_C
                                call compute_conversion_diagnostic_standalone_s( &
                                        params, th1, N_RHO_JULY, r1, dr1, zt7, g7, stg)
                                if (stg /= SHAPE_VALID &
                                        .and. stg /= FOS_ERROR_CONVERGENCE) cycle
                            end if

                            if (g7 >= G_BIN_LO + real(N_G_BINS - 1_ik, rk) * G_BIN_W) cycle

                            s_ref = compute_reference_surface_f(params)
                            call evaluate_shape_quality_diag_s(params, N_POINTS_SWEEP, &
                                    s_ref, dv_rel, ds_rel, rt_max, rp_max, rpp_max, &
                                    zs_d, g_d, stq)
                            if (stq /= SHAPE_VALID &
                                    .and. stq /= FOS_ERROR_CONVERGENCE) cycle

                            gb = int(floor((g7 - G_BIN_LO) / G_BIN_W), ik) + 1_ik
                            gb = min(max(gb, 1_ik), N_G_BINS)
                            g_count(gb) = g_count(gb) + 1_ik
                            if (stq == FOS_ERROR_CONVERGENCE) then
                                g_fail(gb) = g_fail(gb) + 1_ik
                            else
                                g_dv(gb) = max(g_dv(gb), abs(dv_rel))
                                g_ds(gb) = max(g_ds(gb), abs(ds_rel))
                                g_rt(gb) = max(g_rt(gb), rt_max)
                                g_rp(gb) = max(g_rp(gb), rp_max)
                                g_rpp(gb) = max(g_rpp(gb), rpp_max)
                            end if

                            ! f_min bins: rho- and star-passing shapes only.
                            if (g7 <= -STAR_CONVEXITY_MARGIN) then
                                fb = f_bin_f(f_min)
                                br = 1_ik              ! polar (boundary clamp)
                                if (interior_min) br = 2_ik   ! neck (interior)
                                f_count(fb, br) = f_count(fb, br) + 1_ik
                                if (stq == FOS_ERROR_CONVERGENCE) then
                                    f_fail(fb, br) = f_fail(fb, br) + 1_ik
                                else
                                    f_dv(fb, br) = max(f_dv(fb, br), abs(dv_rel))
                                    f_ds(fb, br) = max(f_ds(fb, br), abs(ds_rel))
                                    f_rt(fb, br) = max(f_rt(fb, br), rt_max)
                                    f_rp(fb, br) = max(f_rp(fb, br), rp_max)
                                    f_rpp(fb, br) = max(f_rpp(fb, br), rpp_max)
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        call system_clock(t_end)

        write(*, '(A)') 'g(s*) bins (July parity, diagnostic conversion):'
        write(*, '(A)') '   g_lo      count   fail  worst|dV/V|  worst|dS/S|' // &
                '  worst rt     worst R''     worst R'''''
        do b = 1_ik, N_G_BINS
            if (g_count(b) == 0_ik) cycle
            write(*, '(F8.3,I10,I7,5ES13.4)') G_BIN_LO + real(b - 1_ik, rk) * G_BIN_W, &
                    g_count(b), g_fail(b), g_dv(b), g_ds(b), g_rt(b), g_rp(b), g_rpp(b)
        end do

        write(*, '(A)') 'f_min bins (6/decade, branch: P = polar clamp, N = interior neck):'
        write(*, '(A)') '   f_min_lo   br     count   fail  worst|dV/V|  worst|dS/S|' // &
                '  worst rt     worst R''     worst R'''''
        do b = 0_ik, N_F_BINS
            do br = 1_ik, 2_ik
                if (f_count(b, br) == 0_ik) cycle
                write(*, '(ES11.3,A4,I10,I7,5ES13.4)') f_edge_f(b), &
                        merge('   P', '   N', br == 1_ik), f_count(b, br), &
                        f_fail(b, br), f_dv(b, br), f_ds(b, br), f_rt(b, br), &
                        f_rp(b, br), f_rpp(b, br)
            end do
        end do

        call print_cliff_s('Task-4 baselined tolerances', &
                TOL_VOLUME_B, TOL_SURFACE_REL, TOL_ROUND_TRIP, f_count, f_fail, &
                f_dv, f_ds, f_rt)
        call print_cliff_s('wmmm-inherited tolerances', &
                2.0e-9_rk, 1.0e-3_rk, 1.0e-9_rk, f_count, f_fail, f_dv, f_ds, f_rt)

        write(*, '(A,F10.2,A)') '  Probe wall time: ', &
                real(t_end - t_start, rk) / real(t_rate, rk), ' s'

    end subroutine run_probe_s

    !> Log-spaced bin index: 0 = underflow (< 1e-6), 1..30 over [1e-6, 1e-1],
    !! overflow clamps into bin 30 (6 bins per decade).
    pure function f_bin_f(f) result(bin)

        real(kind = rk), intent(in) :: f
        integer(kind = ik) :: bin

        if (f < 1.0e-6_rk) then
            bin = 0_ik
        else
            bin = int(floor(6.0_rk * (log10(f) + 6.0_rk)), ik) + 1_ik
            bin = min(bin, 30_ik)
        end if

    end function f_bin_f

    !> Lower edge of f bin b (b = 0 prints the underflow marker 0).
    pure function f_edge_f(b) result(edge)

        integer(kind = ik), intent(in) :: b
        real(kind = rk) :: edge

        if (b == 0_ik) then
            edge = 0.0_rk
        else
            edge = 1.0e-6_rk * 10.0_rk**(real(b - 1_ik, rk) / 6.0_rk)
        end if

    end function f_edge_f

    !> Descending-f_min cliff readout per branch against a tolerance triple:
    !! the first bin (from high f_min down) whose worst metrics leave
    !! tolerance or that records conversion failures.
    subroutine print_cliff_s(label, tol_v, tol_s, tol_r, f_count, f_fail, &
            f_dv, f_ds, f_rt)

        character(len = *), intent(in) :: label
        real(kind = rk), intent(in) :: tol_v, tol_s, tol_r
        integer(kind = ik), intent(in) :: f_count(0:, :), f_fail(0:, :)
        real(kind = rk), intent(in) :: f_dv(0:, :), f_ds(0:, :), f_rt(0:, :)

        character(len = 5), parameter :: BR_NAME(2) = ['polar', 'neck ']
        integer(kind = ik) :: bb, brr, cliff, n_bins

        n_bins = ubound(f_count, 1)
        write(*, '(A,A,A,3ES10.2,A)') 'cliff readout [', label, '] (dv,ds,rt <= ', &
                tol_v, tol_s, tol_r, '):'
        do brr = 1_ik, 2_ik
            cliff = -1_ik
            do bb = n_bins, 0_ik, -1_ik
                if (f_count(bb, brr) == 0_ik) cycle
                if (f_fail(bb, brr) > 0_ik .or. f_dv(bb, brr) > tol_v &
                        .or. f_ds(bb, brr) > tol_s .or. f_rt(bb, brr) > tol_r) then
                    cliff = bb
                    exit
                end if
            end do
            if (cliff < 0_ik) then
                write(*, '(A,A,A)') '  cliff(', trim(BR_NAME(brr)), &
                        '): none — all populated bins within tolerance'
            else
                write(*, '(A,A,A,ES11.3,A,ES11.3,A)') '  cliff(', &
                        trim(BR_NAME(brr)), '): f_min in [', f_edge_f(cliff), &
                        ', ', f_edge_f(cliff + 1_ik), ')'
                write(*, '(A,ES11.3)') '  proposed threshold (one bin above): ', &
                        f_edge_f(cliff + 1_ik)
            end if
        end do

    end subroutine print_cliff_s

end program fos_param_geometry_sweep_test
