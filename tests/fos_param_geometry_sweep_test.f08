!> Geometry sweep over the FoS parameter box, ported from WMMM
!! tests/fos_geometry_validation_test.f08 and adapted to the 2.0 tier-1 API.
!!
!! Tier A classifies every grid point through the production gates
!! (`compute_shape_standalone_s`) and pins the verdict counts against golden
!! tallies captured on the default coarse grid. Tier B (accuracy) and the
!! representability probe are added by later tasks.
!!
!! Grid: c in [0.75, 3.50], a3 in [0.00, 0.60], a4 in [-0.20, 0.75],
!! a5/a6 in [-0.15, 0.15]; default step 0.05 on every axis (713,440 points),
!! `--full` refines a5/a6 to 0.01 (~14.0M points).
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
    real(kind = rk), parameter :: A6_LO = -0.15_rk, A6_HI = 0.15_rk

    !> u-grid resolution of every sweep-tier conversion (wmmm N_RHO_INTERNAL
    !! parity).
    integer(kind = ik), parameter :: N_POINTS_SWEEP = 1001_ik

    type :: sweep_config_t
        real(kind = rk) :: dc = 0.05_rk, da3 = 0.05_rk, da4 = 0.05_rk
        real(kind = rk) :: da5 = 0.05_rk, da6 = 0.05_rk
        integer(kind = ik) :: n_c = 0_ik, n_a3 = 0_ik, n_a4 = 0_ik
        integer(kind = ik) :: n_a5 = 0_ik, n_a6 = 0_ik
        logical :: probe = .false., skip_tier_b = .false.
        !> .true. iff no step override and no --full: golden tallies apply.
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
            case default
                write(*, '(A,A)') 'Unknown argument: ', trim(arg)
                write(*, '(A)') 'Usage: fos_param_geometry_sweep_test [--full]' // &
                        ' [--probe] [--skip-tier-b]' // &
                        ' [--dc X] [--da3 X] [--da4 X] [--da5 X] [--da6 X]'
                stop 2
            end select
            i = i + 1_ik
        end do

        cfg%n_c = count_steps_f(C_LO, C_HI, cfg%dc)
        cfg%n_a3 = count_steps_f(A3_LO, A3_HI, cfg%da3)
        cfg%n_a4 = count_steps_f(A4_LO, A4_HI, cfg%da4)
        cfg%n_a5 = count_steps_f(A5_LO, A5_HI, cfg%da5)
        cfg%n_a6 = count_steps_f(A6_LO, A6_HI, cfg%da6)

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
                A5_LO + real(i5 - 1_ik, rk) * cfg%da5, &
                A6_LO + real(i6 - 1_ik, rk) * cfg%da6, &
                0.0_rk, 0.0_rk]

    end function grid_params_f

end module fos_sweep_support_mod

program fos_param_geometry_sweep_test

    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: compute_shape_standalone_s, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
            FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY
    use shape_core_mod, only: SHAPE_VALID
    use fos_sweep_support_mod, only: sweep_config_t, parse_args_s, &
            print_config_s, grid_params_f, N_POINTS_SWEEP
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    !> Tier A verdict counts on the default coarse grid, captured 2026-08-11
    !! (Task 3 first run). Index: 0 = valid, 1..4 = codes 100..103, 5 = other.
    integer(kind = ik), parameter :: GOLDEN_TALLY(0:5) = &
            [468390_ik, 0_ik, 1058_ik, 0_ik, 243992_ik, 0_ik]

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

end program fos_param_geometry_sweep_test
