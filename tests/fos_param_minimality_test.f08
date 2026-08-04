!> Contract family 2: minimality. Every cached compute must recompute exactly
!! the intermediates its parameter diff invalidated AND its output needs — no
!! more, no fewer. The engine's eight per-intermediate recompute counters are
!! the observable; this suite pins the whole table, checkpoint by checkpoint,
!! as an EXACT delta set rather than a spot check on one index.
!!
!! The dependency masks under test (bit j-1 set = depends on params(j)):
!!
!!   #1 a2        p3, p5, p7                 (even coefficients only)
!!   #2 z_shift   p1, p2, p4, p6, p8         (c and the odd coefficients)
!!   #3 f_grid    p2..p8
!!   #4 beak      p2..p8
!!   #5 rho_grid  all
!!   #6 resolve   all
!!   #7 radii     all
!!   #8 radii+dR  all
!!
!! Two needs-lists are asserted alongside them: the cylindrical output stops at
!! [1,2,3,5] (no beak scan, no resolve), and the ungated star-convexity optimum
!! reads the STORED resolve, so after a successful `cache_shape_s` it recomputes
!! nothing at all — but after a 101-rejected one, which returns the engine to
!! cold, it recomputes #1-#6 from scratch.
program fos_param_minimality_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_rho_z_grid_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_star_convexity_optimum_s, &
            cache_recompute_count_f, &
            SHAPE_VALID, FOS_ERROR_NOT_STAR_CONVEX
    use test_utils_mod, only: assert_true, assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_PARAMS = 8_ik
    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 64_ik
    integer(kind = ik), parameter :: N_TRACKED = 8_ik

    real(kind = rk), parameter :: BASE8(N_PARAMS) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Odd-heavy family: raising c drives g(s*) up through the star-convexity
    !! margin, which is how the 101 vector is probed.
    real(kind = rk), parameter :: ODD_HEAVY(N_PARAMS) = &
            [1.0_rk, -0.55_rk, 0.03_rk, -0.02_rk, 0.115_rk, 0.078_rk, 0.085_rk, 0.0_rk]

    character(len = 12), parameter :: NAMES(N_TRACKED) = &
            ['#1 a2       ', '#2 z_shift  ', '#3 f_grid   ', '#4 beak     ', &
             '#5 rho_grid ', '#6 resolve  ', '#7 radii    ', '#8 radii+dR ']

    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: params_101(N_PARAMS)
    integer(kind = ik) :: i

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do

    call probe_101_s(params_101)

    call run_radius_table_s()
    call run_cylindrical_table_s()
    call run_optimum_after_shape_s()
    call run_optimum_after_rejection_s(params_101)

    call test_summary()

contains

    !> The counter table of the R(theta) path, one perturbation class at a time.
    subroutine run_radius_table_s()

        type(cache_t) :: cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii(N_THETA), dr_dtheta(N_THETA)
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: j, status

        call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'radius-table cache init')

        do j = 1_ik, N_PARAMS
            params(j) = BASE8(j)
        end do

        ! Cold radius grid: everything the grid needs runs once; the derivative
        ! table (#8) is a separate intermediate and stays cold.
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold radius grid status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 0_ik]
        call check_counters_s(cache, want, 'cold radius grid')

        ! Repeat on the same vector: nothing is recomputed.
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'repeat radius grid status')
        call check_counters_s(cache, want, 'repeat radius grid')

        ! c-only step: #1 a2 and the f-grid/beak pair do not depend on c.
        params(1) = BASE8(1) + 0.1_rk
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'c-step radius grid status')
        want = [1_ik, 2_ik, 1_ik, 1_ik, 2_ik, 2_ik, 2_ik, 0_ik]
        call check_counters_s(cache, want, 'c-only step')

        ! Even coefficient a4 = params(3): the intrinsic z-shift is blind to it.
        params(3) = BASE8(3) + 0.02_rk
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'a4-step radius grid status')
        want = [2_ik, 2_ik, 2_ik, 2_ik, 3_ik, 3_ik, 3_ik, 0_ik]
        call check_counters_s(cache, want, 'even a4 step')

        ! Odd coefficient a3 = params(2): the volume-fixed a2 is blind to it.
        params(2) = BASE8(2) + 0.02_rk
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'a3-step radius grid status')
        want = [2_ik, 3_ik, 3_ik, 3_ik, 4_ik, 4_ik, 4_ik, 0_ik]
        call check_counters_s(cache, want, 'odd a3 step')

        ! Repeat again: the empty delta set, on a vector three steps from base.
        call cache_radius_grid_s(cache, params, radii, status)
        call assert_int_eq(status, SHAPE_VALID, 'second repeat status')
        call check_counters_s(cache, want, 'repeat after the walk')

        ! Derivative on the same vector: only #8 is missing.
        call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)
        call assert_int_eq(status, SHAPE_VALID, 'warm derivative status')
        want = [2_ik, 3_ik, 3_ik, 3_ik, 4_ik, 4_ik, 4_ik, 1_ik]
        call check_counters_s(cache, want, 'radius_and_derivative after radius_grid')

        call cache_free_s(cache)

    end subroutine run_radius_table_s

    !> The cylindrical output's needs-list, on a FRESH cache: a warm one cannot
    !! show where an entry point STOPS, only that it recomputes nothing.
    subroutine run_cylindrical_table_s()

        type(cache_t) :: cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
        real(kind = rk) :: z_shift
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: j, status

        call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'cylindrical cache init')

        do j = 1_ik, N_PARAMS
            params(j) = BASE8(j)
        end do

        ! Cold: [1, 2, 3, 5]. The beak scan (#4) gates the R(theta) conversion,
        ! not the cylindrical profile, so it is neither computed nor stamped.
        call cache_rho_z_grid_s(cache, params, z, rho, drho_dz, z_shift, status)
        call assert_int_eq(status, SHAPE_VALID, 'cold rho_z_grid status')
        want = [1_ik, 1_ik, 1_ik, 0_ik, 1_ik, 0_ik, 0_ik, 0_ik]
        call check_counters_s(cache, want, 'cold rho_z_grid')

        ! c-step: f(u) does not depend on c, the rho grid does.
        params(1) = BASE8(1) + 0.1_rk
        call cache_rho_z_grid_s(cache, params, z, rho, drho_dz, z_shift, status)
        call assert_int_eq(status, SHAPE_VALID, 'c-step rho_z_grid status')
        want = [1_ik, 2_ik, 1_ik, 0_ik, 2_ik, 0_ik, 0_ik, 0_ik]
        call check_counters_s(cache, want, 'rho_z_grid after a c-step')

        call cache_free_s(cache)

    end subroutine run_cylindrical_table_s

    !> The ungated optimum reads the stored resolve: after a successful shape
    !! call on the same vector it recomputes nothing.
    subroutine run_optimum_after_shape_s()

        type(cache_t) :: cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: z_shift, r_north, r_south, z_shift_total, g_opt
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: j, status

        call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'optimum-after-shape cache init')

        do j = 1_ik, N_PARAMS
            params(j) = BASE8(j)
        end do

        call cache_shape_s(cache, params, z_shift, r_north, r_south, status)
        call assert_int_eq(status, SHAPE_VALID, 'shape call status')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 0_ik, 0_ik]
        call check_counters_s(cache, want, 'cold shape call')

        call cache_star_convexity_optimum_s(cache, params, z_shift_total, g_opt, status)
        call assert_int_eq(status, SHAPE_VALID, 'optimum after shape status')
        call check_counters_s(cache, want, 'optimum after a successful shape')
        call assert_true(g_opt < 0.0_rk, 'accepted shape has a negative g(s*)')

        call cache_free_s(cache)

    end subroutine run_optimum_after_shape_s

    !> A 101 rejection returns the engine to cold, so the very next call — even
    !! on the same vector — recomputes #1-#6 from scratch. The optimum itself is
    !! ungated, so it answers where the shape call refused to.
    subroutine run_optimum_after_rejection_s(params)

        real(kind = rk), intent(in) :: params(N_PARAMS)

        type(cache_t) :: cache
        real(kind = rk) :: z_shift, r_north, r_south, z_shift_total, g_opt
        integer(kind = ik) :: want(N_TRACKED)
        integer(kind = ik) :: status

        call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'optimum-after-101 cache init')

        call cache_shape_s(cache, params, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_ERROR_NOT_STAR_CONVEX, &
                'probe vector is rejected with 101')
        want = [1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 1_ik, 0_ik, 0_ik]
        call check_counters_s(cache, want, 'the 101-rejected shape call')

        call cache_star_convexity_optimum_s(cache, params, z_shift_total, g_opt, status)
        call assert_int_eq(status, SHAPE_VALID, 'optimum reports on the 101 shape')
        want = [2_ik, 2_ik, 2_ik, 2_ik, 2_ik, 2_ik, 0_ik, 0_ik]
        call check_counters_s(cache, want, 'optimum after a 101-rejected shape')

        call cache_free_s(cache)

    end subroutine run_optimum_after_rejection_s

    !> First vector of the odd-heavy family the star-convexity margin rejects.
    subroutine probe_101_s(found_params)

        real(kind = rk), intent(out) :: found_params(N_PARAMS)

        type(cache_t) :: probe
        real(kind = rk) :: p(N_PARAMS), z_shift, r_north, r_south
        integer(kind = ik) :: j, k, status
        logical :: hit

        call cache_init_s(probe, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, '101 probe cache init')

        do j = 1_ik, N_PARAMS
            found_params(j) = ODD_HEAVY(j)
        end do

        hit = .false.
        do k = 0_ik, 60_ik
            do j = 1_ik, N_PARAMS
                p(j) = ODD_HEAVY(j)
            end do
            p(1) = 1.3_rk + 0.05_rk * real(k, rk)
            call cache_shape_s(probe, p, z_shift, r_north, r_south, status)
            if (status == FOS_ERROR_NOT_STAR_CONVEX) then
                do j = 1_ik, N_PARAMS
                    found_params(j) = p(j)
                end do
                hit = .true.
                exit
            end if
        end do

        call assert_true(hit, 'probe found a 101-rejected vector')

        call cache_free_s(probe)

    end subroutine probe_101_s

    !> Assert all eight recompute counters at one checkpoint.
    subroutine check_counters_s(cache, expected, label)

        type(cache_t), intent(in) :: cache
        integer(kind = ik), intent(in) :: expected(N_TRACKED)
        character(len = *), intent(in) :: label

        integer(kind = ik) :: k

        do k = 1_ik, N_TRACKED
            call assert_int_eq(int(cache_recompute_count_f(cache, k), ik), &
                    expected(k), trim(NAMES(k)) // ' count after ' // label)
        end do

    end subroutine check_counters_s

end program fos_param_minimality_test
