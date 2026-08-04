!> Contract family 3: interleaving. Alternating between the five cached outputs
!! must not change any of them. Each output is driven on a warm cache that has
!! just produced the other four, and compared — bit for bit — against a FRESH
!! cache asked for that output alone.
!!
!! This is the family that catches shared-buffer bugs the bitwise family cannot:
!! there, one output is asked for repeatedly, so a scratch array aliased between
!! two outputs is never observed. Here `cache_shape_s`, `cache_radius_grid_s`,
!! `cache_rho_z_grid_s`, `cache_radius_and_derivative_s` and `cache_neck_s` take
!! turns over a six-vector parameter walk, and every one of them must still
!! answer what a cold cache answers.
!!
!! Two further properties, both about the at-thetas form:
!!
!!   - fed the cache's OWN theta nodes it reproduces `cache_radius_and_derivative_s`
!!     bit for bit — same bundle, same cosines, same solver;
!!   - it stamps nothing, so it cannot move a single recompute counter.
!!
!! And the recovery path: a beak-rejected vector in the middle of the walk
!! returns the engine to cold, and the next valid vector still matches cold.
program fos_param_interleave_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_radius_grid_s, cache_rho_z_grid_s, &
            cache_radius_and_derivative_s, cache_radius_and_derivative_at_thetas_s, &
            cache_neck_s, cache_recompute_count_f, &
            SHAPE_VALID, FOS_ERROR_BEAK_SINGULARITY
    use test_utils_mod, only: assert_true, assert_int_eq, assert_bits_eq, &
            test_summary

    implicit none

    integer(kind = ik), parameter :: N_PARAMS = 8_ik
    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 64_ik
    integer(kind = ik), parameter :: N_WALK = 6_ik
    integer(kind = ik), parameter :: N_TRACKED = 8_ik

    real(kind = rk), parameter :: BASE8(N_PARAMS) = &
            [1.6_rk, 0.12_rk, 0.08_rk, 0.05_rk, 0.03_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: walk(N_PARAMS, N_WALK)
    real(kind = rk)    :: beak_params(N_PARAMS)
    integer(kind = ik) :: i, j, w
    logical            :: any_neck

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do

    !---------------------------------------------------------------------------
    ! The walk: one parameter moves per step, so consecutive vectors invalidate
    ! different mask subsets. The last one is a deeply necked symmetric shape,
    ! so the neck output has something to find.
    !---------------------------------------------------------------------------
    do j = 1_ik, N_PARAMS
        walk(j, 1) = BASE8(j)
    end do

    ! Vectors 2-5 each carry every earlier step plus one new one, so consecutive
    ! diffs hit different mask subsets: c, then an odd coefficient, then two
    ! even ones.
    do w = 2_ik, N_WALK - 1_ik
        do j = 1_ik, N_PARAMS
            walk(j, w) = walk(j, w - 1_ik)
        end do
        select case (w)
        case (2_ik)
            walk(1, w) = walk(1, w) + 0.1_rk        ! c
        case (3_ik)
            walk(2, w) = walk(2, w) + 0.03_rk       ! a3, odd
        case (4_ik)
            walk(3, w) = walk(3, w) + 0.02_rk       ! a4, even
        case (5_ik)
            walk(5, w) = walk(5, w) - 0.01_rk       ! a6, even
        end select
    end do

    ! v6: a symmetric, deeply necked shape — the one vector of the walk whose
    ! cylindrical profile actually HAS a neck.
    do j = 1_ik, N_PARAMS
        walk(j, N_WALK) = 0.0_rk
    end do
    walk(1, N_WALK) = 2.0_rk
    walk(3, N_WALK) = 0.4_rk

    call probe_beak_s(beak_params)

    call run_walk_s(any_neck)
    call assert_true(any_neck, 'the walk exercises a shape that HAS a neck')

    call run_at_thetas_s()
    call run_recovery_s()

    call test_summary()

contains

    !> The interleaved walk: one warm cache, all five outputs per vector, each
    !! checked against a fresh cold cache asked for that output alone.
    subroutine run_walk_s(neck_seen)

        logical, intent(out) :: neck_seen

        type(cache_t) :: warm
        integer(kind = ik) :: step, status
        logical :: found
        character(len = 32) :: label

        neck_seen = .false.

        call cache_init_s(warm, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'walk: warm cache init')

        do step = 1_ik, N_WALK
            write(label, '(A,I0)') 'walk vector ', step
            call interleave_step_s(warm, walk(:, step), trim(label), found)
            if (found) neck_seen = .true.
        end do

        call cache_free_s(warm)

    end subroutine run_walk_s

    !> One vector: shape -> radius_grid -> rho_z_grid -> radius_and_derivative
    !! -> neck on the warm cache, each answer compared bitwise with cold.
    subroutine interleave_step_s(warm, params, label, neck_found)

        type(cache_t), intent(inout) :: warm
        real(kind = rk), intent(in) :: params(N_PARAMS)
        character(len = *), intent(in) :: label
        logical, intent(out) :: neck_found

        type(cache_t) :: cold
        real(kind = rk) :: zs_w, rn_w, rs_w, zs_c, rn_c, rs_c
        real(kind = rk) :: radii_w(N_THETA), radii_c(N_THETA)
        real(kind = rk) :: rd_w(N_THETA), rd_c(N_THETA)
        real(kind = rk) :: dr_w(N_THETA), dr_c(N_THETA)
        real(kind = rk) :: z_w(N_POINTS), rho_w(N_POINTS), drho_w(N_POINTS)
        real(kind = rk) :: z_c(N_POINTS), rho_c(N_POINTS), drho_c(N_POINTS)
        real(kind = rk) :: gshift_w, gshift_c
        real(kind = rk) :: zn_w, rhon_w, zn_c, rhon_c
        integer(kind = ik) :: sw, sc, n
        logical :: found_w, found_c

        ! --- shape ---
        call cache_shape_s(warm, params, zs_w, rn_w, rs_w, sw)
        call fresh_cache_s(cold, label // ': cold shape init')
        call cache_shape_s(cold, params, zs_c, rn_c, rs_c, sc)
        call cache_free_s(cold)
        call assert_int_eq(sw, sc, label // ': shape status warm == cold')
        call assert_bits_eq(zs_w, zs_c, label // ': shape z_shift bitwise')
        call assert_bits_eq(rn_w, rn_c, label // ': shape r_north bitwise')
        call assert_bits_eq(rs_w, rs_c, label // ': shape r_south bitwise')

        ! --- radius grid ---
        call cache_radius_grid_s(warm, params, radii_w, sw)
        call fresh_cache_s(cold, label // ': cold radius grid init')
        call cache_radius_grid_s(cold, params, radii_c, sc)
        call cache_free_s(cold)
        call assert_int_eq(sw, sc, label // ': radius grid status warm == cold')
        do n = 1_ik, N_THETA
            call assert_bits_eq(radii_w(n), radii_c(n), &
                    label // ': radius grid node bitwise')
        end do

        ! --- cylindrical grid ---
        call cache_rho_z_grid_s(warm, params, z_w, rho_w, drho_w, gshift_w, sw)
        call fresh_cache_s(cold, label // ': cold rho_z_grid init')
        call cache_rho_z_grid_s(cold, params, z_c, rho_c, drho_c, gshift_c, sc)
        call cache_free_s(cold)
        call assert_int_eq(sw, sc, label // ': rho_z_grid status warm == cold')
        call assert_bits_eq(gshift_w, gshift_c, label // ': rho_z_grid z_shift bitwise')
        do n = 1_ik, N_POINTS
            call assert_bits_eq(z_w(n), z_c(n), label // ': rho_z_grid z node bitwise')
            call assert_bits_eq(rho_w(n), rho_c(n), label // ': rho_z_grid rho node bitwise')
            call assert_bits_eq(drho_w(n), drho_c(n), &
                    label // ': rho_z_grid drho_dz node bitwise')
        end do

        ! --- radius and derivative ---
        call cache_radius_and_derivative_s(warm, params, rd_w, dr_w, sw)
        call fresh_cache_s(cold, label // ': cold derivative init')
        call cache_radius_and_derivative_s(cold, params, rd_c, dr_c, sc)
        call cache_free_s(cold)
        call assert_int_eq(sw, sc, label // ': derivative status warm == cold')
        do n = 1_ik, N_THETA
            call assert_bits_eq(rd_w(n), rd_c(n), label // ': derivative radius bitwise')
            call assert_bits_eq(dr_w(n), dr_c(n), label // ': dR/dtheta bitwise')
        end do

        ! The radius-only grid and the derivative grid are separate
        ! intermediates, but they are the same solve: their radii must agree.
        if (sw == SHAPE_VALID) then
            do n = 1_ik, N_THETA
                call assert_bits_eq(radii_w(n), rd_w(n), &
                        label // ': #7 and #8 radii agree bitwise')
            end do
        end if

        ! --- neck ---
        call cache_neck_s(warm, params, zn_w, rhon_w, found_w, sw)
        call fresh_cache_s(cold, label // ': cold neck init')
        call cache_neck_s(cold, params, zn_c, rhon_c, found_c, sc)
        call cache_free_s(cold)
        call assert_int_eq(sw, sc, label // ': neck status warm == cold')
        call assert_true(found_w .eqv. found_c, label // ': neck verdict warm == cold')
        call assert_bits_eq(zn_w, zn_c, label // ': z_neck bitwise')
        call assert_bits_eq(rhon_w, rhon_c, label // ': rho_neck bitwise')

        neck_found = found_w .and. sw == SHAPE_VALID

    end subroutine interleave_step_s

    !> The at-thetas form on the cache's own nodes: bit-identical to the fixed
    !! grid, and invisible to every counter.
    subroutine run_at_thetas_s()

        type(cache_t) :: cache
        real(kind = rk) :: params(N_PARAMS)
        real(kind = rk) :: radii(N_THETA), dr_dtheta(N_THETA)
        real(kind = rk) :: at_radii(N_THETA), at_dr(N_THETA)
        integer(kind = ik) :: before(N_TRACKED)
        integer(kind = ik) :: j, n, status

        call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'at-thetas cache init')

        do j = 1_ik, N_PARAMS
            params(j) = BASE8(j)
        end do

        call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)
        call assert_int_eq(status, SHAPE_VALID, 'at-thetas: fixed-grid call valid')

        do n = 1_ik, N_TRACKED
            before(n) = int(cache_recompute_count_f(cache, n), ik)
        end do

        call cache_radius_and_derivative_at_thetas_s(cache, params, thetas, &
                at_radii, at_dr, status)
        call assert_int_eq(status, SHAPE_VALID, 'at-thetas call valid')

        do n = 1_ik, N_THETA
            call assert_bits_eq(at_radii(n), radii(n), &
                    'at_thetas radius == fixed grid, bitwise')
            call assert_bits_eq(at_dr(n), dr_dtheta(n), &
                    'at_thetas dR/dtheta == fixed grid, bitwise')
        end do

        do n = 1_ik, N_TRACKED
            call assert_int_eq(int(cache_recompute_count_f(cache, n), ik), before(n), &
                    'at_thetas stamps nothing')
        end do

        call cache_free_s(cache)

    end subroutine run_at_thetas_s

    !> A beak rejection between two valid vectors leaves nothing behind.
    subroutine run_recovery_s()

        type(cache_t) :: warm
        real(kind = rk) :: radii(N_THETA)
        integer(kind = ik) :: status
        logical :: found

        call cache_init_s(warm, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'recovery: warm cache init')

        call interleave_step_s(warm, walk(:, 1), 'recovery, before', found)

        call cache_radius_grid_s(warm, beak_params, radii, status)
        call assert_int_eq(status, FOS_ERROR_BEAK_SINGULARITY, &
                'recovery: interposed vector is beak-rejected')

        call interleave_step_s(warm, walk(:, 1), 'recovery, after', found)
        call interleave_step_s(warm, walk(:, 3), 'recovery, stepped after', found)

        call cache_free_s(warm)

    end subroutine run_recovery_s

    !> A beak-rejected vector, probed rather than assumed: a4 walked up from a
    !! symmetric base until the beak gate fires.
    subroutine probe_beak_s(rejected)

        real(kind = rk), intent(out) :: rejected(N_PARAMS)

        type(cache_t) :: probe
        real(kind = rk) :: p(N_PARAMS), z_shift, r_north, r_south
        integer(kind = ik) :: j, k, status
        logical :: hit

        call cache_init_s(probe, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, 'beak probe cache init')

        do j = 1_ik, N_PARAMS
            rejected(j) = 0.0_rk
        end do
        rejected(1) = 2.0_rk

        hit = .false.
        do k = 0_ik, 300_ik
            do j = 1_ik, N_PARAMS
                p(j) = 0.0_rk
            end do
            p(1) = 2.0_rk
            p(3) = 0.60_rk + 5.0e-4_rk * real(k, rk)
            call cache_shape_s(probe, p, z_shift, r_north, r_south, status)
            if (status == FOS_ERROR_BEAK_SINGULARITY) then
                do j = 1_ik, N_PARAMS
                    rejected(j) = p(j)
                end do
                hit = .true.
                exit
            end if
        end do

        call assert_true(hit, 'probe found a beak-rejected vector')

        call cache_free_s(probe)

    end subroutine probe_beak_s

    !> A fresh cold cache with the suite's init arguments.
    subroutine fresh_cache_s(cold, label)
        type(cache_t), intent(out) :: cold
        character(len = *), intent(in) :: label
        integer(kind = ik) :: status
        call cache_init_s(cold, N_PARAMS, N_POINTS, thetas, status)
        call assert_int_eq(status, SHAPE_VALID, label)
    end subroutine fresh_cache_s

end program fos_param_interleave_test
