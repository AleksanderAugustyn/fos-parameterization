!> Standalone suite: the tier-1 forms — no engine, no cache, no tables handed
!! to the caller. One call in, one answer out, everything discarded.
!!
!! Two properties are the subject:
!!
!!   - **N_max = 50.** Tier-1 sizes its tables per call, so it serves parameter
!!     vectors the 8-parameter cached tier cannot. 50 is accepted, 51 is
!!     SHAPE_ERROR_TOO_MANY_PARAMS.
!!   - **Standalone == cold cache, BITWISE.** Both tiers drive the same kernels
!!     in the same cold order, so at matching (n_points, thetas) the outputs
!!     must agree to the last bit, not to a tolerance. Asserted with a zero
!!     tolerance on every one of the six forms.
!!
!! Degenerate coverage is per entry point: all six forms are driven with an
!! empty parameter vector and with a degenerate c, and each must report 102 and
!! zero-fill every output it owns.
program fos_param_standalone_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: FOS_MAX_PARAMS, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, &
            compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, &
            compute_neck_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            cache_t, cache_init_s, cache_free_s, &
            cache_shape_s, cache_rho_z_grid_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_neck_s, &
            cache_star_convexity_optimum_s, &
            FOS_ERROR_INVALID_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, &
            test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_THETA = 64_ik
    integer(kind = ik), parameter :: N_PARAMS = 7_ik

    real(kind = rk), parameter :: PARAMS7(N_PARAMS) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Symmetric (c, a4) necked shape, from the neck suite's family.
    real(kind = rk), parameter :: NECK7(N_PARAMS) = &
            [2.0_rk, 0.0_rk, 0.4_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]

    !> Written into every buffer before a rejection test: a zero-filled output
    !! is evidence only if the buffer did not already hold zeros.
    real(kind = rk), parameter :: SENTINEL = 7.0_rk

    real(kind = rk), parameter :: NO_PARAMS(0) = [real(kind = rk) ::]

    type(cache_t)      :: cache
    integer(kind = ik) :: status, cstatus, i
    real(kind = rk)    :: thetas(N_THETA)
    real(kind = rk)    :: params50(50), params51(51), params_bad_c(1)
    real(kind = rk)    :: radii(N_THETA), dr(N_THETA)
    real(kind = rk)    :: cradii(N_THETA), cdr(N_THETA)
    real(kind = rk)    :: r50(N_THETA), dr50(N_THETA)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho(N_POINTS)
    real(kind = rk)    :: cz(N_POINTS), crho(N_POINTS), cdrho(N_POINTS)
    real(kind = rk)    :: z_shift, r_north, r_south
    real(kind = rk)    :: cz_shift, cr_north, cr_south
    real(kind = rk)    :: z_neck, rho_neck, cz_neck, crho_neck
    real(kind = rk)    :: g_opt, z_shift_total, cg_opt, cz_shift_total
    logical            :: found, cfound

    do i = 1_ik, N_THETA
        thetas(i) = real(i - 1_ik, rk) * PI_C / real(N_THETA - 1_ik, rk)
    end do

    !---------------------------------------------------------------------------
    ! N_max
    !---------------------------------------------------------------------------
    call assert_int_eq(FOS_MAX_PARAMS, 50_ik, 'FOS_MAX_PARAMS = 50')

    !---------------------------------------------------------------------------
    ! 50 parameters are served; 51 are rejected by every entry point
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== N_max boundary ==='

    params50 = 1.0e-6_rk
    params50(1) = 1.5_rk
    params51(1:50) = params50
    params51(51) = 1.0e-6_rk

    ! A 50-vector is a shape question, not a length error: whatever the shape
    ! gates decide, the answer must not be "too many parameters".
    call compute_shape_standalone_s(params50, N_POINTS, z_shift, r_north, &
            r_south, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: shape not rejected on length')

    call compute_radius_grid_standalone_s(params50, thetas, N_POINTS, r50, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: radius grid not rejected on length')

    call compute_radius_and_derivative_standalone_s(params50, thetas, N_POINTS, &
            r50, dr50, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: radius+derivative not rejected on length')

    call compute_rho_z_grid_standalone_s(params50, N_POINTS, z, rho, drho, &
            z_shift, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: rho(z) grid not rejected on length')

    call compute_neck_standalone_s(params50, N_POINTS, z_neck, rho_neck, found, &
            status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: neck not rejected on length')

    call compute_star_convexity_optimum_standalone_s(params50, N_POINTS, &
            z_shift_total, g_opt, status)
    call assert_true(status /= SHAPE_ERROR_TOO_MANY_PARAMS, &
            '50 params: star optimum not rejected on length')

    call reject_all_forms(params51, SHAPE_ERROR_TOO_MANY_PARAMS, '51 params')

    !---------------------------------------------------------------------------
    ! Degenerate parameter vectors, per entry point
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Degenerate vectors, all six forms ==='

    call reject_all_forms(NO_PARAMS, FOS_ERROR_INVALID_C, 'empty params')

    params_bad_c(1) = 1.0e-11_rk
    call reject_all_forms(params_bad_c, FOS_ERROR_INVALID_C, 'degenerate c')

    !---------------------------------------------------------------------------
    ! Standalone == cold cache, bitwise, on every form
    !---------------------------------------------------------------------------
    write(*, '(A)') '=== Standalone vs cold cache (bitwise) ==='

    ! R(theta) on the fixed grid
    call compute_radius_grid_standalone_s(PARAMS7, thetas, N_POINTS, radii, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone radius grid valid')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache init valid')
    call cache_radius_grid_s(cache, PARAMS7, cradii, cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache radius grid valid')
    call assert_bits_eq(radii, cradii, 'radius grid')

    ! R(theta) and dR/dtheta
    call compute_radius_and_derivative_standalone_s(PARAMS7, thetas, N_POINTS, &
            radii, dr, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone radius+derivative valid')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call cache_radius_and_derivative_s(cache, PARAMS7, cradii, cdr, cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache radius+derivative valid')
    call assert_bits_eq(radii, cradii, 'radius (derivative form)')
    call assert_bits_eq(dr, cdr, 'dR/dtheta')

    ! Resolved shape
    call compute_shape_standalone_s(PARAMS7, N_POINTS, z_shift, r_north, &
            r_south, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone shape valid')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call cache_shape_s(cache, PARAMS7, cz_shift, cr_north, cr_south, cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache shape valid')
    call assert_abs_close(z_shift, cz_shift, 0.0_rk, 'shape z_shift bitwise')
    call assert_abs_close(r_north, cr_north, 0.0_rk, 'shape r_north bitwise')
    call assert_abs_close(r_south, cr_south, 0.0_rk, 'shape r_south bitwise')

    ! Cylindrical profile
    call compute_rho_z_grid_standalone_s(PARAMS7, N_POINTS, z, rho, drho, &
            z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone rho(z) grid valid')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call cache_rho_z_grid_s(cache, PARAMS7, cz, crho, cdrho, cz_shift, cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache rho(z) grid valid')
    call assert_bits_eq(z, cz, 'rho(z) grid z')
    call assert_bits_eq(rho, crho, 'rho(z) grid rho')
    call assert_bits_eq(drho, cdrho, 'rho(z) grid drho_dz')
    call assert_abs_close(z_shift, cz_shift, 0.0_rk, 'rho(z) grid z_shift bitwise')

    ! Neck (on a vector that has one)
    call compute_neck_standalone_s(NECK7, N_POINTS, z_neck, rho_neck, found, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone neck valid')
    call assert_true(found, 'standalone neck found')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call cache_neck_s(cache, NECK7, cz_neck, crho_neck, cfound, cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache neck valid')
    call assert_true(found .eqv. cfound, 'neck found matches')
    call assert_abs_close(z_neck, cz_neck, 0.0_rk, 'z_neck bitwise')
    call assert_abs_close(rho_neck, crho_neck, 0.0_rk, 'rho_neck bitwise')

    ! Ungated star-convexity optimum
    call compute_star_convexity_optimum_standalone_s(PARAMS7, N_POINTS, &
            z_shift_total, g_opt, status)
    call assert_int_eq(status, SHAPE_VALID, 'standalone star optimum valid')
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, cstatus)
    call cache_star_convexity_optimum_s(cache, PARAMS7, cz_shift_total, cg_opt, &
            cstatus)
    call cache_free_s(cache)
    call assert_int_eq(cstatus, SHAPE_VALID, 'cold cache star optimum valid')
    call assert_abs_close(z_shift_total, cz_shift_total, 0.0_rk, &
            'star optimum z_shift_total bitwise')
    call assert_abs_close(g_opt, cg_opt, 0.0_rk, 'star optimum g_opt bitwise')

    call test_summary()

contains

    !> Drives all six standalone forms with one parameter vector and asserts
    !! that each reports `expect` and leaves every output it owns zero-filled.
    subroutine reject_all_forms(params, expect, tag)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: expect
        character(len = *), intent(in) :: tag

        real(kind = rk) :: rr(N_THETA), dd(N_THETA)
        real(kind = rk) :: zz(N_POINTS), pp(N_POINTS), ss(N_POINTS)
        real(kind = rk) :: shift, north, south, zn, rn, shift_tot, g
        integer(kind = ik) :: st
        logical :: fnd

        rr = SENTINEL
        call compute_radius_grid_standalone_s(params, thetas, N_POINTS, rr, st)
        call assert_int_eq(st, expect, tag // ': radius grid status')
        call assert_true(all_zero_f(rr), tag // ': radius grid zero-filled')

        rr = SENTINEL
        dd = SENTINEL
        call compute_radius_and_derivative_standalone_s(params, thetas, N_POINTS, &
                rr, dd, st)
        call assert_int_eq(st, expect, tag // ': radius+derivative status')
        call assert_true(all_zero_f(rr) .and. all_zero_f(dd), &
                tag // ': radius+derivative zero-filled')

        shift = SENTINEL
        north = SENTINEL
        south = SENTINEL
        call compute_shape_standalone_s(params, N_POINTS, shift, north, south, st)
        call assert_int_eq(st, expect, tag // ': shape status')
        call assert_true(is_zero_f(shift) .and. is_zero_f(north) &
                .and. is_zero_f(south), tag // ': shape zero-filled')

        zz = SENTINEL
        pp = SENTINEL
        ss = SENTINEL
        shift = SENTINEL
        call compute_rho_z_grid_standalone_s(params, N_POINTS, zz, pp, ss, shift, st)
        call assert_int_eq(st, expect, tag // ': rho(z) grid status')
        call assert_true(all_zero_f(zz) .and. all_zero_f(pp) .and. all_zero_f(ss) &
                .and. is_zero_f(shift), tag // ': rho(z) grid zero-filled')

        zn = SENTINEL
        rn = SENTINEL
        fnd = .true.
        call compute_neck_standalone_s(params, N_POINTS, zn, rn, fnd, st)
        call assert_int_eq(st, expect, tag // ': neck status')
        call assert_true(is_zero_f(zn) .and. is_zero_f(rn) .and. .not. fnd, &
                tag // ': neck zero-filled')

        shift_tot = SENTINEL
        g = SENTINEL
        call compute_star_convexity_optimum_standalone_s(params, N_POINTS, &
                shift_tot, g, st)
        call assert_int_eq(st, expect, tag // ': star optimum status')
        call assert_true(is_zero_f(shift_tot) .and. is_zero_f(g), &
                tag // ': star optimum zero-filled')

    end subroutine reject_all_forms

    !> Element-wise zero-tolerance comparison of two same-length arrays.
    subroutine assert_bits_eq(got, want, label)

        real(kind = rk), intent(in) :: got(:)
        real(kind = rk), intent(in) :: want(:)
        character(len = *), intent(in) :: label

        integer(kind = ik) :: j
        logical :: same

        same = .true.
        do j = 1_ik, size(got, kind = ik)
            if (abs(got(j) - want(j)) > 0.0_rk) same = .false.
        end do

        call assert_true(same, label // ': standalone == cold cache bitwise')

    end subroutine assert_bits_eq

    pure logical function all_zero_f(a) result(ok)

        real(kind = rk), intent(in) :: a(:)

        ok = all(abs(a) <= 0.0_rk)

    end function all_zero_f

    pure logical function is_zero_f(x) result(ok)

        real(kind = rk), intent(in) :: x

        ok = abs(x) <= 0.0_rk

    end function is_zero_f

end program fos_param_standalone_test
