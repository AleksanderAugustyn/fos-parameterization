!> Analytic dR/dtheta: elemental core, shape split, batch evaluator.
program fos_param_derivative_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: fos_shape_t, make_fos_shape_f, &
            compute_fos_radius_and_derivative_s, compute_radius_fos_with_zshift_s, &
            compute_fos_radius_grid_s, compute_rho_at_z_s, FOS_VALID, &
            compute_fos_shape_s, FOS_ERROR_INVALID_C, &
            compute_fos_radius_and_derivative_at_thetas_s
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    call test_sphere_and_poles()
    call test_wrapper_parity()
    call test_derivative_vs_fd()
    call test_derivative_vs_implicit_formula()
    call test_shape_split()
    call test_batch_evaluator()
    call test_summary()

contains

    !> Sphere (c = 1): R == 1 and dR/dtheta == 0 everywhere; pole branch
    !! returns the exact pole extents with zero derivative.
    subroutine test_sphere_and_poles()
        type(fos_shape_t) :: shape
        real(kind = rk)   :: x(5), r(5), dr(5), params(1)
        integer(kind = ik) :: i

        params = [1.0_rk]
        shape = make_fos_shape_f(params, 0.0_rk)
        x = [-0.9_rk, -0.4_rk, 0.0_rk, 0.5_rk, 0.8_rk]
        call compute_fos_radius_and_derivative_s(shape, x, r, dr)
        do i = 1_ik, 5_ik
            call assert_close(r(i), 1.0_rk, 1.0e-10_rk, 'sphere: R == 1')
            call assert_abs_close(dr(i), 0.0_rk, 1.0e-10_rk, 'sphere: dR == 0')
        end do

        ! Pole branch (|x| > POLE_THRESH): exact extents, zero derivative
        shape = make_fos_shape_f([1.4_rk], 0.05_rk)
        call compute_fos_radius_and_derivative_s(shape, 1.0_rk, r(1), dr(1))
        call assert_close(r(1), 1.45_rk, 1.0e-15_rk, 'pole: r_north = c + z_shift')
        call assert_abs_close(dr(1), 0.0_rk, 1.0e-15_rk, 'pole: dR == 0 north')
        call compute_fos_radius_and_derivative_s(shape, -1.0_rk, r(1), dr(1))
        call assert_close(r(1), 1.35_rk, 1.0e-15_rk, 'pole: r_south = |z_shift - c|')
        call assert_abs_close(dr(1), 0.0_rk, 1.0e-15_rk, 'pole: dR == 0 south')

        ! Degenerate params: unit-sphere fallback with zero derivative
        shape = make_fos_shape_f([0.0_rk], 0.0_rk)
        call compute_fos_radius_and_derivative_s(shape, 0.3_rk, r(1), dr(1))
        call assert_close(r(1), 1.0_rk, 1.0e-15_rk, 'degenerate: sphere fallback')
        call assert_abs_close(dr(1), 0.0_rk, 1.0e-15_rk, 'degenerate: dR == 0')
    end subroutine test_sphere_and_poles

    !> The refactored wrapper must reproduce the elemental core bit-for-bit
    !! (they are the same code path), on an asymmetric valid shape.
    subroutine test_wrapper_parity()
        type(fos_shape_t)    :: shape
        real(kind = rk)      :: params(3), radii_grid(91), z_shift
        real(kind = rk)      :: x, r_wrapper, r_core, dr
        logical              :: is_valid
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        params = [1.5_rk, 0.08_rk, 0.05_rk]
        call compute_fos_radius_grid_s(params, 91_ik, radii_grid, z_shift, &
                is_valid, message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'wrapper: shape is valid')
        call assert_true(is_valid, 'wrapper: is_valid set')

        shape = make_fos_shape_f(params, z_shift)
        do i = 1_ik, 7_ik
            x = -0.9_rk + real(i - 1_ik, rk) * 0.3_rk
            call compute_radius_fos_with_zshift_s(params, x, z_shift, r_wrapper)
            call compute_fos_radius_and_derivative_s(shape, x, r_core, dr)
            call assert_close(r_core, r_wrapper, 1.0e-15_rk, 'wrapper: same R as core')
        end do
    end subroutine test_wrapper_parity

    !> dR/dtheta against 5-point central FD in theta of the core's own R,
    !! tol 1e-9 (spec gate). Newton residual 1e-12 / (12h) stays below tol.
    subroutine test_derivative_vs_fd()
        type(fos_shape_t)    :: shape
        real(kind = rk), parameter :: h = 1.0e-3_rk
        real(kind = rk), parameter :: test_thetas(5) = &
                [0.3_rk, 0.9_rk, PI_C / 2.0_rk, 2.2_rk, 2.9_rk]
        real(kind = rk)      :: params(3), radii_grid(91), z_shift
        real(kind = rk)      :: theta, r, dr, rm2, rm1, rp1, rp2, dr_unused, fd
        logical              :: is_valid
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        params = [1.5_rk, 0.08_rk, 0.05_rk]
        call compute_fos_radius_grid_s(params, 91_ik, radii_grid, z_shift, &
                is_valid, message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'fd: shape is valid')
        shape = make_fos_shape_f(params, z_shift)

        do i = 1_ik, 5_ik
            theta = test_thetas(i)
            call compute_fos_radius_and_derivative_s(shape, cos(theta), r, dr)
            call compute_fos_radius_and_derivative_s(shape, cos(theta - 2.0_rk * h), rm2, dr_unused)
            call compute_fos_radius_and_derivative_s(shape, cos(theta - h), rm1, dr_unused)
            call compute_fos_radius_and_derivative_s(shape, cos(theta + h), rp1, dr_unused)
            call compute_fos_radius_and_derivative_s(shape, cos(theta + 2.0_rk * h), rp2, dr_unused)
            fd = (rm2 - 8.0_rk * rm1 + 8.0_rk * rp1 - rp2) / (12.0_rk * h)
            call assert_abs_close(dr, fd, 1.0e-9_rk, 'fd: analytic dR/dtheta vs 5-point FD')
        end do
    end subroutine test_derivative_vs_fd

    !> Independent re-evaluation of the implicit-differentiation formula at the
    !! returned root (same formula Tier A uses in WMMM's fos_geometry_validation):
    !! dR/dtheta = -(r cos + drho_dz r sin) / (sin - drho_dz cos).
    subroutine test_derivative_vs_implicit_formula()
        type(fos_shape_t)    :: shape
        real(kind = rk), parameter :: test_thetas(5) = &
                [0.3_rk, 0.9_rk, PI_C / 2.0_rk, 2.2_rk, 2.9_rk]
        real(kind = rk)      :: params(3), radii_grid(91), z_shift
        real(kind = rk)      :: theta, x, sin_th, r, dr, rho, drho_dz, expected
        logical              :: is_valid
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        params = [1.5_rk, 0.08_rk, 0.05_rk]
        call compute_fos_radius_grid_s(params, 91_ik, radii_grid, z_shift, &
                is_valid, message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'implicit: shape is valid')
        shape = make_fos_shape_f(params, z_shift)

        do i = 1_ik, 5_ik
            theta = test_thetas(i)
            x = cos(theta)
            sin_th = sqrt(max(1.0_rk - x**2, 0.0_rk))
            call compute_fos_radius_and_derivative_s(shape, x, r, dr)
            call compute_rho_at_z_s(params, r * x, z_shift, rho, drho_dz)
            expected = -(r * x + drho_dz * r * sin_th) / (sin_th - drho_dz * x)
            call assert_close(dr, expected, 1.0e-13_rk, 'implicit: formula parity at root')
        end do
    end subroutine test_derivative_vs_implicit_formula

    !> compute_fos_shape_s must reproduce compute_fos_radius_grid_s's validity,
    !! z_shift, and messages exactly (same code, relocated), and its pole radii
    !! must equal the grid's endpoint values (the grid's theta = 0 and pi entries
    !! come from the Newton routine's pole branch — the same expressions).
    subroutine test_shape_split()
        real(kind = rk)      :: params(3), bad_params(1), radii_grid(91)
        real(kind = rk)      :: z_shift_grid, z_shift, r_north, r_south
        logical              :: is_valid
        integer(kind = ik)   :: code
        character(len = 256) :: message

        params = [1.5_rk, 0.08_rk, 0.05_rk]
        call compute_fos_radius_grid_s(params, 91_ik, radii_grid, z_shift_grid, &
                is_valid, message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'split: reference valid')

        ! compute_fos_radius_grid_s with no n_rho_grid override uses n_grid (91)
        ! as its internal rho-grid size; pass the same for exact parity.
        call compute_fos_shape_s(params, 91_ik, z_shift, r_north, r_south, &
                is_valid, message, code)
        call assert_int_eq(code, FOS_VALID, 'split: shape valid')
        call assert_true(is_valid, 'split: is_valid set')
        call assert_close(z_shift, z_shift_grid, 1.0e-15_rk, 'split: z_shift parity')
        call assert_close(r_north, radii_grid(1), 1.0e-15_rk, 'split: r_north == R(0)')
        call assert_close(r_south, radii_grid(91), 1.0e-15_rk, 'split: r_south == R(pi)')

        bad_params = [-1.0_rk]
        call compute_fos_shape_s(bad_params, 91_ik, z_shift, r_north, r_south, &
                is_valid, message, code)
        call assert_int_eq(code, FOS_ERROR_INVALID_C, 'split: invalid c rejected')
        call assert_true(.not. is_valid, 'split: invalid flag')
        call assert_close(r_north, 0.0_rk, 1.0e-15_rk, 'split: r_north zero-filled')
        call assert_close(r_south, 0.0_rk, 1.0e-15_rk, 'split: r_south zero-filled')
    end subroutine test_shape_split

    !> End-to-end split flow: shape resolve + batch evaluation at the uniform
    !! grid's own thetas must reproduce compute_fos_radius_grid_s bit-for-bit
    !! (identical Newton path underneath).
    subroutine test_batch_evaluator()
        real(kind = rk)      :: params(3), radii_grid(91), z_shift_grid
        real(kind = rk)      :: thetas(91), radii(91), dr(91)
        real(kind = rk)      :: z_shift, r_north, r_south
        logical              :: is_valid
        integer(kind = ik)   :: code, i
        character(len = 256) :: message

        params = [1.5_rk, 0.08_rk, 0.05_rk]
        call compute_fos_radius_grid_s(params, 91_ik, radii_grid, z_shift_grid, &
                is_valid, message, error_code = code)
        call assert_int_eq(code, FOS_VALID, 'batch: reference valid')

        call compute_fos_shape_s(params, 91_ik, z_shift, r_north, r_south, &
                is_valid, message, code)
        do i = 1_ik, 91_ik
            thetas(i) = real(i - 1_ik, rk) * PI_C / 90.0_rk
        end do
        call compute_fos_radius_and_derivative_at_thetas_s(params, thetas, z_shift, radii, dr)

        do i = 1_ik, 91_ik
            call assert_close(radii(i), radii_grid(i), 1.0e-15_rk, 'batch: R parity with grid API')
        end do
        ! Uniform grid includes the poles: derivative must be exactly 0 there
        call assert_abs_close(dr(1), 0.0_rk, 1.0e-15_rk, 'batch: dR == 0 at theta = 0')
        call assert_abs_close(dr(91), 0.0_rk, 1.0e-15_rk, 'batch: dR == 0 at theta = pi')
    end subroutine test_batch_evaluator

end program fos_param_derivative_test
