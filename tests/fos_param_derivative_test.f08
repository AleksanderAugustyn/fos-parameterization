!> Analytic dR/dtheta on the 2.0 tier-1 surface: pole branches, the shape
!! resolve/evaluate split, and the two forms' agreement.
!!
!! 1.x drove the elemental Newton core directly through `fos_shape_t` +
!! `make_fos_shape_f`, both privatized at 2.0. Every probe below now goes
!! through `compute_radius_and_derivative_standalone_s`, which resolves the
!! shift itself — so the cases that 1.x expressed by handing the core an
!! arbitrary z_shift are expressed here against the resolved shift the library
!! reports (`compute_shape_standalone_s`).
program fos_param_derivative_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, compute_shape_standalone_s, &
            compute_rho_at_z_s, FOS_VALID, FOS_ERROR_INVALID_C
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    !> n_points for every resolve here. The tier-1 forms floor the u-grid at
    !! FOS_N_POINTS_FLOOR = 100, so 1.x's 91 is no longer a legal resolution.
    integer(kind = ik), parameter :: N_GRID = 101_ik

    !> Mild asymmetric shape used by every parity/derivative probe.
    real(kind = rk), parameter :: SHAPE3(3) = [1.5_rk, 0.08_rk, 0.05_rk]

    call test_sphere_and_poles()
    call test_form_parity()
    call test_derivative_vs_fd()
    call test_derivative_vs_implicit_formula()
    call test_shape_split()
    call test_batch_evaluator()
    call test_summary()

contains

    !> Sphere (c = 1): R == 1 and dR/dtheta == 0 everywhere, and the pole nodes
    !! take the analytic branch (exact pole extents, zero derivative). A
    !! degenerate vector is a rejection now, not a silent unit-sphere fallback.
    subroutine test_sphere_and_poles()
        real(kind = rk)    :: thetas(5), r(5), dr(5), params(1)
        real(kind = rk)    :: pole_t(2), pole_r(2), pole_dr(2)
        real(kind = rk)    :: z_shift, r_north, r_south
        integer(kind = ik) :: i, status

        params = [1.0_rk]
        thetas = [acos(-0.9_rk), acos(-0.4_rk), acos(0.0_rk), acos(0.5_rk), acos(0.8_rk)]
        call compute_radius_and_derivative_standalone_s(params, thetas, N_GRID, &
                r, dr, status)
        call assert_int_eq(status, FOS_VALID, 'sphere: valid')
        do i = 1_ik, 5_ik
            call assert_close(r(i), 1.0_rk, 1.0e-10_rk, 'sphere: R == 1')
            call assert_abs_close(dr(i), 0.0_rk, 1.0e-10_rk, 'sphere: dR == 0')
        end do

        ! Pole branch: R(0) and R(pi) are the analytic pole extents the resolve
        ! reports, with zero derivative.
        call compute_shape_standalone_s([1.4_rk], N_GRID, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_VALID, 'pole: resolve valid')
        pole_t = [0.0_rk, PI_C]
        call compute_radius_and_derivative_standalone_s([1.4_rk], pole_t, N_GRID, &
                pole_r, pole_dr, status)
        call assert_int_eq(status, FOS_VALID, 'pole: evaluation valid')
        call assert_close(pole_r(1), r_north, 0.0_rk, 'pole: R(0) = r_north')
        call assert_close(pole_r(2), r_south, 0.0_rk, 'pole: R(pi) = r_south')
        call assert_abs_close(pole_dr(1), 0.0_rk, 1.0e-15_rk, 'pole: dR == 0 north')
        call assert_abs_close(pole_dr(2), 0.0_rk, 1.0e-15_rk, 'pole: dR == 0 south')

        ! Degenerate params: 1.x returned r = 1 silently; 2.0 rejects and
        ! zero-fills (spec: the unit-sphere fallback is withdrawn).
        call compute_radius_and_derivative_standalone_s([0.0_rk], pole_t, N_GRID, &
                pole_r, pole_dr, status)
        call assert_int_eq(status, FOS_ERROR_INVALID_C, 'degenerate: rejected with 102')
        call assert_abs_close(pole_r(1), 0.0_rk, 0.0_rk, 'degenerate: R zero-filled')
        call assert_abs_close(pole_dr(1), 0.0_rk, 0.0_rk, 'degenerate: dR zero-filled')
    end subroutine test_sphere_and_poles

    !> The radius-only and radius+derivative forms share one solve, so their
    !! radii must agree bit-for-bit on identical input.
    subroutine test_form_parity()
        real(kind = rk)    :: thetas(N_GRID), radii(N_GRID)
        real(kind = rk)    :: radii_d(N_GRID), dr(N_GRID)
        integer(kind = ik) :: status, i

        call uniform_thetas_s(thetas)
        call compute_radius_grid_standalone_s(SHAPE3, thetas, N_GRID, radii, status)
        call assert_int_eq(status, FOS_VALID, 'parity: radius grid valid')
        call compute_radius_and_derivative_standalone_s(SHAPE3, thetas, N_GRID, &
                radii_d, dr, status)
        call assert_int_eq(status, FOS_VALID, 'parity: radius+derivative valid')

        do i = 1_ik, N_GRID
            call assert_abs_close(radii_d(i), radii(i), 0.0_rk, &
                    'parity: both forms return the same R bit-for-bit')
        end do
    end subroutine test_form_parity

    !> dR/dtheta against 5-point central FD in theta of the form's own R,
    !! tol 1e-9 (spec gate). Newton residual 1e-12 / (12h) stays below tol.
    subroutine test_derivative_vs_fd()
        real(kind = rk), parameter :: h = 1.0e-3_rk
        real(kind = rk), parameter :: test_thetas(5) = &
                [0.3_rk, 0.9_rk, PI_C / 2.0_rk, 2.2_rk, 2.9_rk]
        real(kind = rk)    :: r(5), dr(5), fd
        real(kind = rk)    :: rm2(5), rm1(5), rp1(5), rp2(5), sink(5)
        integer(kind = ik) :: status, i

        call compute_radius_and_derivative_standalone_s(SHAPE3, test_thetas, N_GRID, &
                r, dr, status)
        call assert_int_eq(status, FOS_VALID, 'fd: shape is valid')

        call compute_radius_and_derivative_standalone_s(SHAPE3, &
                test_thetas - 2.0_rk * h, N_GRID, rm2, sink, status)
        call compute_radius_and_derivative_standalone_s(SHAPE3, test_thetas - h, &
                N_GRID, rm1, sink, status)
        call compute_radius_and_derivative_standalone_s(SHAPE3, test_thetas + h, &
                N_GRID, rp1, sink, status)
        call compute_radius_and_derivative_standalone_s(SHAPE3, &
                test_thetas + 2.0_rk * h, N_GRID, rp2, sink, status)

        do i = 1_ik, 5_ik
            fd = (rm2(i) - 8.0_rk * rm1(i) + 8.0_rk * rp1(i) - rp2(i)) / (12.0_rk * h)
            call assert_abs_close(dr(i), fd, 1.0e-9_rk, 'fd: analytic dR/dtheta vs 5-point FD')
        end do
    end subroutine test_derivative_vs_fd

    !> Independent re-evaluation of the implicit-differentiation formula at the
    !! returned root (same formula Tier A uses in WMMM's fos_geometry_validation):
    !! dR/dtheta = -(r cos + drho_dz r sin) / (sin - drho_dz cos).
    subroutine test_derivative_vs_implicit_formula()
        real(kind = rk), parameter :: test_thetas(5) = &
                [0.3_rk, 0.9_rk, PI_C / 2.0_rk, 2.2_rk, 2.9_rk]
        real(kind = rk)    :: r(5), dr(5)
        real(kind = rk)    :: z_shift, r_north, r_south
        real(kind = rk)    :: x, sin_th, rho, drho_dz, expected
        integer(kind = ik) :: status, i

        call compute_shape_standalone_s(SHAPE3, N_GRID, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_VALID, 'implicit: shape is valid')
        call compute_radius_and_derivative_standalone_s(SHAPE3, test_thetas, N_GRID, &
                r, dr, status)
        call assert_int_eq(status, FOS_VALID, 'implicit: evaluation valid')

        do i = 1_ik, 5_ik
            x = cos(test_thetas(i))
            sin_th = sqrt(max(1.0_rk - x**2, 0.0_rk))
            call compute_rho_at_z_s(SHAPE3, r(i) * x, z_shift, rho, drho_dz)
            expected = -(r(i) * x + drho_dz * r(i) * sin_th) / (sin_th - drho_dz * x)
            call assert_close(dr(i), expected, 1.0e-13_rk, 'implicit: formula parity at root')
        end do
    end subroutine test_derivative_vs_implicit_formula

    !> The resolve form's pole radii must equal the R(theta) grid's endpoint
    !! values — the grid's theta = 0 and pi entries come from the Newton
    !! routine's pole branch, which evaluates the same expressions.
    subroutine test_shape_split()
        real(kind = rk)    :: thetas(N_GRID), radii(N_GRID)
        real(kind = rk)    :: z_shift, r_north, r_south
        integer(kind = ik) :: status

        call uniform_thetas_s(thetas)
        call compute_radius_grid_standalone_s(SHAPE3, thetas, N_GRID, radii, status)
        call assert_int_eq(status, FOS_VALID, 'split: reference valid')

        call compute_shape_standalone_s(SHAPE3, N_GRID, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_VALID, 'split: shape valid')
        call assert_close(r_north, radii(1), 0.0_rk, 'split: r_north == R(0)')
        call assert_close(r_south, radii(N_GRID), 0.0_rk, 'split: r_south == R(pi)')
        call assert_close(r_north, SHAPE3(1) + z_shift, 0.0_rk, 'split: r_north = c + z_shift')
        call assert_close(r_south, abs(-SHAPE3(1) + z_shift), 0.0_rk, &
                'split: r_south = |-c + z_shift|')

        call compute_shape_standalone_s([-1.0_rk], N_GRID, z_shift, r_north, r_south, status)
        call assert_int_eq(status, FOS_ERROR_INVALID_C, 'split: invalid c rejected')
        call assert_abs_close(z_shift, 0.0_rk, 0.0_rk, 'split: z_shift zero-filled')
        call assert_abs_close(r_north, 0.0_rk, 0.0_rk, 'split: r_north zero-filled')
        call assert_abs_close(r_south, 0.0_rk, 0.0_rk, 'split: r_south zero-filled')
    end subroutine test_shape_split

    !> Batch evaluation at the uniform grid's own thetas must reproduce the
    !! radius-grid form bit-for-bit (one solve loop underneath), and the
    !! derivative must be exactly 0 at the two poles the grid includes.
    subroutine test_batch_evaluator()
        real(kind = rk)    :: thetas(N_GRID), radii_grid(N_GRID)
        real(kind = rk)    :: radii(N_GRID), dr(N_GRID)
        integer(kind = ik) :: status, i

        call uniform_thetas_s(thetas)
        call compute_radius_grid_standalone_s(SHAPE3, thetas, N_GRID, radii_grid, status)
        call assert_int_eq(status, FOS_VALID, 'batch: reference valid')

        call compute_radius_and_derivative_standalone_s(SHAPE3, thetas, N_GRID, &
                radii, dr, status)
        call assert_int_eq(status, FOS_VALID, 'batch: evaluation valid')

        do i = 1_ik, N_GRID
            call assert_close(radii(i), radii_grid(i), 1.0e-15_rk, &
                    'batch: R parity with the grid form')
        end do
        ! Uniform grid includes the poles: derivative must be exactly 0 there
        call assert_abs_close(dr(1), 0.0_rk, 1.0e-15_rk, 'batch: dR == 0 at theta = 0')
        call assert_abs_close(dr(N_GRID), 0.0_rk, 1.0e-15_rk, 'batch: dR == 0 at theta = pi')
    end subroutine test_batch_evaluator

    !> Uniform theta nodes on [0, pi], endpoints pinned: under -ffast-math the
    !! last product can land one ulp above PI_C, which the library rejects as
    !! an out-of-domain theta.
    pure subroutine uniform_thetas_s(thetas)
        real(kind = rk), intent(out) :: thetas(N_GRID)
        integer(kind = ik) :: i
        do i = 1_ik, N_GRID
            thetas(i) = real(i - 1_ik, rk) * PI_C / real(N_GRID - 1_ik, rk)
        end do
        thetas(1) = 0.0_rk
        thetas(N_GRID) = PI_C
    end subroutine uniform_thetas_s

end program fos_param_derivative_test
