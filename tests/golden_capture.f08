!> Prints golden values for the regression shapes as paste-ready Python.
!! Output goes into python/tests/test_api.py.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: compute_fos_radius_grid_s, FOS_VALID, &
            compute_fos_shape_s, compute_fos_radius_and_derivative_at_thetas_s
    use test_utils_mod, only: assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID = 181_ik
    integer(kind = ik), parameter :: IDX(7) = [1_ik, 31_ik, 61_ik, 91_ik, 121_ik, 151_ik, 181_ik]

    ! F2: asymmetric actinide-like saddle; F3: mild symmetric deformation;
    ! F4: deep-necked symmetric shape (same as the smoke test's neck case).
    call capture('F2', [1.80_rk, 0.20_rk, 0.30_rk, 0.01_rk, -0.02_rk, 0.0_rk, 0.0_rk])
    call capture('F3', [1.50_rk, 0.10_rk, 0.20_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    call capture('F4', [2.00_rk, 0.00_rk, 0.50_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    call test_summary()

contains

    subroutine capture(name, params)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:)

        real(kind = rk)      :: radii(N_GRID), z_shift
        real(kind = rk)      :: shp_z_shift, r_north, r_south
        real(kind = rk)      :: eval_thetas(3), eval_r(3), eval_dr(3)
        logical              :: is_valid
        integer(kind = ik)   :: code, i
        character(len = 256) :: msg

        call compute_fos_radius_grid_s(params, N_GRID, radii, z_shift, is_valid, &
                msg, error_code = code)
        call assert_int_eq(code, FOS_VALID, name // ': valid')
        if (code /= FOS_VALID) return

        write(*, '(A,A)') name, '_EXPECTED = ['
        do i = 1_ik, 7_ik
            write(*, '(A,ES24.16E3,A)') '    ', radii(IDX(i)), ','
        end do
        write(*, '(A)') ']'
        write(*, '(A,A,ES24.16E3)') name, '_Z_SHIFT = ', z_shift

        ! Shape split + derivative goldens (n_rho_grid = N_GRID, used verbatim,
        ! so shp_z_shift must equal the grid API's z_shift above).
        call compute_fos_shape_s(params, N_GRID, shp_z_shift, r_north, r_south, &
                is_valid, msg, code)
        call assert_int_eq(code, FOS_VALID, name // ': shape valid')
        if (code /= FOS_VALID) return

        ! Thetas as pre-rounded literals (pi/8, pi/2, 7pi/8), NOT arithmetic:
        ! Release builds use -ffast-math, which may reassociate 7*PI_C/8 during
        ! constant folding (~1 ulp in theta, ~1e-15 in dR/dtheta). Literals parse
        ! identically here and in Python, so the goldens stay bit-comparable.
        eval_thetas = [0.39269908169872414_rk, 1.5707963267948966_rk, &
                2.748893571891069_rk]
        call compute_fos_radius_and_derivative_at_thetas_s(params, eval_thetas, &
                shp_z_shift, eval_r, eval_dr)

        write(*, '(A,A,ES24.16E3)') name, '_R_NORTH = ', r_north
        write(*, '(A,A,ES24.16E3)') name, '_R_SOUTH = ', r_south
        write(*, '(A,A)') name, '_DERIV_R = ['
        do i = 1_ik, 3_ik
            write(*, '(A,ES24.16E3,A)') '    ', eval_r(i), ','
        end do
        write(*, '(A)') ']'
        write(*, '(A,A)') name, '_DERIV_DR = ['
        do i = 1_ik, 3_ik
            write(*, '(A,ES24.16E3,A)') '    ', eval_dr(i), ','
        end do
        write(*, '(A)') ']'
    end subroutine capture

end program golden_capture
