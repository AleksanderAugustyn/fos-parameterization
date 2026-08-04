!> Prints golden values for the regression shapes as paste-ready Python.
!!
!! ## Recapture procedure
!!
!! The goldens in `python/tests/test_api.py` are bit-comparisons against this
!! library. They must be regenerated whenever the arithmetic behind the C API
!! changes — a new evaluator, a new resolve, a new internal resolution — and
!! NOT otherwise: a golden that drifts without a deliberate recapture is a
!! regression, not a refresh.
!!
!!   1. Build the capture program:
!!        `cmake --build build-debug --target fos_param_golden_capture`
!!   2. Run it and keep the output:
!!        `./build-debug/fos_param_golden_capture`
!!      Every line is valid Python. A non-zero exit or a FAIL line means a
!!      regression shape stopped being accepted — investigate, do not paste.
!!   3. Replace the F2/F3/F4/F5 blocks at the top of `python/tests/test_api.py`
!!      with the printed lines, leaving the `*_PARAMS`, `GOLDEN_IDX`,
!!      `DERIV_THETAS` and tolerance definitions alone.
!!   4. Run the Python suite: `pytest python/tests -q`. It must pass with the
!!      tolerances UNCHANGED. Loosening a tolerance to make a recapture fit is
!!      never the fix.
!!
!! Every value is captured through the SAME entry points the Python bindings
!! call, at the same internal resolution — `fos_compute_radius_grid`,
!! `fos_compute_shape` and `fos_compute_radius_and_derivative_at_thetas` from
!! the C API module. Capturing through the Fortran tier-1 forms directly would
!! be a different internal resolution and the goldens would not reproduce.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use c_bindings_mod, only: ik_c, rk_c, c_char
    use fos_parameterization_c_api_mod, only: fos_compute_radius_grid, &
            fos_compute_shape, fos_compute_radius_and_derivative_at_thetas
    use fos_parameterization_mod, only: FOS_VALID
    use test_utils_mod, only: assert_int_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_GRID = 181_ik
    integer(kind = ik), parameter :: IDX(7) = [1_ik, 31_ik, 61_ik, 91_ik, 121_ik, 151_ik, 181_ik]

    ! F2: asymmetric actinide-like saddle; F3: mild symmetric deformation;
    ! F4: deep-necked symmetric shape (same as the smoke test's neck case).
    call capture('F2', [1.80_rk, 0.20_rk, 0.30_rk, 0.01_rk, -0.02_rk, 0.0_rk, 0.0_rk])
    call capture('F3', [1.50_rk, 0.10_rk, 0.20_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    call capture('F4', [2.00_rk, 0.00_rk, 0.50_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    ! F5: marginal star-convex shape recovered by the golden-section shift search
    ! (best origin gives max-T = -0.103 R0, just past the -0.1 margin).
    call capture('F5', [2.00_rk, 0.40_rk, 0.66_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
    call test_summary()

contains

    subroutine capture(name, params)
        character(len = *), intent(in) :: name
        real(kind = rk),    intent(in) :: params(:)

        real(kind = rk_c)      :: c_params(size(params))
        real(kind = rk_c)      :: radii(N_GRID), z_shift
        real(kind = rk_c)      :: shp_z_shift, r_north, r_south
        real(kind = rk_c)      :: eval_thetas(3), eval_r(3), eval_dr(3)
        character(kind = c_char) :: msg_buf(256)
        integer(kind = ik_c)   :: code, n_params
        integer(kind = ik)     :: i

        c_params(:) = real(params(:), rk_c)
        n_params = int(size(params, kind = ik), ik_c)

        code = fos_compute_radius_grid(c_params, n_params, int(N_GRID, ik_c), &
                radii, z_shift, 256_ik_c, msg_buf)
        call assert_int_eq(int(code, ik), FOS_VALID, name // ': valid')
        if (int(code, ik) /= FOS_VALID) return

        write(*, '(A,A)') name, '_EXPECTED = ['
        do i = 1_ik, 7_ik
            write(*, '(A,ES24.16E3,A)') '    ', radii(IDX(i)), ','
        end do
        write(*, '(A)') ']'
        write(*, '(A,A,ES24.16E3)') name, '_Z_SHIFT = ', z_shift

        ! Shape split + derivative goldens (n_rho_grid = N_GRID, used verbatim,
        ! so shp_z_shift must equal the grid API's z_shift above).
        code = fos_compute_shape(c_params, n_params, int(N_GRID, ik_c), &
                shp_z_shift, r_north, r_south, 256_ik_c, msg_buf)
        call assert_int_eq(int(code, ik), FOS_VALID, name // ': shape valid')
        if (int(code, ik) /= FOS_VALID) return

        ! Thetas as pre-rounded literals (pi/8, pi/2, 7pi/8), NOT arithmetic:
        ! Release builds use -ffast-math, which may reassociate 7*PI_C/8 during
        ! constant folding (~1 ulp in theta, ~1e-15 in dR/dtheta). Literals parse
        ! identically here and in Python, so the goldens stay bit-comparable.
        eval_thetas = [0.39269908169872414_rk_c, 1.5707963267948966_rk_c, &
                2.748893571891069_rk_c]
        code = fos_compute_radius_and_derivative_at_thetas(c_params, n_params, &
                eval_thetas, 3_ik_c, shp_z_shift, eval_r, eval_dr)
        call assert_int_eq(int(code, ik), FOS_VALID, name // ': derivative eval valid')
        if (int(code, ik) /= FOS_VALID) return

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
