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
!! ## The stored tolerances must absorb a codegen gap
!!
!! Capture is Debug (step 1), but the PUBLISHED wheel is a different codegen:
!! `ci/build-wheel.sh` sets no `CMAKE_BUILD_TYPE` — so CMakeLists defaults to
!! Release (`-O3 -ffast-math -flto`) — and builds with `GCC_OPTS_MARCH=nehalem`
!! under manylinux2014's devtoolset-10 GCC 10.2. Local Release and a
!! native-march wheel are two further variants. The same goldens are compared
!! against ALL of them, so a tolerance is not a Debug-vs-Debug margin: it must
!! cover the Debug-capture-vs-shipped-library difference in FMA/contraction and
!! constant folding.
!!
!! Measured after the 2.0.0 recapture: worst case 2.2e-15 on `F2_DERIV_DR(1)`
!! (Debug capture vs local Release) against an unchanged `DERIV_DR_ATOL` of
!! 5e-15 — roughly 2.8e-15 of headroom left, and about double what 1.x consumed
!! (~1.2e-15). Spend that budget deliberately: an algorithm change that widens
!! the gap has to be measured against the SHIPPED configuration, not only
!! Debug. If it no longer fits, the honest fix is a per-configuration golden
!! set, never a loosened tolerance. CI does run the Python suite against the
!! repaired wheel before publishing, so a bust surfaces at tag time.
!!
!! Every value is captured through the SAME entry points the Python bindings
!! call, at the same internal resolution — the flat tier-1 C entries
!! `fos_param_radius_grid`, `fos_param_shape` and
!! `fos_param_radius_and_derivative`. Capturing through the Fortran forms
!! directly would bypass the C-level marshalling and the goldens would not be
!! the values Python sees.
program golden_capture

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use c_bindings_mod, only: ik_c, rk_c
    use fos_parameterization_c_api_mod, only: fos_param_radius_grid, &
            fos_param_shape, fos_param_radius_and_derivative
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
        real(kind = rk_c)      :: grid_thetas(N_GRID), radii(N_GRID)
        real(kind = rk_c)      :: shp_z_shift, r_north, r_south
        real(kind = rk_c)      :: eval_thetas(3), eval_r(3), eval_dr(3)
        integer(kind = ik_c)   :: code, n_params
        integer(kind = ik)     :: i

        c_params(:) = real(params(:), rk_c)
        n_params = int(size(params, kind = ik), ik_c)

        ! Uniform nodes with the endpoints pinned: under -ffast-math the product
        ! (n-1)*PI_C/(n-1) can land one ulp ABOVE PI_C, which the theta-domain
        ! check then rejects.
        do i = 1_ik, N_GRID
            grid_thetas(i) = real(i - 1_ik, rk_c) * real(PI_C, rk_c) &
                    / real(N_GRID - 1_ik, rk_c)
        end do
        grid_thetas(1) = 0.0_rk_c
        grid_thetas(N_GRID) = real(PI_C, rk_c)

        call fos_param_radius_grid(c_params, n_params, grid_thetas, &
                int(N_GRID, ik_c), int(N_GRID, ik_c), radii, code)
        call assert_int_eq(int(code, ik), FOS_VALID, name // ': valid')
        if (int(code, ik) /= FOS_VALID) return

        ! Shape split (n_points = N_GRID, the resolution the grid above used),
        ! so shp_z_shift is the total shift baked into those radii.
        call fos_param_shape(c_params, n_params, int(N_GRID, ik_c), &
                shp_z_shift, r_north, r_south, code)
        call assert_int_eq(int(code, ik), FOS_VALID, name // ': shape valid')
        if (int(code, ik) /= FOS_VALID) return

        write(*, '(A,A)') name, '_EXPECTED = ['
        do i = 1_ik, 7_ik
            write(*, '(A,ES24.16E3,A)') '    ', radii(IDX(i)), ','
        end do
        write(*, '(A)') ']'
        write(*, '(A,A,ES24.16E3)') name, '_Z_SHIFT = ', shp_z_shift

        ! Thetas as pre-rounded literals (pi/8, pi/2, 7pi/8), NOT arithmetic:
        ! Release builds use -ffast-math, which may reassociate 7*PI_C/8 during
        ! constant folding (~1 ulp in theta, ~1e-15 in dR/dtheta). Literals parse
        ! identically here and in Python, so the goldens stay bit-comparable.
        eval_thetas = [0.39269908169872414_rk_c, 1.5707963267948966_rk_c, &
                2.748893571891069_rk_c]
        call fos_param_radius_and_derivative(c_params, n_params, eval_thetas, &
                3_ik_c, int(N_GRID, ik_c), eval_r, eval_dr, code)
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
