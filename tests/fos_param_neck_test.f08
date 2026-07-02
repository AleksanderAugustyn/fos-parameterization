!> Validates neck finding on the symmetric (c, a4) family where the neck
!! has a closed form: a3 = a5 = a6 = 0 gives a2 = a4/3 (volume constraint),
!! the shape is reflection-symmetric, and the neck sits at u = 0 with
!!   rho_neck = sqrt(f(0)/c) = sqrt((1 - 4 a4 / 3) / c).
!! A neck (local minimum of f at u = 0) exists when
!!   f''(0) = -2 + a4 (pi^2/12 + 9 pi^2/4) > 0, i.e. a4 > 6/(7 pi^2) ~ 0.0868.
!! Note f(0) -> 0 at a4 = 0.75: the upper a4 bound is the scission line
!! of this family.
program fos_param_neck_test

    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: compute_fos_neck_s
    use test_utils_mod, only: assert_true, assert_abs_close, test_summary

    implicit none

    real(kind = rk), parameter :: TOL_NECK = 1.0e-9_rk
    real(kind = rk), parameter :: NECK_C(5) = [1.0_rk, 1.5_rk, 2.0_rk, 2.5_rk, 3.0_rk]

    real(kind = rk) :: params(7), a4, z_neck, rho_neck, rho_exact
    logical :: found
    integer(kind = ik) :: i, j
    character(len = 64) :: label

    write(*, '(A)') '=== Neck validation (symmetric c x a4 family) ==='

    do i = 1_ik, 5_ik
        do j = 0_ik, 12_ik
            a4 = 0.10_rk + 0.05_rk * real(j, rk)
            params = 0.0_rk
            params(1) = NECK_C(i)
            params(3) = a4
            write(label, '(A,F4.2,A,F4.2)') 'neck c=', NECK_C(i), ' a4=', a4

            call compute_fos_neck_s(params, z_neck, rho_neck, found)
            call assert_true(found, trim(label) // ': found')
            if (.not. found) cycle

            rho_exact = sqrt((1.0_rk - 4.0_rk * a4 / 3.0_rk) / NECK_C(i))
            call assert_abs_close(z_neck, 0.0_rk, TOL_NECK, trim(label) // ': z_neck = 0')
            call assert_abs_close(rho_neck, rho_exact, TOL_NECK, &
                    trim(label) // ': rho_neck analytic')
        end do
    end do

    call test_summary()

end program fos_param_neck_test
