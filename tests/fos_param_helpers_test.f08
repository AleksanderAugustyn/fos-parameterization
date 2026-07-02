!> Helper unit tests: the a2 volume identity, the intrinsic z-shift as -COM,
!! and coefficient indexing (a_k for k >= 3 lives at params(k - 1)).
program fos_param_helpers_test

    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: compute_fos_f_and_derivatives_s, &
            compute_fos_a2_f, compute_fos_z_shift_f, get_fos_coefficient_f
    use fos_test_reference_mod, only: init_quadrature_s, gl_ref_x, gl_ref_w, N_GL_REF
    use test_utils_mod, only: assert_abs_close, test_summary

    implicit none

    ! Representative parameter sets [c, a3, a4, a5, a6, a7, a8]
    real(kind = rk), parameter :: PARAM_SETS(7, 5) = reshape([ &
            1.00_rk, 0.00_rk, 0.00_rk, 0.00_rk, 0.00_rk, 0.00_rk, 0.00_rk, &
                    1.80_rk, 0.20_rk, 0.30_rk, 0.01_rk, -0.02_rk, 0.00_rk, 0.00_rk, &
                    2.50_rk, 0.40_rk, 0.50_rk, -0.03_rk, 0.04_rk, 0.01_rk, -0.01_rk, &
                    0.90_rk, 0.10_rk, -0.15_rk, 0.05_rk, 0.05_rk, 0.02_rk, 0.01_rk, &
                    3.20_rk, 0.55_rk, 0.70_rk, -0.05_rk, 0.05_rk, -0.02_rk, 0.02_rk], &
            [7, 5])

    real(kind = rk) :: params(7), f_val, integral_f, integral_uf, com, z_sh
    integer(kind = ik) :: i, j
    character(len = 64) :: label

    call init_quadrature_s()

    write(*, '(A)') '=== Helper unit tests ==='

    do j = 1_ik, 5_ik
        params = PARAM_SETS(:, j)

        integral_f = 0.0_rk
        integral_uf = 0.0_rk
        do i = 1_ik, N_GL_REF
            call compute_fos_f_and_derivatives_s(params, gl_ref_x(i), f_val)
            integral_f = integral_f + gl_ref_w(i) * f_val
            integral_uf = integral_uf + gl_ref_w(i) * gl_ref_x(i) * f_val
        end do

        ! a2 volume identity: with a2 from compute_fos_a2_f, int f du = 4/3,
        ! which makes V = pi int rho^2 dz = pi int f du = 4 pi / 3.
        write(label, '(A,I0)') 'a2 identity: param set ', j
        call assert_abs_close(integral_f, 4.0_rk / 3.0_rk, 1.0e-12_rk, trim(label))

        ! z-shift: COM_z of the unshifted shape is
        !   pi int z rho^2 dz / V = (3 c / 4) int u f(u) du
        ! and the intrinsic shift must move the COM to the origin: z_sh = -COM_z.
        com = 0.75_rk * params(1) * integral_uf
        z_sh = compute_fos_z_shift_f(params)
        write(label, '(A,I0)') 'z_shift = -COM: param set ', j
        call assert_abs_close(z_sh, -com, 1.0e-12_rk, trim(label))
    end do

    ! Coefficient indexing: a_k for k >= 3 lives at params(k - 1)
    params = [1.5_rk, 0.2_rk, 0.3_rk, 0.04_rk, 0.05_rk, 0.06_rk, 0.07_rk]
    call assert_abs_close(get_fos_coefficient_f(params, 3_ik), 0.2_rk, 0.0_rk, &
            'indexing: a3 = params(2)')
    call assert_abs_close(get_fos_coefficient_f(params, 4_ik), 0.3_rk, 0.0_rk, &
            'indexing: a4 = params(3)')
    call assert_abs_close(get_fos_coefficient_f(params, 8_ik), 0.07_rk, 0.0_rk, &
            'indexing: a8 = params(7)')
    call assert_abs_close(get_fos_coefficient_f(params, 2_ik), compute_fos_a2_f(params), &
            0.0_rk, 'indexing: a2 from volume constraint')
    call assert_abs_close(get_fos_coefficient_f(params, 1_ik), 0.0_rk, 0.0_rk, &
            'indexing: a1 = 0')
    call assert_abs_close(get_fos_coefficient_f(params, 9_ik), 0.0_rk, 0.0_rk, &
            'indexing: a9 beyond array = 0')

    call test_summary()

end program fos_param_helpers_test
