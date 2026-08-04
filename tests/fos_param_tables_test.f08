program fos_param_tables_test
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: tables_t, tables_init_s, tables_free_s, &
            compute_fos_f_and_derivatives_s, FOS_VALID
    use shape_core_mod, only: SHAPE_ERROR_INVALID_GRID
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, test_summary
    implicit none

    type(tables_t) :: tables
    integer(kind = ik) :: status, i
    real(kind = rk) :: thetas(64), f_ref, fp_ref, f_tab, fp_tab
    real(kind = rk), parameter :: params(7) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]
    real(kind = rk) :: sum_f, sum_fp
    ! NOTE (-Werror hygiene): declare ONLY what this program uses; the
    ! contained helper declares its own loop index k.

    do i = 1, 64
        thetas(i) = real(i, rk) * PI_C / 65.0_rk
    end do

    ! Rejections, validation order: n_points floor, theta count, theta domain
    call tables_init_s(tables, 99_ik, thetas, status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_INVALID_GRID), 'n_points 99 rejected')
    call tables_init_s(tables, 100_ik, thetas(1:0), status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_INVALID_GRID), 'empty thetas rejected')
    call tables_init_s(tables, 100_ik, [PI_C + 0.1_rk], status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_INVALID_GRID), 'theta > pi rejected')
    call tables_init_s(tables, 100_ik, [-0.1_rk], status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_INVALID_GRID), 'theta < 0 rejected')

    call tables_init_s(tables, 501_ik, thetas, status)
    call assert_int_eq(int(status), int(FOS_VALID), 'valid init accepted')
    call assert_true(tables%initialized, 'initialized flag set')
    call assert_abs_close(tables%u(1), -1.0_rk, 1.0e-15_rk, 'u(1) = -1')
    call assert_abs_close(tables%u(501), 1.0_rk, 1.0e-15_rk, 'u(n) = +1')

    ! Trig-basis correctness: f from tables must match the live evaluator
    ! at an interior node. Reproduce the table-based sum exactly as the
    ! worker will compute it (a2 folded in as coefficient index 2).
    do i = 100, 400, 150
        call compute_fos_f_and_derivatives_s(params, tables%u(i), f_ref, fp_ref)
        sum_f = 0.0_rk; sum_fp = 0.0_rk
        ! (test-local replica of the tabled sum; see step 3 worker for the real one)
        call table_f_sum(tables, params, i, sum_f, sum_fp)
        f_tab = 1.0_rk - tables%u(i)**2 - sum_f
        fp_tab = -2.0_rk * tables%u(i) - sum_fp
        call assert_abs_close(f_tab, f_ref, 1.0e-14_rk, 'tabled f matches live f')
        call assert_abs_close(fp_tab, fp_ref, 1.0e-13_rk, 'tabled fp matches live fp')
    end do

    call tables_free_s(tables)
    call assert_true(.not. tables%initialized, 'freed')
    call test_summary()
contains
    subroutine table_f_sum(t, p, i, sf, sfp)
        use fos_parameterization_mod, only: get_fos_coefficient_f
        type(tables_t), intent(in) :: t
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(in) :: i
        real(kind = rk), intent(out) :: sf, sfp
        integer(kind = ik) :: k
        real(kind = rk) :: ae, ao
        sf = 0.0_rk; sfp = 0.0_rk
        do k = 1_ik, t%k_max
            ae = get_fos_coefficient_f(p, 2_ik * k)
            ao = get_fos_coefficient_f(p, 2_ik * k + 1_ik)
            sf = sf + ae * t%cos_even(i, k) + ao * t%sin_odd(i, k)
            sfp = sfp - ae * t%omega(k) * t%sin_even(i, k) + ao * t%psi(k) * t%cos_odd(i, k)
        end do
    end subroutine table_f_sum
end program fos_param_tables_test
