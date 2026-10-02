!> Worker-kernel suite: status-reporting a2/z_shift, tabled f-grid and beak
!! scan, rho scaling, and the bounded-bracket Newton radius core.
!!
!! Parity anchors survive the 1.x removal in three forms: the raw evaluators
!! that stayed public (`get_fos_coefficient_f`, whose k = 2 branch IS the
!! volume-constraint a2, and `compute_fos_f_and_derivatives_s`), and — where
!! the 1.x routine was deleted outright — a literal frozen from that surface
!! before it went (`ZS_PARAMS7_1X`, `BEAK_A4_1X`). The kernels must reproduce
!! them node-for-node, bitwise where the arithmetic is identical.
program fos_param_workers_test

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: compute_a2_s, compute_z_shift_s, &
            get_fos_coefficient_f, compute_fos_f_and_derivatives_s, &
            FOS_ERROR_INVALID_C
    use fos_parameterization_workers_mod, only: tables_t, &
            tables_free_s, build_tables_s, fos_bundle_t, compute_f_grid_s, &
            beak_scan_f_min_s, scale_rho_grid_s, newton_radius_s, &
            active_length_f, table_order_f, shifted_origin_f, refine_neck_s
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_WRONG_PARAM_COUNT
    use test_utils_mod, only: assert_true, assert_int_eq, assert_abs_close, &
            assert_bits_eq, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    real(kind = rk), parameter :: PARAMS7(7) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    !> Frozen 1.x anchors, captured from `compute_fos_z_shift_f` and the 1.x
    !! beak probe at commit 648428c, the last commit that carried that surface.
    !! `compute_fos_z_shift_f` has no 2.0 counterpart to compare against (the
    !! status form `compute_z_shift_s` IS the code under test), so the number
    !! itself is the anchor.
    real(kind = rk), parameter :: ZS_PARAMS7_1X = 6.565141402540683E-002_rk

    !> a4 of the symmetric (c = 2) family used as the below-threshold beak
    !! fixture: f(0) = 1 - 4 a4 / 3 = 2.0e-4 < F_MIN_THRESHOLD = 5.0e-4.
    !! Historical note: the 1.x-frozen probe value was 0.7495 (f(0) =
    !! 6.67e-4), a beak only under the 1.x-era threshold 1e-3; the 2026-08-11
    !! retune (see F_MIN_THRESHOLD) made that shape representable, so the
    !! fixture deepened.
    real(kind = rk), parameter :: BEAK_A4_1X = 0.74985_rk

    !> A short vector, and a necked symmetric shape whose neck is analytic:
    !! rho_neck = sqrt((1 - 4 a4 / 3) / c) at z = 0.
    real(kind = rk), parameter :: SHORT3(3) = [1.5_rk, 0.1_rk, 0.05_rk]
    real(kind = rk), parameter :: NECK3(3) = [2.0_rk, 0.0_rk, 0.4_rk]

    !> Bit pattern of -0.0. A `-0.0_rk` literal may be folded to +0.0 in
    !! Release, so the value is built from its bits into a volatile variable.
    integer(kind = ikl), parameter :: NEG_ZERO_BITS = int(z'8000000000000000', kind = ikl)

    type(tables_t)     :: tables, tables_small, tables_big
    type(fos_bundle_t) :: bundle
    integer(kind = ik) :: status, i
    real(kind = rk)    :: thetas(4), params51(51), params_beak(3)
    real(kind = rk)    :: a2_new, zs, f_ref, fp_ref, f_min
    real(kind = rk)    :: r, dr, rho_max, c_oblate
    real(kind = rk)    :: f_grid(N_POINTS), fp_grid(N_POINTS)
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
    real(kind = rk)    :: f_big(N_POINTS), fp_big(N_POINTS)
    real(kind = rk)    :: padded(6), a2_pad, zs_pad, f_min_big, z_neck, rho_neck
    real(kind = rk), volatile :: neg_zero
    logical            :: beak_ok, rho_positive, converged, found
    ! NOTE (-Werror hygiene): declare ONLY what this program uses.

    do i = 1_ik, 4_ik
        thetas(i) = real(i, rk) * PI_C / 5.0_rk
    end do

    ! Order 6: what a max_params = 8 cache builds.
    call build_tables_s(tables, N_POINTS, thetas, table_order_f(8_ik), status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'tables built')

    !---------------------------------------------------------------------------
    ! a2: parity with the raw coefficient reader, empty ok, oversize rejected
    !---------------------------------------------------------------------------
    call compute_a2_s(PARAMS7, a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'a2 valid')
    ! get_fos_coefficient_f(p, 2) is the surviving public spelling of the 1.x
    ! volume-constraint a2 — same code, so the comparison is exact.
    call assert_abs_close(a2_new, get_fos_coefficient_f(PARAMS7, 2_ik), 0.0_rk, &
            'a2 bitwise-parity')

    call compute_a2_s(PARAMS7(1:0), a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'a2 empty ok')
    call assert_abs_close(a2_new, 0.0_rk, 0.0_rk, 'a2 empty = 0')

    do i = 1_ik, 51_ik
        params51(i) = 0.01_rk * real(i, rk)
    end do
    call compute_a2_s(params51, a2_new, status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_TOO_MANY_PARAMS), 'a2 51 params rejected')

    !---------------------------------------------------------------------------
    ! z_shift: parity, empty -> 4, degenerate c -> 102
    !---------------------------------------------------------------------------
    call compute_z_shift_s(PARAMS7, zs, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'z_shift valid')
    call assert_abs_close(zs, ZS_PARAMS7_1X, 0.0_rk, 'z_shift bitwise-parity with 1.x')

    call compute_z_shift_s(PARAMS7(1:0), zs, status)
    call assert_int_eq(int(status), int(SHAPE_ERROR_WRONG_PARAM_COUNT), 'z_shift empty -> 4')
    call compute_z_shift_s([1.0e-11_rk, 0.1_rk], zs, status)
    call assert_int_eq(int(status), int(FOS_ERROR_INVALID_C), 'z_shift degenerate c -> 102')

    !---------------------------------------------------------------------------
    ! f-grid: parity against the live evaluator at every interior node
    !---------------------------------------------------------------------------
    call compute_f_grid_s(tables, PARAMS7, f_grid, fp_grid)
    do i = 2_ik, tables%n_points - 1_ik
        call compute_fos_f_and_derivatives_s(PARAMS7, tables%u(i), f_ref, fp_ref)
        call assert_abs_close(f_grid(i), f_ref, 1.0e-14_rk, 'f_grid node parity')
        call assert_abs_close(fp_grid(i), fp_ref, 1.0e-13_rk, 'fp_grid node parity')
    end do

    !---------------------------------------------------------------------------
    ! Beak scan: sphere passes; a vector the 1.x surface rejects as a beak fails
    !---------------------------------------------------------------------------
    ! The vector is not a guess: BEAK_A4_1X is the below-threshold beak
    ! fixture (see the declaration for its 1.x provenance and 2.1.0 retune).
    params_beak(1) = 2.0_rk
    params_beak(2) = 0.0_rk
    params_beak(3) = BEAK_A4_1X

    call beak_scan_f_min_s(tables, [1.0_rk], f_min, beak_ok)
    call assert_true(beak_ok, 'sphere passes beak scan')
    call assert_true(f_min > 0.0_rk, 'sphere f_min positive')

    call beak_scan_f_min_s(tables, params_beak, f_min, beak_ok)
    call assert_true(.not. beak_ok, 'beak vector fails scan')

    !---------------------------------------------------------------------------
    ! rho scaling: sphere geometry, tip convention, shift, rho-positivity verdict
    !---------------------------------------------------------------------------
    call compute_f_grid_s(tables, [1.0_rk], f_grid, fp_grid)
    call scale_rho_grid_s(tables, 1.0_rk, 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_true(rho_positive, 'sphere rho grid positive')
    call assert_abs_close(rho_max, 1.0_rk, 1.0e-12_rk, 'sphere rho_max = 1')
    call assert_abs_close(rho(1), 0.0_rk, 0.0_rk, 'south tip rho = 0')
    call assert_abs_close(drho_dz(N_POINTS), 0.0_rk, 0.0_rk, 'north tip drho/dz = 0')
    i = (N_POINTS + 1_ik) / 2_ik
    call assert_abs_close(z(i), 0.0_rk, 1.0e-15_rk, 'equator z = 0')
    call assert_abs_close(rho(i), 1.0_rk, 1.0e-15_rk, 'equator rho = 1')

    ! z(i) = c*u(i) + z_shift_intrinsic, rho scales as 1/sqrt(c)
    call scale_rho_grid_s(tables, 2.0_rk, 0.5_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_abs_close(z(1), -1.5_rk, 1.0e-14_rk, 'z(1) = -c + z_shift')
    call assert_abs_close(z(N_POINTS), 2.5_rk, 1.0e-14_rk, 'z(n) = c + z_shift')
    call assert_abs_close(rho_max, 1.0_rk / sqrt(2.0_rk), 1.0e-12_rk, 'rho_max = 1/sqrt(c)')

    ! f(0) = 0 exactly for a4 = 0.75: the interior node at u = 0 pinches shut
    call compute_f_grid_s(tables, [1.0_rk, 0.0_rk, 0.75_rk], f_grid, fp_grid)
    call scale_rho_grid_s(tables, 1.0_rk, 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call assert_true(.not. rho_positive, 'pinched shape flagged rho <= 0')

    !---------------------------------------------------------------------------
    ! Newton with the analytic bracket: extreme-oblate correctness
    !---------------------------------------------------------------------------
    ! c = 2e-10, single param: rho(z = 0) = sqrt(f(0)/c) = 1/sqrt(c) ~ 7.07e4.
    ! The 1.x doubling bracket silently returns ~1e-7 here.
    c_oblate = 2.0e-10_rk
    bundle%n_params = 1_ik
    bundle%params(1) = c_oblate
    bundle%z_shift = 0.0_rk
    bundle%r_hi_bound = 2.0_rk * sqrt((1.0_rk / sqrt(c_oblate))**2 + c_oblate**2)
    call newton_radius_s(bundle, 0.0_rk, r, dr, converged)
    call assert_true(converged, 'extreme oblate converges')
    call assert_abs_close(r, 1.0_rk / sqrt(c_oblate), 1.0e-6_rk * 7.1e4_rk, &
            'equatorial radius ~ 1/sqrt(c)')

    ! Sphere sanity: r = 1 everywhere, converged
    bundle%n_params = 1_ik
    bundle%params(1) = 1.0_rk
    bundle%z_shift = 0.0_rk
    bundle%r_hi_bound = 2.0_rk * sqrt(2.0_rk)
    call newton_radius_s(bundle, cos(1.1_rk), r, dr, converged)
    call assert_true(converged, 'sphere converges')
    call assert_abs_close(r, 1.0_rk, 1.0e-10_rk, 'sphere r = 1')
    call assert_abs_close(dr, 0.0_rk, 1.0e-8_rk, 'sphere dR/dtheta = 0')

    ! Poles are analytic and always converged
    call newton_radius_s(bundle, 1.0_rk, r, dr, converged)
    call assert_true(converged, 'north pole converged')
    call assert_abs_close(r, 1.0_rk, 0.0_rk, 'north pole r = c + z_shift')
    call newton_radius_s(bundle, -1.0_rk, r, dr, converged)
    call assert_true(converged, 'south pole converged')
    call assert_abs_close(r, 1.0_rk, 0.0_rk, 'south pole r = |-c + z_shift|')

    !---------------------------------------------------------------------------
    ! Active length: trailing zeros of either sign are dropped, interior kept
    !---------------------------------------------------------------------------
    neg_zero = transfer(NEG_ZERO_BITS, 1.0_rk)

    call assert_int_eq(active_length_f(PARAMS7), 7_ik, 'active length: full vector')
    do i = 1_ik, 6_ik
        padded(i) = 0.0_rk
    end do
    do i = 1_ik, 3_ik
        padded(i) = SHORT3(i)
    end do
    call assert_int_eq(active_length_f(padded), 3_ik, 'active length: +0 padding dropped')
    padded(5) = neg_zero
    padded(6) = neg_zero
    call assert_true(transfer(padded(6), 0_ikl) == NEG_ZERO_BITS, &
            'the -0.0 padding reached memory')
    call assert_int_eq(active_length_f(padded), 3_ik, 'active length: -0 padding dropped')

    !---------------------------------------------------------------------------
    ! a2 and z_shift trim: short == zero-padded, bit for bit
    !---------------------------------------------------------------------------
    call compute_a2_s(SHORT3, a2_new, status)
    call compute_a2_s(padded, a2_pad, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'a2 padded valid')
    call assert_bits_eq(a2_pad, a2_new, 'a2: short == zero-padded')
    call compute_z_shift_s(SHORT3, zs, status)
    call compute_z_shift_s(padded, zs_pad, status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'z_shift padded valid')
    call assert_bits_eq(zs_pad, zs, 'z_shift: short == zero-padded')

    padded(2) = 0.0_rk
    call assert_int_eq(active_length_f(padded), 3_ik, 'active length: interior zero kept')
    padded(1) = 0.0_rk
    padded(3) = 0.0_rk
    call assert_int_eq(active_length_f(padded), 0_ik, 'active length: all-zero -> 0')
    call assert_int_eq(active_length_f(PARAMS7(1:0)), 0_ik, 'active length: empty -> 0')

    !---------------------------------------------------------------------------
    ! Table order: (n + 2)/2 + 1
    !---------------------------------------------------------------------------
    call assert_int_eq(table_order_f(1_ik), 2_ik, 'table order: 1 param -> 2')
    call assert_int_eq(table_order_f(3_ik), 3_ik, 'table order: 3 params -> 3')
    call assert_int_eq(table_order_f(8_ik), 6_ik, 'table order: 8 params -> 6')
    call assert_int_eq(table_order_f(50_ik), 27_ik, 'table order: 50 params -> 27')

    !---------------------------------------------------------------------------
    ! Kernel loop bound: tables of different order give the same bits, because
    ! the kernels sum the orders the VECTOR needs
    !---------------------------------------------------------------------------
    call build_tables_s(tables_small, N_POINTS, thetas, table_order_f(3_ik), status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'order-3 tables built')
    call build_tables_s(tables_big, N_POINTS, thetas, table_order_f(50_ik), status)
    call assert_int_eq(int(status), int(SHAPE_VALID), 'order-27 tables built')

    call compute_f_grid_s(tables_small, SHORT3, f_grid, fp_grid)
    call compute_f_grid_s(tables_big, SHORT3, f_big, fp_big)
    do i = 1_ik, N_POINTS
        call assert_bits_eq(f_big(i), f_grid(i), 'f grid: order-27 tables == order-3 tables')
        call assert_bits_eq(fp_big(i), fp_grid(i), 'fp grid: order-27 tables == order-3 tables')
    end do

    call beak_scan_f_min_s(tables_small, SHORT3, f_min, beak_ok)
    call beak_scan_f_min_s(tables_big, SHORT3, f_min_big, beak_ok)
    call assert_bits_eq(f_min_big, f_min, 'f_min: order-27 tables == order-3 tables')

    call tables_free_s(tables_small)
    call tables_free_s(tables_big)

    !---------------------------------------------------------------------------
    ! Neck kernel: analytic neck of the symmetric family; the sphere has none
    !---------------------------------------------------------------------------
    call compute_f_grid_s(tables, NECK3, f_grid, fp_grid)
    call scale_rho_grid_s(tables, NECK3(1), 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call refine_neck_s(tables, NECK3, rho, 0.0_rk, z_neck, rho_neck, found)
    call assert_true(found, 'neck kernel: necked shape has a neck')
    call assert_abs_close(z_neck, 0.0_rk, 1.0e-9_rk, 'neck kernel: z_neck = 0')
    call assert_abs_close(rho_neck, &
            sqrt((1.0_rk - 4.0_rk * NECK3(3) / 3.0_rk) / NECK3(1)), 1.0e-9_rk, &
            'neck kernel: rho_neck analytic')

    call compute_f_grid_s(tables, [1.0_rk], f_grid, fp_grid)
    call scale_rho_grid_s(tables, 1.0_rk, 0.0_rk, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)
    call refine_neck_s(tables, [1.0_rk], rho, 0.0_rk, z_neck, rho_neck, found)
    call assert_true(.not. found, 'neck kernel: the sphere has no neck')
    call assert_abs_close(rho_neck, 0.0_rk, 0.0_rk, 'neck kernel: no neck -> rho_neck = 0')

    !---------------------------------------------------------------------------
    ! Shifted origin: one exact addition
    !---------------------------------------------------------------------------
    call assert_bits_eq(shifted_origin_f(0.25_rk, 0.5_rk), 0.75_rk, &
            'shifted origin: 0.25 + 0.5 = 0.75')

    call tables_free_s(tables)
    call test_summary()

end program fos_param_workers_test
