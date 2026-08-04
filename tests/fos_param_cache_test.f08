!> Cache suite: `cache_t` lifecycle, `shape_engine_t` wiring, and the first
!! cached output (`cache_rho_z_grid_s`).
!!
!! Parity anchor is the OLD 1.x grid workflow (`compute_rho_z_grid_s` +
!! `compute_fos_z_shift_f`): the cached cylindrical output must reproduce it
!! node-for-node, including `drho_dz`. The two paths derive u differently
!! (tabled u_i vs z*(1/c)), so the comparison is to a tight tolerance, not
!! bitwise.
!!
!! Minimality is asserted through `cache_recompute_count_f`: a repeat call
!! recomputes nothing, and a pure c-step leaves the f-grid (#3) alone while
!! rebuilding the rho grid (#5).
program fos_param_cache_test

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_init_shared_s, &
            cache_free_s, cache_rho_z_grid_s, cache_recompute_count_f, &
            tables_t, tables_init_s, tables_free_s, &
            rho_z_grid_t, compute_rho_z_grid_s, validate_rho_grid_s, &
            compute_fos_z_shift_f, &
            FOS_ERROR_INVALID_C, FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_BUFFER_MISMATCH
    use fos_parameterization_workers_mod, only: fos_masks_f
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, &
            SHAPE_ERROR_INVALID_INIT, SHAPE_ERROR_TABLES_NOT_INITIALIZED
    use test_utils_mod, only: assert_true, assert_int_eq, assert_close, &
            assert_abs_close, test_summary

    implicit none

    integer(kind = ik), parameter :: N_POINTS = 501_ik
    integer(kind = ik), parameter :: N_PARAMS = 7_ik
    real(kind = rk), parameter :: PARAMS7(N_PARAMS) = &
            [1.5_rk, 0.1_rk, 0.05_rk, 0.02_rk, 0.01_rk, 0.005_rk, 0.002_rk]

    type(cache_t)      :: cache
    type(tables_t), target :: shared_tables
    type(tables_t)     :: dead_tables
    type(rho_z_grid_t) :: grid
    integer(kind = ik) :: status, i, err_code
    integer(kind = ik) :: masks(8)
    real(kind = rk)    :: thetas(4)
    real(kind = rk)    :: params_step(N_PARAMS), params_bad_c(N_PARAMS)
    real(kind = rk)    :: params_rho(N_PARAMS), a3
    real(kind = rk)    :: z(N_POINTS), rho(N_POINTS), drho_dz(N_POINTS)
    real(kind = rk)    :: z_shift, zs_ref
    real(kind = rk)    :: big_z(N_POINTS + 3_ik), big_rho(N_POINTS + 3_ik)
    real(kind = rk)    :: big_drho(N_POINTS + 3_ik)
    logical            :: all_zero
    character(len = 256) :: message

    do i = 1_ik, 4_ik
        thetas(i) = real(i, rk) * PI_C / 5.0_rk
    end do

    !---------------------------------------------------------------------------
    ! Masks: full-width values truncated to the declared parameter count
    !---------------------------------------------------------------------------
    masks = fos_masks_f(8_ik)
    call assert_int_eq(masks(1), 84_ik, 'mask #1 a2 = 84')
    call assert_int_eq(masks(2), 171_ik, 'mask #2 z_shift = 171')
    call assert_int_eq(masks(3), 254_ik, 'mask #3 f_grid = 254')
    call assert_int_eq(masks(4), 254_ik, 'mask #4 beak = 254')
    call assert_int_eq(masks(5), 255_ik, 'mask #5 rho_grid = 255')
    call assert_int_eq(masks(8), 255_ik, 'mask #8 radius_and_derivative = 255')

    masks = fos_masks_f(N_PARAMS)
    call assert_int_eq(masks(2), 43_ik, 'mask #2 truncated to 7 bits')
    call assert_int_eq(masks(5), 127_ik, 'mask #5 truncated to 7 bits')

    !---------------------------------------------------------------------------
    ! cache_init_s: the engine owns the n_params contract
    !---------------------------------------------------------------------------
    call cache_init_s(cache, 0_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'init n_params 0 -> 5')
    call cache_free_s(cache)

    ! 9 parameters must be rejected by the shared cap: the fixed tables carry
    ! k_max = (8+2)/2+1 = 6 orders, exactly enough for 8 parameters and no more.
    call cache_init_s(cache, 9_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'init n_params 9 -> 1')
    call cache_free_s(cache)

    call cache_init_s(cache, 8_ik, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'init n_params 8 -> 0')
    call cache_free_s(cache)

    call cache_init_s(cache, N_PARAMS, 10_ik, thetas, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_GRID, 'init n_points below floor -> 3')
    call cache_free_s(cache)

    !---------------------------------------------------------------------------
    ! cache_init_shared_s: uninitialized tables rejected, shared mode works
    !---------------------------------------------------------------------------
    call cache_init_shared_s(cache, dead_tables, N_PARAMS, status)
    call assert_int_eq(status, SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            'shared init on unbuilt tables -> 6')
    call cache_free_s(cache)

    call tables_init_s(shared_tables, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared tables built')

    call cache_init_shared_s(cache, shared_tables, 0_ik, status)
    call assert_int_eq(status, SHAPE_ERROR_INVALID_INIT, 'shared init n_params 0 -> 5')
    call cache_free_s(cache)

    call cache_init_shared_s(cache, shared_tables, 9_ik, status)
    call assert_int_eq(status, SHAPE_ERROR_TOO_MANY_PARAMS, 'shared init n_params 9 -> 1')
    call cache_free_s(cache)

    call cache_init_shared_s(cache, shared_tables, N_PARAMS, status)
    call assert_int_eq(status, SHAPE_VALID, 'shared init valid')

    ! Freeing a shared-mode cache must leave the caller's tables intact
    call cache_free_s(cache)
    call assert_true(shared_tables%initialized, 'shared tables survive cache_free_s')

    !---------------------------------------------------------------------------
    ! Private-mode cache for the compute cases
    !---------------------------------------------------------------------------
    call cache_init_s(cache, N_PARAMS, N_POINTS, thetas, status)
    call assert_int_eq(status, SHAPE_VALID, 'private init valid')

    !---------------------------------------------------------------------------
    ! Wrong parameter count wins over everything and zero-fills
    !---------------------------------------------------------------------------
    z = 1.0_rk
    rho = 1.0_rk
    drho_dz = 1.0_rk
    z_shift = 1.0_rk
    call cache_rho_z_grid_s(cache, PARAMS7(1:3), z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'short params -> 4')
    call assert_true(all_zero_f(z) .and. all_zero_f(rho) .and. all_zero_f(drho_dz), &
            'short params zero-fills arrays')
    call assert_abs_close(z_shift, 0.0_rk, 0.0_rk, 'short params zero-fills z_shift')

    z = 1.0_rk
    call cache_rho_z_grid_s(cache, PARAMS7(1:0), z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_ERROR_WRONG_PARAM_COUNT, 'empty params -> 4')
    call assert_true(all_zero_f(z), 'empty params zero-fills')

    !---------------------------------------------------------------------------
    ! Happy path: node-for-node parity with the 1.x grid workflow
    !---------------------------------------------------------------------------
    call cache_rho_z_grid_s(cache, PARAMS7, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'rho_z_grid valid')

    call compute_rho_z_grid_s(PARAMS7, N_POINTS, grid, err_code, message)
    call assert_int_eq(err_code, SHAPE_VALID, 'reference grid built')
    zs_ref = compute_fos_z_shift_f(PARAMS7)
    call assert_abs_close(z_shift, zs_ref, 0.0_rk, 'z_shift bitwise-parity')

    all_zero = .true.
    do i = 1_ik, N_POINTS
        call assert_abs_close(z(i), grid%z(i) + zs_ref, 1.0e-14_rk, 'z node parity')
        call assert_abs_close(rho(i), grid%rho(i), 1.0e-14_rk, 'rho node parity')
        ! drho_dz is the u-derivative amplified by 1/sqrt(f) near the tips, so
        ! the ~1 ulp difference in u between the two paths is magnified there:
        ! relative, not absolute, tolerance (observed worst case ~1e-15 relative).
        call assert_close(drho_dz(i), grid%drho_dz(i), 1.0e-13_rk, 'drho_dz node parity')
        if (i > 1_ik .and. i < N_POINTS) then
            if (abs(drho_dz(i)) > 0.0_rk) all_zero = .false.
        end if
    end do
    call assert_true(.not. all_zero, 'drho_dz is not identically zero')

    !---------------------------------------------------------------------------
    ! Minimality: a repeat call recomputes nothing
    !---------------------------------------------------------------------------
    call assert_int_eq(count_f(1_ik), 1_ik, 'a2 computed once')
    call assert_int_eq(count_f(2_ik), 1_ik, 'z_shift computed once')
    call assert_int_eq(count_f(3_ik), 1_ik, 'f_grid computed once')
    call assert_int_eq(count_f(4_ik), 0_ik, 'beak never computed by this path')
    call assert_int_eq(count_f(5_ik), 1_ik, 'rho_grid computed once')

    call cache_rho_z_grid_s(cache, PARAMS7, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'repeat call valid')
    call assert_int_eq(count_f(1_ik), 1_ik, 'repeat: a2 unchanged')
    call assert_int_eq(count_f(2_ik), 1_ik, 'repeat: z_shift unchanged')
    call assert_int_eq(count_f(3_ik), 1_ik, 'repeat: f_grid unchanged')
    call assert_int_eq(count_f(5_ik), 1_ik, 'repeat: rho_grid unchanged')

    !---------------------------------------------------------------------------
    ! c-step: f(u) does not depend on c, the rho grid does
    !---------------------------------------------------------------------------
    params_step = PARAMS7
    params_step(1) = 1.6_rk
    call cache_rho_z_grid_s(cache, params_step, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'c-step valid')
    call assert_int_eq(count_f(1_ik), 1_ik, 'c-step: a2 unchanged')
    call assert_int_eq(count_f(3_ik), 1_ik, 'c-step: f_grid unchanged')
    call assert_int_eq(count_f(2_ik), 2_ik, 'c-step: z_shift recomputed')
    call assert_int_eq(count_f(5_ik), 2_ik, 'c-step: rho_grid recomputed')

    !---------------------------------------------------------------------------
    ! Buffer mismatch: 105, entire actual extent zeroed
    !---------------------------------------------------------------------------
    big_z = 1.0_rk
    big_rho = 1.0_rk
    big_drho = 1.0_rk
    call cache_rho_z_grid_s(cache, PARAMS7, big_z, big_rho, big_drho, z_shift, status)
    call assert_int_eq(status, FOS_ERROR_BUFFER_MISMATCH, 'oversized buffers -> 105')
    call assert_true(all_zero_f(big_z) .and. all_zero_f(big_rho) &
            .and. all_zero_f(big_drho), '105 zero-fills the full extent')

    !---------------------------------------------------------------------------
    ! Degenerate c: 102 before any intermediate runs
    !---------------------------------------------------------------------------
    params_bad_c = 0.0_rk
    params_bad_c(1) = 1.0e-11_rk
    z = 1.0_rk
    call cache_rho_z_grid_s(cache, params_bad_c, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, FOS_ERROR_INVALID_C, 'degenerate c -> 102')
    call assert_true(all_zero_f(z) .and. all_zero_f(rho) .and. all_zero_f(drho_dz), &
            '102 zero-fills')

    !---------------------------------------------------------------------------
    ! rho <= 0 in the interior: 100. Never hardcode an unverified vector —
    ! confirm the rejection against the OLD validate_rho_grid_s path first, and
    ! probe a3 upward if the first guess is accepted.
    !---------------------------------------------------------------------------
    params_rho = 0.0_rk
    params_rho(1) = 1.0_rk
    params_rho(2) = 0.9_rk
    call probe_old_rho_s(params_rho, err_code)
    if (err_code /= FOS_ERROR_RHO_NEGATIVE) then
        do i = 0_ik, 120_ik
            a3 = 0.90_rk + 5.0e-3_rk * real(i, rk)
            params_rho(2) = a3
            call probe_old_rho_s(params_rho, err_code)
            if (err_code == FOS_ERROR_RHO_NEGATIVE) exit
        end do
    end if
    call assert_int_eq(err_code, FOS_ERROR_RHO_NEGATIVE, &
            'probe vector is rho-rejected by the 1.x surface')

    z = 1.0_rk
    call cache_rho_z_grid_s(cache, params_rho, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, FOS_ERROR_RHO_NEGATIVE, 'rho <= 0 -> 100')
    call assert_true(all_zero_f(z) .and. all_zero_f(rho) .and. all_zero_f(drho_dz), &
            '100 zero-fills')

    ! A rejection invalidates: the next good call recomputes every intermediate
    call cache_rho_z_grid_s(cache, PARAMS7, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_VALID, 'recovery call valid')
    call assert_true(count_f(3_ik) > 1_ik, 'rejection invalidated the f-grid')

    !---------------------------------------------------------------------------
    ! Post-free compute: the engine is back to pristine, so 2, not a crash
    !---------------------------------------------------------------------------
    call cache_free_s(cache)
    z = 1.0_rk
    call cache_rho_z_grid_s(cache, PARAMS7, z, rho, drho_dz, z_shift, status)
    call assert_int_eq(status, SHAPE_ERROR_CACHE_NOT_INITIALIZED, &
            'compute after free -> 2')
    call assert_true(all_zero_f(z), 'post-free call zero-fills')
    call assert_int_eq(count_f(5_ik), 0_ik, 'post-free counters reset')

    ! Freeing twice is safe
    call cache_free_s(cache)

    call tables_free_s(shared_tables)
    call test_summary()

contains

    !> Recompute counter as an `ik` integer, for the assertion helpers.
    function count_f(intermediate) result(n)
        integer(kind = ik), intent(in) :: intermediate
        integer(kind = ik) :: n
        n = int(cache_recompute_count_f(cache, intermediate), ik)
    end function count_f

    !> .true. iff every element is exactly zero.
    pure function all_zero_f(a) result(zeroed)
        real(kind = rk), intent(in) :: a(:)
        logical :: zeroed
        integer(kind = ik) :: j
        zeroed = .true.
        do j = 1_ik, size(a, kind = ik)
            if (abs(a(j)) > 0.0_rk) zeroed = .false.
        end do
    end function all_zero_f

    !> Runs the 1.x grid + validation path purely for its rejection code.
    subroutine probe_old_rho_s(p, code)
        real(kind = rk), intent(in) :: p(:)
        integer(kind = ik), intent(out) :: code

        type(rho_z_grid_t) :: probe_grid
        character(len = 256) :: probe_message
        logical :: is_valid

        call compute_rho_z_grid_s(p, N_POINTS, probe_grid, code, probe_message)
        if (code /= 0_ik) return
        call validate_rho_grid_s(probe_grid, p, is_valid, code, probe_message)

    end subroutine probe_old_rho_s

end program fos_param_cache_test
