!> C-interop layer for the FoS parameterization library (`fos_parameterization.h`).
!!
!! Two opaque handles, each a `c_loc` of a heap-allocated derived type:
!!   - `fos_param_tables_create` -> `tables_t` (shared, immutable after create)
!!   - `fos_param_cache_create*` -> `cache_t`  (THREAD-CONFINED, mutated by
!!                                              every compute call)
!!
!! Every `_create` returns a null handle on failure. The cached computes return
!! the status code directly; the flat tier-1 entries report through a trailing
!! nullable `int* status` (an absent `optional` dummy when C passes NULL). No
!! entry point stops, allocates for the caller, or formats a message:
!! diagnostics are the static strings behind `fos_param_status_message`.
!!
!! This layer only marshals. The buffer-length and parameter-count contracts
!! are checked once, in the Fortran tiers, and their codes pass through
!! untouched: an output buffer whose stated C size differs from the handle's
!! own extent becomes a Fortran array of that size, so the tier's own size
!! comparison rejects it with `FOS_ERROR_BUFFER_MISMATCH` (105) and zero-fills
!! exactly the caller's stated extent. The checks that live HERE are the ones
!! the Fortran layer cannot see: the NULL-handle guard and the negative-size
!! guard on the C integer arguments.
module fos_parameterization_c_api_mod

    use, intrinsic :: iso_c_binding, only: &
            c_int, c_double, c_char, c_ptr, c_loc, c_f_pointer, c_associated, &
            c_null_ptr, c_null_char
    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: &
            tables_t, tables_init_s, tables_free_s, &
            cache_t, cache_init_s, cache_init_shared_s, cache_free_s, &
            cache_shape_s, cache_rho_z_grid_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_radius_and_derivative_at_thetas_s, &
            cache_neck_s, cache_star_convexity_optimum_s, &
            compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, compute_neck_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            compute_a2_s, compute_z_shift_s, compute_rho_at_z_s, &
            STATUS_MESSAGE_LEN, &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE, FOS_ERROR_BUFFER_MISMATCH

    implicit none

    private

    public :: fos_param_status_message
    public :: fos_param_tables_create, fos_param_tables_destroy
    public :: fos_param_cache_create, fos_param_cache_create_shared
    public :: fos_param_cache_destroy
    public :: fos_param_cache_radius_grid
    public :: fos_param_cache_radius_and_derivative
    public :: fos_param_cache_radius_and_derivative_at_thetas
    public :: fos_param_cache_shape
    public :: fos_param_cache_rho_z_grid
    public :: fos_param_cache_neck
    public :: fos_param_cache_star_convexity_optimum
    public :: fos_param_radius_grid
    public :: fos_param_radius_and_derivative
    public :: fos_param_shape
    public :: fos_param_rho_z_grid
    public :: fos_param_neck
    public :: fos_param_star_convexity_optimum
    public :: fos_param_z_shift
    public :: fos_param_a2
    public :: fos_param_rho_at_z

    !---------------------------------------------------------------------------
    ! Static status strings
    !---------------------------------------------------------------------------
    !> One null-terminated string per status code, laid out as fixed columns of a
    !! module array that is initialized at compile time and never written again.
    !! `fos_param_status_message` hands out `c_loc` of a column, so the pointer
    !! is valid forever and every thread gets the same read-only bytes — a shared
    !! mutable buffer would make the message racy.
    !!
    !! The texts are the `status_message` texts: keep the two in sync (a C caller
    !! and a Fortran caller must not read different words for one code).
    integer(kind = ik), parameter :: MSG_LEN = STATUS_MESSAGE_LEN + 1_ik
    integer(kind = ik), parameter :: N_MSG   = 14_ik

    !> Codes in column order; column N_MSG is the unknown-code fallback and has
    !! no entry here.
    integer(kind = ik), parameter :: MSG_CODE(N_MSG - 1_ik) = [ &
            SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, &
            FOS_ERROR_INVALID_C, FOS_ERROR_BEAK_SINGULARITY, &
            FOS_ERROR_CONVERGENCE, FOS_ERROR_BUFFER_MISMATCH]

    character(kind = c_char, len = MSG_LEN), parameter :: MSG_TEXT(N_MSG) = &
            [character(kind = c_char, len = MSG_LEN) :: &
                    'valid' // c_null_char, &
                    'too many parameters for this tier' // c_null_char, &
                    'cache not initialized' // c_null_char, &
                    'theta grid below minimum size (2)' // c_null_char, &
                    'params length differs from n_params' // c_null_char, &
                    'invalid init arguments' // c_null_char, &
                    'tables not initialized' // c_null_char, &
                    'rho <= 0 away from the poles' // c_null_char, &
                    'shape not star-convex from any origin' // c_null_char, &
                    'elongation c below the minimum' // c_null_char, &
                    'f_min below the beak threshold' // c_null_char, &
                    'iteration did not converge' // c_null_char, &
                    'output buffer size mismatch' // c_null_char, &
                    'unknown status code' // c_null_char]

    character(kind = c_char), save, target :: MSG_CHARS(MSG_LEN, N_MSG) = &
            reshape(transfer(MSG_TEXT, c_null_char, MSG_LEN * N_MSG), &
                    [MSG_LEN, N_MSG])

contains

    !===========================================================================
    ! DIAGNOSTICS
    !===========================================================================

    !> Pointer to the static, null-terminated description of `status`.
    !! Unknown codes get the fallback column. Never returns NULL.
    function fos_param_status_message(status) result(msg) &
            bind(c, name = 'fos_param_status_message')
        integer(c_int), value, intent(in) :: status
        type(c_ptr) :: msg
        integer(kind = ik) :: i, idx

        idx = N_MSG
        do i = 1_ik, N_MSG - 1_ik
            if (status == int(MSG_CODE(i), c_int)) then
                idx = i
                exit
            end if
        end do
        msg = c_loc(MSG_CHARS(1_ik, idx))
    end function fos_param_status_message

    !===========================================================================
    ! TABLES LIFECYCLE
    !===========================================================================

    !> Build the shared immutable level. Null handle on failure.
    function fos_param_tables_create(n_points, thetas, n_theta) result(handle) &
            bind(c, name = 'fos_param_tables_create')
        integer(c_int), value, intent(in) :: n_points, n_theta
        real(c_double), intent(in) :: thetas(n_theta)
        type(c_ptr) :: handle

        type(tables_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: thetas_f(max(n_theta, 0_c_int))

        handle = c_null_ptr
        if (n_theta < 1_c_int) return

        thetas_f = real(thetas, rk)
        allocate(p)
        call tables_init_s(p, int(n_points, ik), thetas_f, st)
        if (st == SHAPE_VALID) then
            handle = c_loc(p)
        else
            call tables_free_s(p)
            deallocate(p)
        end if
    end function fos_param_tables_create

    !> Release a tables handle. NULL-safe.
    subroutine fos_param_tables_destroy(handle) &
            bind(c, name = 'fos_param_tables_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(tables_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call tables_free_s(p)
        deallocate(p)
    end subroutine fos_param_tables_destroy

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Create a cache owning private tables. Null on failure.
    function fos_param_cache_create(n_params, n_points, thetas, n_theta) &
            result(handle) bind(c, name = 'fos_param_cache_create')
        integer(c_int), value, intent(in) :: n_params, n_points, n_theta
        real(c_double), intent(in) :: thetas(n_theta)
        type(c_ptr) :: handle

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: thetas_f(max(n_theta, 0_c_int))

        handle = c_null_ptr
        if (n_theta < 0_c_int) return

        thetas_f = real(thetas, rk)
        allocate(p)
        call cache_init_s(p, int(n_params, ik), int(n_points, ik), thetas_f, st)
        if (st == SHAPE_VALID) then
            handle = c_loc(p)
        else
            call cache_free_s(p)   ! a failed init may still own tables
            deallocate(p)
        end if
    end function fos_param_cache_create

    !> Create a cache bound to caller-owned shared tables, which must outlive it.
    !! Null on failure.
    function fos_param_cache_create_shared(tables, n_params) result(handle) &
            bind(c, name = 'fos_param_cache_create_shared')
        type(c_ptr), value, intent(in) :: tables
        integer(c_int), value, intent(in) :: n_params
        type(c_ptr) :: handle

        type(cache_t),  pointer :: p
        type(tables_t), pointer :: tp
        integer(kind = ik) :: st

        handle = c_null_ptr
        if (.not. c_associated(tables)) return

        call c_f_pointer(tables, tp)
        allocate(p)
        call cache_init_shared_s(p, tp, int(n_params, ik), st)
        if (st == SHAPE_VALID) then
            handle = c_loc(p)
        else
            call cache_free_s(p)
            deallocate(p)
        end if
    end function fos_param_cache_create_shared

    !> Release a cache. NULL-safe. Shared tables are left to their owner.
    subroutine fos_param_cache_destroy(handle) &
            bind(c, name = 'fos_param_cache_destroy')
        type(c_ptr), value, intent(in) :: handle
        type(cache_t), pointer :: p
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, p)
        call cache_free_s(p)
        deallocate(p)
    end subroutine fos_param_cache_destroy

    !===========================================================================
    ! CACHED COMPUTES (tier 2)
    !===========================================================================

    !> R(theta) at the cache's own thetas. n_radii must equal the handle's n_theta.
    function fos_param_cache_radius_grid(cache, params, n_params, radii, n_radii) &
            result(status) bind(c, name = 'fos_param_cache_radius_grid')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: radii_f(max(n_radii, 0_c_int))

        radii = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_radius_grid_s(p, params_f, radii_f, st)
        radii = real(radii_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_radius_grid

    !> R and dR/dtheta at the cache's own thetas.
    function fos_param_cache_radius_and_derivative(cache, params, n_params, &
            radii, dr_dtheta, n_radii) &
            result(status) bind(c, name = 'fos_param_cache_radius_and_derivative')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_radii
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: radii(n_radii), dr_dtheta(n_radii)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: radii_f(max(n_radii, 0_c_int)), dr_f(max(n_radii, 0_c_int))

        radii     = 0.0_c_double
        dr_dtheta = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_radius_and_derivative_s(p, params_f, radii_f, dr_f, st)
        radii     = real(radii_f, c_double)
        dr_dtheta = real(dr_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_radius_and_derivative

    !> R and dR/dtheta at caller thetas, against the cache's resolved shape.
    function fos_param_cache_radius_and_derivative_at_thetas(cache, params, &
            n_params, thetas, n_thetas, radii, dr_dtheta) result(status) &
            bind(c, name = 'fos_param_cache_radius_and_derivative_at_thetas')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_thetas
        real(c_double), intent(in)  :: params(n_params), thetas(n_thetas)
        real(c_double), intent(out) :: radii(n_thetas), dr_dtheta(n_thetas)
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: thetas_f(max(n_thetas, 0_c_int))
        real(kind = rk) :: radii_f(max(n_thetas, 0_c_int))
        real(kind = rk) :: dr_f(max(n_thetas, 0_c_int))

        radii     = 0.0_c_double
        dr_dtheta = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        thetas_f = real(thetas, rk)
        call cache_radius_and_derivative_at_thetas_s(p, params_f, thetas_f, &
                radii_f, dr_f, st)
        radii     = real(radii_f, c_double)
        dr_dtheta = real(dr_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_radius_and_derivative_at_thetas

    !> Resolve step: total z-shift and the analytic pole radii.
    function fos_param_cache_shape(cache, params, n_params, z_shift, r_north, r_south) &
            result(status) bind(c, name = 'fos_param_cache_shape')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_shift, r_north, r_south
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_shift_f, r_north_f, r_south_f

        z_shift = 0.0_c_double
        r_north = 0.0_c_double
        r_south = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_shape_s(p, params_f, z_shift_f, r_north_f, r_south_f, st)
        z_shift = real(z_shift_f, c_double)
        r_north = real(r_north_f, c_double)
        r_south = real(r_south_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_shape

    !> Cylindrical rho(z) in the COM frame. n_z must equal the handle's n_points.
    function fos_param_cache_rho_z_grid(cache, params, n_params, z, rho, drho_dz, &
            n_z, z_shift) result(status) bind(c, name = 'fos_param_cache_rho_z_grid')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params, n_z
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z(n_z), rho(n_z), drho_dz(n_z)
        real(c_double), intent(out) :: z_shift
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_f(max(n_z, 0_c_int)), rho_f(max(n_z, 0_c_int))
        real(kind = rk) :: drho_f(max(n_z, 0_c_int))
        real(kind = rk) :: z_shift_f

        z       = 0.0_c_double
        rho     = 0.0_c_double
        drho_dz = 0.0_c_double
        z_shift = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_rho_z_grid_s(p, params_f, z_f, rho_f, drho_f, z_shift_f, st)
        z       = real(z_f, c_double)
        rho     = real(rho_f, c_double)
        drho_dz = real(drho_f, c_double)
        z_shift = real(z_shift_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_rho_z_grid

    !> Neck in the COM frame; `found` is 0 or 1.
    function fos_param_cache_neck(cache, params, n_params, z_neck, rho_neck, found) &
            result(status) bind(c, name = 'fos_param_cache_neck')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_neck, rho_neck
        integer(c_int), intent(out) :: found
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_neck_f, rho_neck_f
        logical :: found_f

        z_neck   = 0.0_c_double
        rho_neck = 0.0_c_double
        found    = 0_c_int
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_neck_s(p, params_f, z_neck_f, rho_neck_f, found_f, st)
        z_neck   = real(z_neck_f, c_double)
        rho_neck = real(rho_neck_f, c_double)
        if (found_f) found = 1_c_int
        status = int(st, c_int)
    end function fos_param_cache_neck

    !> Star-convexity optimum: total shift and g(s*).
    function fos_param_cache_star_convexity_optimum(cache, params, n_params, &
            z_shift_total, g_opt) result(status) &
            bind(c, name = 'fos_param_cache_star_convexity_optimum')
        type(c_ptr), value, intent(in) :: cache
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_shift_total, g_opt
        integer(c_int) :: status

        type(cache_t), pointer :: p
        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_shift_f, g_opt_f

        z_shift_total = 0.0_c_double
        g_opt         = 0.0_c_double
        if (.not. c_associated(cache)) then
            status = int(SHAPE_ERROR_CACHE_NOT_INITIALIZED, c_int)
            return
        end if
        call c_f_pointer(cache, p)
        params_f = real(params, rk)
        call cache_star_convexity_optimum_s(p, params_f, z_shift_f, g_opt_f, st)
        z_shift_total = real(z_shift_f, c_double)
        g_opt         = real(g_opt_f, c_double)
        status = int(st, c_int)
    end function fos_param_cache_star_convexity_optimum

    !===========================================================================
    ! FLAT TIER-1 COMPUTES — trailing nullable status
    !===========================================================================

    !> One-shot R(theta) at caller thetas; builds and discards its own tables.
    subroutine fos_param_radius_grid(params, n_params, thetas, n_thetas, &
            n_points, radii, status) bind(c, name = 'fos_param_radius_grid')
        integer(c_int), value, intent(in) :: n_params, n_thetas, n_points
        real(c_double), intent(in)  :: params(n_params), thetas(n_thetas)
        real(c_double), intent(out) :: radii(n_thetas)
        integer(c_int), intent(out), optional :: status   ! NULL-able from C

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: thetas_f(max(n_thetas, 0_c_int))
        real(kind = rk) :: radii_f(max(n_thetas, 0_c_int))

        radii = 0.0_c_double
        if (n_params < 0_c_int .or. n_thetas < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        thetas_f = real(thetas, rk)
        call compute_radius_grid_standalone_s(params_f, thetas_f, &
                int(n_points, ik), radii_f, st)
        radii = real(radii_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_radius_grid

    !> One-shot R and dR/dtheta at caller thetas.
    subroutine fos_param_radius_and_derivative(params, n_params, thetas, n_thetas, &
            n_points, radii, dr_dtheta, status) &
            bind(c, name = 'fos_param_radius_and_derivative')
        integer(c_int), value, intent(in) :: n_params, n_thetas, n_points
        real(c_double), intent(in)  :: params(n_params), thetas(n_thetas)
        real(c_double), intent(out) :: radii(n_thetas), dr_dtheta(n_thetas)
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: thetas_f(max(n_thetas, 0_c_int))
        real(kind = rk) :: radii_f(max(n_thetas, 0_c_int))
        real(kind = rk) :: dr_f(max(n_thetas, 0_c_int))

        radii     = 0.0_c_double
        dr_dtheta = 0.0_c_double
        if (n_params < 0_c_int .or. n_thetas < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        thetas_f = real(thetas, rk)
        call compute_radius_and_derivative_standalone_s(params_f, thetas_f, &
                int(n_points, ik), radii_f, dr_f, st)
        radii     = real(radii_f, c_double)
        dr_dtheta = real(dr_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_radius_and_derivative

    !> One-shot resolve step: total z-shift and the analytic pole radii.
    subroutine fos_param_shape(params, n_params, n_points, z_shift, r_north, &
            r_south, status) bind(c, name = 'fos_param_shape')
        integer(c_int), value, intent(in) :: n_params, n_points
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_shift, r_north, r_south
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_shift_f, r_north_f, r_south_f

        z_shift = 0.0_c_double
        r_north = 0.0_c_double
        r_south = 0.0_c_double
        if (n_params < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_shape_standalone_s(params_f, int(n_points, ik), z_shift_f, &
                r_north_f, r_south_f, st)
        z_shift = real(z_shift_f, c_double)
        r_north = real(r_north_f, c_double)
        r_south = real(r_south_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_shape

    !> One-shot cylindrical rho(z) in the COM frame; buffers are n_points long.
    subroutine fos_param_rho_z_grid(params, n_params, n_points, z, rho, drho_dz, &
            z_shift, status) bind(c, name = 'fos_param_rho_z_grid')
        integer(c_int), value, intent(in) :: n_params, n_points
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z(n_points), rho(n_points), drho_dz(n_points)
        real(c_double), intent(out) :: z_shift
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_f(max(n_points, 0_c_int)), rho_f(max(n_points, 0_c_int))
        real(kind = rk) :: drho_f(max(n_points, 0_c_int))
        real(kind = rk) :: z_shift_f

        z       = 0.0_c_double
        rho     = 0.0_c_double
        drho_dz = 0.0_c_double
        z_shift = 0.0_c_double
        if (n_params < 0_c_int .or. n_points < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_rho_z_grid_standalone_s(params_f, int(n_points, ik), z_f, &
                rho_f, drho_f, z_shift_f, st)
        z       = real(z_f, c_double)
        rho     = real(rho_f, c_double)
        drho_dz = real(drho_f, c_double)
        z_shift = real(z_shift_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_rho_z_grid

    !> One-shot neck in the COM frame; `found` is 0 or 1.
    subroutine fos_param_neck(params, n_params, n_points, z_neck, rho_neck, &
            found, status) bind(c, name = 'fos_param_neck')
        integer(c_int), value, intent(in) :: n_params, n_points
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_neck, rho_neck
        integer(c_int), intent(out) :: found
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_neck_f, rho_neck_f
        logical :: found_f

        z_neck   = 0.0_c_double
        rho_neck = 0.0_c_double
        found    = 0_c_int
        if (n_params < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_neck_standalone_s(params_f, int(n_points, ik), z_neck_f, &
                rho_neck_f, found_f, st)
        z_neck   = real(z_neck_f, c_double)
        rho_neck = real(rho_neck_f, c_double)
        if (found_f) found = 1_c_int
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_neck

    !> One-shot star-convexity optimum: total shift and g(s*).
    subroutine fos_param_star_convexity_optimum(params, n_params, n_points, &
            z_shift_total, g_opt, status) &
            bind(c, name = 'fos_param_star_convexity_optimum')
        integer(c_int), value, intent(in) :: n_params, n_points
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_shift_total, g_opt
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_shift_f, g_opt_f

        z_shift_total = 0.0_c_double
        g_opt         = 0.0_c_double
        if (n_params < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_star_convexity_optimum_standalone_s(params_f, &
                int(n_points, ik), z_shift_f, g_opt_f, st)
        z_shift_total = real(z_shift_f, c_double)
        g_opt         = real(g_opt_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_star_convexity_optimum

    !> Intrinsic COM z-shift (closed form).
    subroutine fos_param_z_shift(params, n_params, z_shift, status) &
            bind(c, name = 'fos_param_z_shift')
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: z_shift
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: z_shift_f

        z_shift = 0.0_c_double
        if (n_params < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_z_shift_s(params_f, z_shift_f, st)
        z_shift = real(z_shift_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_z_shift

    !> a2 from the volume-conservation constraint.
    subroutine fos_param_a2(params, n_params, a2, status) &
            bind(c, name = 'fos_param_a2')
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), intent(out) :: a2
        integer(c_int), intent(out), optional :: status

        integer(kind = ik) :: st
        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: a2_f

        a2 = 0.0_c_double
        if (n_params < 0_c_int) then
            if (present(status)) status = int(SHAPE_ERROR_INVALID_INIT, c_int)
            return
        end if
        params_f = real(params, rk)
        call compute_a2_s(params_f, a2_f, st)
        a2 = real(a2_f, c_double)
        if (present(status)) status = int(st, c_int)
    end subroutine fos_param_a2

    !===========================================================================
    ! RAW EVALUATOR — outside the tier rules
    !===========================================================================

    !> rho and drho/dz at one z in the SHIFTED frame. No validation, no status.
    subroutine fos_param_rho_at_z(params, n_params, z, z_shift, rho, drho_dz) &
            bind(c, name = 'fos_param_rho_at_z')
        integer(c_int), value, intent(in) :: n_params
        real(c_double), intent(in)  :: params(n_params)
        real(c_double), value, intent(in) :: z, z_shift
        real(c_double), intent(out) :: rho, drho_dz

        real(kind = rk) :: params_f(max(n_params, 0_c_int))
        real(kind = rk) :: rho_f, drho_f

        rho     = 0.0_c_double
        drho_dz = 0.0_c_double
        if (n_params < 0_c_int) return
        params_f = real(params, rk)
        call compute_rho_at_z_s(params_f, real(z, rk), real(z_shift, rk), rho_f, drho_f)
        rho     = real(rho_f, c_double)
        drho_dz = real(drho_f, c_double)
    end subroutine fos_param_rho_at_z

end module fos_parameterization_c_api_mod
