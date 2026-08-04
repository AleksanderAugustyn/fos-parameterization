!> C-interop API for the FoS parameterization library.
!!
!! Standalone functions only — this layer still exposes the 1.x `fos_*` symbol
!! set, now implemented on the 2.0 tier-1 (standalone) forms. Status codes 100+
!! match the Fortran FOS_* parameters and the C FOS_* macros;
!! FOS_ERROR_INVALID_ARGUMENTS (5) exists only at this layer (bad n_grid / n_z
!! — the Fortran API takes assumed-size arrays and cannot receive these).
!!
!! Internal rho(z) resolution: the tier-1 forms enforce a floor of
!! FOS_N_POINTS_FLOOR nodes, so every entry point here raises the caller's
!! request to `max(FOS_N_POINTS_FLOOR, n)`. 1.x used the caller's value
!! verbatim; below the floor the resolved z-shift therefore moves slightly.
!!
!! `fos_compute_radius_and_derivative_at_thetas` is the one call that takes an
!! externally resolved z_shift, which no tier-1 form does. It is built directly
!! on the worker evaluator (`fos_bundle_t` + `newton_radius_s`), with the
!! Newton bracket bound taken from a rho(z) grid at N_RHO_AT_THETAS nodes.
!! That grid is a GATED form, so the call now validates the shape where 1.x
!! evaluated unconditionally: a pinched interior propagates 100 and an
!! over-length vector propagates 1, both with the outputs zero-filled.
module fos_parameterization_c_api_mod

    use c_bindings_mod, only: ik_c, rk_c, c_char, c_null_char
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: &
            compute_rho_at_z_s, compute_a2_s, compute_z_shift_s, &
            compute_radius_grid_standalone_s, compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, compute_neck_standalone_s, &
            status_message, FOS_N_POINTS_FLOOR, &
            FOS_VALID, FOS_ERROR_INVALID_C, C_MIN
    use fos_parameterization_workers_mod, only: fos_bundle_t, newton_radius_s, &
            FOS_MAX_K

    implicit none

    private

    public :: fos_compute_radius_grid
    public :: fos_compute_rho_profile
    public :: fos_compute_neck
    public :: fos_z_shift
    public :: fos_a2
    public :: fos_compute_shape
    public :: fos_compute_radius_and_derivative_at_thetas

    !> C-API-level status for invalid grid/profile sizes (see module docstring).
    integer(kind = ik), parameter :: FOS_ERROR_INVALID_ARGUMENTS = 5_ik

    !> rho(z) resolution behind `fos_compute_radius_and_derivative_at_thetas`,
    !! used only for the Newton bracket bound (rho_max).
    integer(kind = ik), parameter :: N_RHO_AT_THETAS = FOS_N_POINTS_FLOOR

    !> Neck scan resolution — the 1.x scanner's default, kept unchanged.
    integer(kind = ik), parameter :: NECK_N_POINTS = 1001_ik

contains

    !===========================================================================
    ! INTERNAL HELPER — Fortran string -> C buffer (null-terminated, truncated)
    !===========================================================================

    pure subroutine marshal_message_to_c(f_message, c_buf, c_buf_len)
        character(len = *),         intent(in)  :: f_message
        integer(kind = ik_c),       intent(in)  :: c_buf_len
        character(kind = c_char),   intent(out) :: c_buf(c_buf_len)

        integer(kind = ik) :: i, msg_len, max_copy

        if (c_buf_len < 1_ik_c) return
        msg_len  = len_trim(f_message)
        max_copy = min(int(msg_len, ik), int(c_buf_len, ik) - 1_ik)
        do i = 1_ik, max_copy
            c_buf(i) = f_message(i:i)
        end do
        c_buf(max_copy + 1_ik) = c_null_char
    end subroutine marshal_message_to_c

    !===========================================================================
    ! RADIUS GRID
    !===========================================================================

    function fos_compute_radius_grid( &
            params, n_params, n_grid, radii, z_shift, message_buf_len, message_buf) &
            result(status) bind(c, name='fos_compute_radius_grid')

        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        integer(kind = ik_c),     intent(in), value :: n_grid
        real(kind = rk_c),        intent(out)       :: radii(n_grid)
        real(kind = rk_c),        intent(out)       :: z_shift
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:), f_radii(:), f_thetas(:)
        real(kind = rk)      :: f_z_shift, r_north, r_south
        integer(kind = ik)   :: error_code, n, n_rho, i

        radii = 0.0_rk_c
        z_shift = 0.0_rk_c

        if (n_grid < 2_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            call marshal_message_to_c('n_grid must be >= 2', message_buf, message_buf_len)
            return
        end if

        n = int(n_grid, ik)
        n_rho = max(FOS_N_POINTS_FLOOR, n)

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(n), source = 0.0_rk)
        allocate(f_thetas(n))
        f_params(:) = real(params(:), rk)
        ! Uniform nodes with the endpoints pinned: under -ffast-math the
        ! product (n-1)*PI_C/(n-1) can land one ulp ABOVE PI_C, which the
        ! theta-domain check then rejects with SHAPE_ERROR_INVALID_GRID.
        do i = 1_ik, n
            f_thetas(i) = real(i - 1_ik, rk) * PI_C / real(n - 1_ik, rk)
        end do
        f_thetas(1) = 0.0_rk
        f_thetas(n) = PI_C

        ! Resolve first: the z_shift is reported even though the radius form
        ! resolves again internally (identical inputs, identical answer).
        call compute_shape_standalone_s(f_params, n_rho, f_z_shift, r_north, &
                r_south, error_code)
        if (error_code /= FOS_VALID) then
            status = int(error_code, ik_c)
            call marshal_message_to_c(status_message(error_code), message_buf, &
                    message_buf_len)
            return
        end if

        call compute_radius_grid_standalone_s(f_params, f_thetas, n_rho, f_radii, &
                error_code)
        if (error_code /= FOS_VALID) then
            status = int(error_code, ik_c)
            call marshal_message_to_c(status_message(error_code), message_buf, &
                    message_buf_len)
            return
        end if

        radii(:) = real(f_radii(:), rk_c)
        z_shift = real(f_z_shift, rk_c)
        status = int(FOS_VALID, ik_c)
        call marshal_message_to_c('', message_buf, message_buf_len)
    end function fos_compute_radius_grid

    !===========================================================================
    ! RHO(Z) PROFILE (COM frame)
    !===========================================================================

    function fos_compute_rho_profile( &
            params, n_params, n_z, z, rho, drho_dz, message_buf_len, message_buf) &
            result(status) bind(c, name='fos_compute_rho_profile')

        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        integer(kind = ik_c),     intent(in), value :: n_z
        real(kind = rk_c),        intent(out)       :: z(n_z)
        real(kind = rk_c),        intent(out)       :: rho(n_z)
        real(kind = rk_c),        intent(out)       :: drho_dz(n_z)
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:)
        real(kind = rk)    :: c, z_sh, dz, f_z, f_rho, f_drho
        integer(kind = ik) :: i, n, code

        z = 0.0_rk_c
        rho = 0.0_rk_c
        drho_dz = 0.0_rk_c

        if (n_z < 2_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            call marshal_message_to_c('n_z must be >= 2', message_buf, message_buf_len)
            return
        end if

        if (n_params < 1_ik_c) then
            status = int(FOS_ERROR_INVALID_C, ik_c)
            call marshal_message_to_c('Empty parameter array', message_buf, message_buf_len)
            return
        end if

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)
        c = f_params(1)

        if (c <= C_MIN) then
            status = int(FOS_ERROR_INVALID_C, ik_c)
            call marshal_message_to_c('Elongation c must be positive', &
                    message_buf, message_buf_len)
            return
        end if

        n = int(n_z, ik)
        call compute_z_shift_s(f_params, z_sh, code)
        if (code /= FOS_VALID) then
            status = int(code, ik_c)
            call marshal_message_to_c(status_message(code), message_buf, message_buf_len)
            return
        end if
        dz = 2.0_rk * c / real(n - 1_ik, rk)

        do i = 1_ik, n
            f_z = -c + z_sh + real(i - 1_ik, rk) * dz
            call compute_rho_at_z_s(f_params, f_z, z_sh, f_rho, f_drho)
            z(i) = real(f_z, rk_c)
            rho(i) = real(f_rho, rk_c)
            drho_dz(i) = real(f_drho, rk_c)
        end do

        status = int(FOS_VALID, ik_c)
        call marshal_message_to_c('', message_buf, message_buf_len)
    end function fos_compute_rho_profile

    !===========================================================================
    ! NECK
    !===========================================================================

    function fos_compute_neck(params, n_params, z_neck, rho_neck, found) &
            result(status) bind(c, name='fos_compute_neck')

        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c),    intent(in)        :: params(n_params)
        real(kind = rk_c),    intent(out)       :: z_neck
        real(kind = rk_c),    intent(out)       :: rho_neck
        integer(kind = ik_c), intent(out)       :: found
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:)
        real(kind = rk)    :: f_z_neck, f_rho_neck
        integer(kind = ik) :: code
        logical            :: l_found

        z_neck = 0.0_rk_c
        rho_neck = 0.0_rk_c
        found = 0_ik_c

        if (n_params < 1_ik_c) then
            status = int(FOS_ERROR_INVALID_C, ik_c)
            return
        end if

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)

        if (f_params(1) <= C_MIN) then
            status = int(FOS_ERROR_INVALID_C, ik_c)
            return
        end if

        call compute_neck_standalone_s(f_params, NECK_N_POINTS, f_z_neck, &
                f_rho_neck, l_found, code)
        if (code /= FOS_VALID) then
            status = int(code, ik_c)
            return
        end if

        z_neck = real(f_z_neck, rk_c)
        rho_neck = real(f_rho_neck, rk_c)
        if (l_found) found = 1_ik_c
        status = int(FOS_VALID, ik_c)
    end function fos_compute_neck

    !===========================================================================
    ! SHAPE SPLIT + DERIVATIVE EVALUATION
    !===========================================================================

    function fos_compute_shape( &
            params, n_params, n_rho_grid, z_shift, r_north, r_south, &
            message_buf_len, message_buf) &
            result(status) bind(c, name='fos_compute_shape')

        integer(kind = ik_c),     intent(in), value :: n_params
        real(kind = rk_c),        intent(in)        :: params(n_params)
        integer(kind = ik_c),     intent(in), value :: n_rho_grid
        real(kind = rk_c),        intent(out)       :: z_shift
        real(kind = rk_c),        intent(out)       :: r_north
        real(kind = rk_c),        intent(out)       :: r_south
        integer(kind = ik_c),     intent(in), value :: message_buf_len
        character(kind = c_char), intent(out)       :: message_buf(message_buf_len)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:)
        real(kind = rk)      :: f_z_shift, f_r_north, f_r_south
        integer(kind = ik)   :: error_code

        z_shift = 0.0_rk_c
        r_north = 0.0_rk_c
        r_south = 0.0_rk_c

        if (n_rho_grid < 2_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            call marshal_message_to_c('n_rho_grid must be >= 2', message_buf, message_buf_len)
            return
        end if

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)

        call compute_shape_standalone_s(f_params, &
                max(FOS_N_POINTS_FLOOR, int(n_rho_grid, ik)), f_z_shift, &
                f_r_north, f_r_south, error_code)

        z_shift = real(f_z_shift, rk_c)
        r_north = real(f_r_north, rk_c)
        r_south = real(f_r_south, rk_c)
        status  = int(error_code, ik_c)
        if (error_code == FOS_VALID) then
            call marshal_message_to_c('', message_buf, message_buf_len)
        else
            call marshal_message_to_c(status_message(error_code), message_buf, &
                    message_buf_len)
        end if
    end function fos_compute_shape

    function fos_compute_radius_and_derivative_at_thetas( &
            params, n_params, thetas, n_thetas, z_shift, radii, dr_dtheta) &
            result(status) bind(c, name='fos_compute_radius_and_derivative_at_thetas')

        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c),    intent(in)        :: params(n_params)
        integer(kind = ik_c), intent(in), value :: n_thetas
        real(kind = rk_c),    intent(in)        :: thetas(n_thetas)
        real(kind = rk_c),    intent(in), value :: z_shift
        real(kind = rk_c),    intent(out)       :: radii(n_thetas)
        real(kind = rk_c),    intent(out)       :: dr_dtheta(n_thetas)
        integer(kind = ik_c) :: status

        real(kind = rk), allocatable :: f_params(:), f_thetas(:), f_radii(:), f_dr(:)
        real(kind = rk) :: gz(N_RHO_AT_THETAS), grho(N_RHO_AT_THETAS)
        real(kind = rk) :: gdrho(N_RHO_AT_THETAS)
        real(kind = rk) :: grid_shift, rho_max, shift, c
        type(fos_bundle_t) :: bundle
        logical, allocatable :: converged(:)
        integer(kind = ik) :: code, n, i

        if (n_thetas < 1_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            return
        end if

        radii(:)     = 0.0_rk_c
        dr_dtheta(:) = 0.0_rk_c

        n = int(n_thetas, ik)
        allocate(f_params(int(n_params, ik)), f_thetas(n))
        allocate(f_radii(n), f_dr(n), converged(n))
        f_params(:) = real(params(:), rk)
        f_thetas(:) = real(thetas(:), rk)
        shift = real(z_shift, rk)

        ! The caller supplies the shift, so no tier-1 form applies; the Newton
        ! bracket bound still needs rho_max, which comes from a rho(z) grid.
        call compute_rho_z_grid_standalone_s(f_params, N_RHO_AT_THETAS, gz, grho, &
                gdrho, grid_shift, code)
        if (code /= FOS_VALID) then
            status = int(code, ik_c)
            return
        end if

        rho_max = 0.0_rk
        do i = 1_ik, N_RHO_AT_THETAS
            rho_max = max(rho_max, grho(i))
        end do

        c = f_params(1)
        bundle%n_params = min(int(n_params, ik), FOS_MAX_K)
        bundle%params(:) = 0.0_rk
        bundle%params(1:bundle%n_params) = f_params(1:bundle%n_params)
        bundle%z_shift = shift
        bundle%r_hi_bound = 2.0_rk * sqrt(rho_max**2 &
                + max(c + shift, abs(-c + shift))**2)

        call newton_radius_s(bundle, cos(f_thetas), f_radii, f_dr, converged)

        radii(:)     = real(f_radii(:), rk_c)
        dr_dtheta(:) = real(f_dr(:), rk_c)
        status = int(FOS_VALID, ik_c)
    end function fos_compute_radius_and_derivative_at_thetas

    !===========================================================================
    ! SCALAR HELPERS
    !===========================================================================

    function fos_z_shift(params, n_params) result(z_sh) bind(c, name='fos_z_shift')
        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c),    intent(in)        :: params(n_params)
        real(kind = rk_c) :: z_sh

        real(kind = rk), allocatable :: f_params(:)
        real(kind = rk)    :: f_z_sh
        integer(kind = ik) :: code

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)
        call compute_z_shift_s(f_params, f_z_sh, code)
        if (code /= FOS_VALID) f_z_sh = 0.0_rk
        z_sh = real(f_z_sh, rk_c)
    end function fos_z_shift

    function fos_a2(params, n_params) result(a2) bind(c, name='fos_a2')
        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c),    intent(in)        :: params(n_params)
        real(kind = rk_c) :: a2

        real(kind = rk), allocatable :: f_params(:)
        real(kind = rk)    :: f_a2
        integer(kind = ik) :: code

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)
        call compute_a2_s(f_params, f_a2, code)
        if (code /= FOS_VALID) f_a2 = 0.0_rk
        a2 = real(f_a2, rk_c)
    end function fos_a2

end module fos_parameterization_c_api_mod
