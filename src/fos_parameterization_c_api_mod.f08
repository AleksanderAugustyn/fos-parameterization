!> C-interop API for the FoS parameterization library.
!!
!! Standalone functions only — FoS has no shape-independent precompute, so
!! there is no cache/handle tier. Status codes 0-4 match the Fortran FOS_*
!! parameters and the C FOS_* macros; FOS_ERROR_INVALID_ARGUMENTS (5) exists
!! only at this layer (bad n_grid / n_z — the Fortran API takes assumed-size
!! arrays and cannot receive these).
module fos_parameterization_c_api_mod

    use c_bindings_mod, only: ik_c, rk_c, c_char, c_null_char
    use precision_utilities_mod, only: ik, rk
    use fos_parameterization_mod, only: &
            compute_fos_radius_grid_s, compute_rho_at_z_s, compute_fos_neck_s, &
            compute_fos_z_shift_f, compute_fos_a2_f, &
            compute_fos_shape_s, compute_fos_radius_and_derivative_at_thetas_s, &
            FOS_VALID, FOS_ERROR_INVALID_C, C_MIN

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

        real(kind = rk), allocatable :: f_params(:), f_radii(:)
        real(kind = rk)      :: f_z_shift
        logical              :: is_valid
        integer(kind = ik)   :: error_code
        character(len = 256) :: f_message

        radii = 0.0_rk_c
        z_shift = 0.0_rk_c

        if (n_grid < 2_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            call marshal_message_to_c('n_grid must be >= 2', message_buf, message_buf_len)
            return
        end if

        allocate(f_params(int(n_params, ik)))
        allocate(f_radii(int(n_grid, ik)), source = 0.0_rk)
        f_params(:) = real(params(:), rk)
        f_z_shift = 0.0_rk
        error_code = FOS_VALID
        f_message = ''

        call compute_fos_radius_grid_s(f_params, int(n_grid, ik), f_radii, f_z_shift, &
                is_valid, f_message, error_code = error_code)

        radii(:) = real(f_radii(:), rk_c)
        z_shift = real(f_z_shift, rk_c)
        status = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
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
        integer(kind = ik) :: i, n

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
        z_sh = compute_fos_z_shift_f(f_params)
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
        real(kind = rk) :: f_z_neck, f_rho_neck
        logical         :: l_found

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

        call compute_fos_neck_s(f_params, f_z_neck, f_rho_neck, l_found)

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
        logical              :: is_valid
        integer(kind = ik)   :: error_code
        character(len = 256) :: f_message

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

        call compute_fos_shape_s(f_params, int(n_rho_grid, ik), f_z_shift, &
                f_r_north, f_r_south, is_valid, f_message, error_code)

        z_shift = real(f_z_shift, rk_c)
        r_north = real(f_r_north, rk_c)
        r_south = real(f_r_south, rk_c)
        status  = int(error_code, ik_c)
        call marshal_message_to_c(f_message, message_buf, message_buf_len)
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

        if (n_thetas < 1_ik_c) then
            status = int(FOS_ERROR_INVALID_ARGUMENTS, ik_c)
            return
        end if

        radii(:)     = 0.0_rk_c
        dr_dtheta(:) = 0.0_rk_c

        allocate(f_params(int(n_params, ik)), f_thetas(int(n_thetas, ik)))
        allocate(f_radii(int(n_thetas, ik)), f_dr(int(n_thetas, ik)))
        f_params(:) = real(params(:), rk)
        f_thetas(:) = real(thetas(:), rk)

        call compute_fos_radius_and_derivative_at_thetas_s( &
                f_params, f_thetas, real(z_shift, rk), f_radii, f_dr)

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

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)
        z_sh = real(compute_fos_z_shift_f(f_params), rk_c)
    end function fos_z_shift

    function fos_a2(params, n_params) result(a2) bind(c, name='fos_a2')
        integer(kind = ik_c), intent(in), value :: n_params
        real(kind = rk_c),    intent(in)        :: params(n_params)
        real(kind = rk_c) :: a2

        real(kind = rk), allocatable :: f_params(:)

        allocate(f_params(int(n_params, ik)))
        f_params(:) = real(params(:), rk)
        a2 = real(compute_fos_a2_f(f_params), rk_c)
    end function fos_a2

end module fos_parameterization_c_api_mod
