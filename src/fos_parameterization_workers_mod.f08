!> Computational kernels for the Fourier-over-Spheroid (FoS) parameterization.
!!
!! This module holds the parameter-independent precomputation (`tables_t`) and,
!! as the library grows, the inner-loop math that consumes it. No derived-type
!! I/O, no `iso_c_binding`, no allocation in the hot path: everything a shape
!! evaluation needs is allocated once, at table-build time.
!!
!! `tables_t` caches what depends only on the grid, not on the shape
!! parameters: the u-grid, the Fourier basis sampled on it, the same basis on
!! the clamped beak-scan grid, and the caller's theta nodes.
module fos_parameterization_workers_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_INVALID_GRID

    implicit none

    private

    ! Types
    public :: tables_t

    ! Procedures
    public :: tables_init_s
    public :: tables_free_s
    public :: build_tables_s

    !---------------------------------------------------------------------------
    ! Table-construction parameters
    !---------------------------------------------------------------------------
    !! Minimum u-grid resolution. Below this the rho(z) grid is too coarse for
    !! the neck/star-convexity scans that read it.
    integer(kind = ik), parameter, public :: FOS_N_POINTS_FLOOR = 100_ik

    !! Fourier orders tabulated by `tables_init_s`: (8 + 2) / 2 + 1 for the
    !! 8-parameter FoS surface, matching `compute_fos_f_and_derivatives_s`.
    integer(kind = ik), parameter, public :: FOS_TABLES_K_MAX = 6_ik

    !! Sample count of the beak (f_min) scan grid.
    integer(kind = ik), parameter, public :: FOS_BEAK_SCAN_POINTS = 1001_ik

    !! Beak-scan clamp. f(±1) = 0 analytically, so an unclamped scan would find
    !! f_min = 0 for every shape — including the sphere — and reject it all.
    real(kind = rk), parameter :: U_BEAK_CLAMP = 0.999_rk

    !> Parameter-independent trigonometric tables for FoS shape evaluation.
    type, public :: tables_t
        integer(kind = ik) :: n_points = 0_ik, n_theta = 0_ik, k_max = 0_ik
        real(kind = rk), allocatable :: u(:)                              !! n_points
        real(kind = rk), allocatable :: cos_even(:, :), sin_even(:, :)    !! (n_points, k_max)
        real(kind = rk), allocatable :: cos_odd(:, :), sin_odd(:, :)      !! (n_points, k_max)
        real(kind = rk), allocatable :: bk_cos_even(:, :), bk_sin_odd(:, :) !! (FOS_BEAK_SCAN_POINTS, k_max)
        real(kind = rk), allocatable :: u_beak(:)                         !! FOS_BEAK_SCAN_POINTS
        real(kind = rk), allocatable :: thetas(:), cos_thetas(:)          !! n_theta
        real(kind = rk), allocatable :: omega(:), psi(:)                  !! k_max
        logical :: initialized = .false.
    end type tables_t

contains

    !> Builds the standard tables for a caller-supplied theta set.
    !!
    !! Requires at least one theta node; the theta-less tier-1 forms call
    !! `build_tables_s` directly instead.
    !!
    !! @param[out] tables    Freshly built tables (reset on any rejection)
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  thetas    Polar angles in [0, pi], at least one
    !! @param[out] status    SHAPE_VALID, or SHAPE_ERROR_INVALID_GRID
    subroutine tables_init_s(tables, n_points, thetas, status)

        type(tables_t), intent(out) :: tables
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        if (n_points < FOS_N_POINTS_FLOOR) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        if (size(thetas, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call build_tables_s(tables, n_points, thetas, FOS_TABLES_K_MAX, status)

    end subroutine tables_init_s

    !> Releases every table array and clears the initialized flag. Infallible.
    !!
    !! @param[in,out] tables  Tables to reset; safe on an already-free instance
    subroutine tables_free_s(tables)

        type(tables_t), intent(inout) :: tables

        if (allocated(tables%u)) deallocate(tables%u)
        if (allocated(tables%cos_even)) deallocate(tables%cos_even)
        if (allocated(tables%sin_even)) deallocate(tables%sin_even)
        if (allocated(tables%cos_odd)) deallocate(tables%cos_odd)
        if (allocated(tables%sin_odd)) deallocate(tables%sin_odd)
        if (allocated(tables%bk_cos_even)) deallocate(tables%bk_cos_even)
        if (allocated(tables%bk_sin_odd)) deallocate(tables%bk_sin_odd)
        if (allocated(tables%u_beak)) deallocate(tables%u_beak)
        if (allocated(tables%thetas)) deallocate(tables%thetas)
        if (allocated(tables%cos_thetas)) deallocate(tables%cos_thetas)
        if (allocated(tables%omega)) deallocate(tables%omega)
        if (allocated(tables%psi)) deallocate(tables%psi)

        tables%n_points = 0_ik
        tables%n_theta = 0_ik
        tables%k_max = 0_ik
        tables%initialized = .false.

    end subroutine tables_free_s

    !> Table-construction worker: u-grid, Fourier bases, beak grid, theta nodes.
    !!
    !! Unlike `tables_init_s` this accepts an empty theta set — the theta-less
    !! tier-1 forms (shape, rho_z_grid, neck, optimum) need no theta nodes, and
    !! then `thetas`/`cos_thetas` stay unallocated with `n_theta = 0`.
    !!
    !! @param[out] tables    Freshly built tables (reset on any rejection)
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  thetas    Polar angles in [0, pi]; may be empty
    !! @param[in]  k_max     Number of Fourier orders to tabulate
    !! @param[out] status    SHAPE_VALID, or SHAPE_ERROR_INVALID_GRID
    subroutine build_tables_s(tables, n_points, thetas, k_max, status)

        type(tables_t), intent(out) :: tables
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: k_max
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i, k, n_theta
        real(kind = rk) :: u_raw

        if (n_points < FOS_N_POINTS_FLOOR) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        n_theta = size(thetas, kind = ik)
        do i = 1_ik, n_theta
            if (thetas(i) < 0.0_rk .or. thetas(i) > PI_C) then
                status = SHAPE_ERROR_INVALID_GRID
                return
            end if
        end do

        tables%n_points = n_points
        tables%n_theta = n_theta
        tables%k_max = k_max

        allocate(tables%omega(k_max), tables%psi(k_max))
        do k = 1_ik, k_max
            tables%omega(k) = real(2_ik * k - 1_ik, rk) * PI_C / 2.0_rk
            tables%psi(k) = real(k, rk) * PI_C
        end do

        ! u-grid and its Fourier basis
        allocate(tables%u(n_points))
        do i = 1_ik, n_points
            tables%u(i) = -1.0_rk + 2.0_rk * real(i - 1_ik, rk) / real(n_points - 1_ik, rk)
        end do

        allocate(tables%cos_even(n_points, k_max), tables%sin_even(n_points, k_max))
        allocate(tables%cos_odd(n_points, k_max), tables%sin_odd(n_points, k_max))
        do k = 1_ik, k_max
            do i = 1_ik, n_points
                tables%cos_even(i, k) = cos(tables%omega(k) * tables%u(i))
                tables%sin_even(i, k) = sin(tables%omega(k) * tables%u(i))
                tables%cos_odd(i, k) = cos(tables%psi(k) * tables%u(i))
                tables%sin_odd(i, k) = sin(tables%psi(k) * tables%u(i))
            end do
        end do

        ! Beak-scan grid: same construction, clamped away from the poles
        allocate(tables%u_beak(FOS_BEAK_SCAN_POINTS))
        do i = 1_ik, FOS_BEAK_SCAN_POINTS
            u_raw = -1.0_rk + 2.0_rk * real(i - 1_ik, rk) &
                    / real(FOS_BEAK_SCAN_POINTS - 1_ik, rk)
            tables%u_beak(i) = max(-U_BEAK_CLAMP, min(U_BEAK_CLAMP, u_raw))
        end do

        allocate(tables%bk_cos_even(FOS_BEAK_SCAN_POINTS, k_max))
        allocate(tables%bk_sin_odd(FOS_BEAK_SCAN_POINTS, k_max))
        do k = 1_ik, k_max
            do i = 1_ik, FOS_BEAK_SCAN_POINTS
                tables%bk_cos_even(i, k) = cos(tables%omega(k) * tables%u_beak(i))
                tables%bk_sin_odd(i, k) = sin(tables%psi(k) * tables%u_beak(i))
            end do
        end do

        ! Theta nodes (absent for the theta-less tier-1 forms)
        if (n_theta > 0_ik) then
            allocate(tables%thetas(n_theta), tables%cos_thetas(n_theta))
            do i = 1_ik, n_theta
                tables%thetas(i) = thetas(i)
                tables%cos_thetas(i) = cos(thetas(i))
            end do
        end if

        tables%initialized = .true.
        status = SHAPE_VALID

    end subroutine build_tables_s

end module fos_parameterization_workers_mod
