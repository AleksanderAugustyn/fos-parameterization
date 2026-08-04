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
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_TOO_MANY_PARAMS

    implicit none

    private

    ! Types
    public :: tables_t
    public :: fos_bundle_t

    ! Procedures
    public :: tables_init_s
    public :: tables_free_s
    public :: build_tables_s
    public :: compute_a2_s
    public :: compute_z_shift_s
    public :: compute_f_grid_s
    public :: beak_scan_f_min_s
    public :: scale_rho_grid_s
    public :: newton_radius_s

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

    !---------------------------------------------------------------------------
    ! Kernel constants
    !---------------------------------------------------------------------------
    ! MIGRATION NOTE: every constant below duplicates one in
    ! fos_parameterization_mod. The workers module sits UNDER the main module in
    ! the dependency order and cannot import from it, and the 1.x surface must
    ! stay bit-for-bit untouched until it is retired. The duplicates disappear
    ! with that surface; this module is their permanent home. Values must match
    ! their 1.x twin exactly — the kernels are parity-tested against it.

    !> Largest Fourier coefficient index the evaluator reads.
    integer(kind = ik), parameter, public :: FOS_MAX_K = 50_ik

    !> Coefficient pairs below this magnitude contribute nothing and are skipped
    !! — skipping is load-bearing for bitwise parity with the 1.x evaluator.
    real(kind = rk), parameter :: COEFF_NEGLIGIBLE = 1.0e-30_rk

    !> Smallest elongation the evaluator treats as a shape.
    real(kind = rk), parameter, public :: C_MIN = 1.0e-10_rk

    !> Beak-singularity threshold: shapes whose f_min falls below it are rejected.
    real(kind = rk), parameter, public :: F_MIN_THRESHOLD = 1.0e-3_rk

    !> Tip detection tolerance: f(±1) = 0 analytically, so u within roundoff of
    !! a tip is treated AS the tip (rho = 0, drho/dz = 0).
    real(kind = rk), parameter :: U_TIP_TOL = 4.0_rk * epsilon(1.0_rk)

    !> Interior nodes at or below this rho are a pinched (invalid) shape.
    real(kind = rk), parameter :: RHO_TOLERANCE = 1.0e-12_rk

    !> Newton residual tolerance, scaled by max(1, r).
    real(kind = rk), parameter :: NR_TOLERANCE = 1.0e-12_rk
    integer(kind = ik), parameter :: NR_MAX_ITER = 50_ik

    !> |cos(theta)| above which the radius is taken from the analytic pole.
    real(kind = rk), parameter :: POLE_THRESH = 1.0_rk - 1.0e-10_rk

    !> Newton bracket floor and the derivative magnitude below which a Newton
    !! step is abandoned for bisection.
    real(kind = rk), parameter :: R_LO_FLOOR = 1.0e-10_rk
    real(kind = rk), parameter :: DF_DR_FLOOR = 1.0e-14_rk

    !> Degenerate or missing elongation. Duplicated from fos_parameterization_mod
    !! (see the migration note above).
    integer(kind = ik), parameter, public :: FOS_ERROR_INVALID_C = 102_ik

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

    !> Scalar evaluator bundle: parameters, resolved z-shift, and the analytic
    !! Newton bracket bound. Exists so `newton_radius_s` can be elemental —
    !! Fortran requires every elemental dummy to be scalar, which an
    !! assumed-shape params(:) can never be.
    !!
    !! `r_hi_bound` is the caller's analytic outer bracket:
    !! `2*sqrt(rho_max**2 + max(z_max, abs(z_min))**2)` from the resolved rho(z)
    !! grid. It replaces the 1.x doubling loop, which silently returned the
    !! initial guess for extreme-oblate shapes whose equatorial radius exceeded
    !! 256x the polar extent.
    type :: fos_bundle_t
        integer(kind = ik) :: n_params = 0_ik
        real(kind = rk)    :: params(FOS_MAX_K) = 0.0_rk
        real(kind = rk)    :: z_shift = 0.0_rk
        real(kind = rk)    :: r_hi_bound = 0.0_rk
    end type fos_bundle_t

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

    !===========================================================================
    ! COEFFICIENTS
    !===========================================================================

    !> Computes a2 from the volume constraint, reporting oversize vectors.
    !!
    !! c plays no part, so any length from 0 to FOS_MAX_K is valid; an empty
    !! vector is the sphere, a2 = 0.
    !!
    !! @param[in]  params  FoS parameters, at most FOS_MAX_K entries
    !! @param[out] a2      Volume-constraint coefficient (0 on rejection)
    !! @param[out] status  SHAPE_VALID, or SHAPE_ERROR_TOO_MANY_PARAMS
    pure subroutine compute_a2_s(params, a2, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: a2
        integer(kind = ik), intent(out) :: status

        a2 = 0.0_rk

        if (size(params, kind = ik) > FOS_MAX_K) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        a2 = a2_f(params)
        status = SHAPE_VALID

    end subroutine compute_a2_s

    !> Computes the intrinsic z-shift that places the COM at the origin.
    !!
    !! @param[in]  params   FoS parameters; params(1) = c is mandatory
    !! @param[out] z_shift  Intrinsic shift (0 on rejection)
    !! @param[out] status   SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, or
    !!                      FOS_ERROR_INVALID_C (empty vector or c <= C_MIN)
    pure subroutine compute_z_shift_s(params, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: c, sum_term, a_odd, sign_factor
        integer(kind = ik) :: n

        z_shift = 0.0_rk

        if (size(params, kind = ik) > FOS_MAX_K) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        if (size(params, kind = ik) < 1_ik) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        c = params(1)
        if (c <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        sum_term = 0.0_rk
        do n = 1_ik, FOS_MAX_K
            a_odd = coefficient_f(params, 2_ik * n + 1_ik, 0.0_rk)
            if (abs(a_odd) < COEFF_NEGLIGIBLE) cycle
            sign_factor = merge(-1.0_rk, 1.0_rk, mod(n, 2_ik) == 0_ik)
            sum_term = sum_term + sign_factor * a_odd / real(n, rk)
        end do

        z_shift = (3.0_rk / (2.0_rk * PI_C)) * c * sum_term
        status = SHAPE_VALID

    end subroutine compute_z_shift_s

    !===========================================================================
    ! TABLED SHAPE FUNCTION
    !===========================================================================

    !> Evaluates f(u) and f'(u) at every u-grid node from the trig tables.
    !!
    !! f depends on the Fourier coefficients only — never on c — so the grid
    !! survives a pure elongation step.
    !!
    !! Preconditions (caller's, unchecked): `tables` initialized,
    !! `size(f_grid) = size(fp_grid) >= tables%n_points`, and the table's k_max
    !! covers the vector, k_max >= (size(params) + 2)/2 + 1.
    !!
    !! @param[in]  tables   Initialized trig tables
    !! @param[in]  params   FoS parameters (params(1) = c is not read)
    !! @param[out] f_grid   f at tables%u
    !! @param[out] fp_grid  df/du at tables%u
    subroutine compute_f_grid_s(tables, params, f_grid, fp_grid)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: f_grid(:)
        real(kind = rk), intent(out) :: fp_grid(:)

        real(kind = rk) :: a_even(tables%k_max), a_odd(tables%k_max)
        logical :: active(tables%k_max)
        real(kind = rk) :: sum_f, sum_fp, u
        integer(kind = ik) :: i, k

        call pair_coefficients_s(params, tables%k_max, a_even, a_odd, active)

        do i = 1_ik, tables%n_points
            sum_f = 0.0_rk
            sum_fp = 0.0_rk
            do k = 1_ik, tables%k_max
                if (.not. active(k)) cycle
                sum_f = sum_f + a_even(k) * tables%cos_even(i, k) &
                        + a_odd(k) * tables%sin_odd(i, k)
                sum_fp = sum_fp - a_even(k) * tables%omega(k) * tables%sin_even(i, k) &
                        + a_odd(k) * tables%psi(k) * tables%cos_odd(i, k)
            end do
            u = tables%u(i)
            f_grid(i) = 1.0_rk - u**2 - sum_f
            fp_grid(i) = -2.0_rk * u - sum_fp
        end do

    end subroutine compute_f_grid_s

    !> Scans f(u) over the clamped beak grid and applies the beak threshold.
    !!
    !! f -> 0 in the interior is a cusp: drho/dz diverges, Newton stops
    !! converging, and surface/Coulomb integrals blow up. The scan grid is
    !! clamped away from the poles, where f = 0 by construction.
    !!
    !! @param[in]  tables   Initialized trig tables
    !! @param[in]  params   FoS parameters (params(1) = c is not read)
    !! @param[out] f_min    Smallest f over the scan grid
    !! @param[out] beak_ok  .true. iff f_min > F_MIN_THRESHOLD
    subroutine beak_scan_f_min_s(tables, params, f_min, beak_ok)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: f_min
        logical, intent(out) :: beak_ok

        real(kind = rk) :: a_even(tables%k_max), a_odd(tables%k_max)
        logical :: active(tables%k_max)
        real(kind = rk) :: sum_f, f_val, u
        integer(kind = ik) :: i, k

        call pair_coefficients_s(params, tables%k_max, a_even, a_odd, active)

        f_min = huge(1.0_rk)
        do i = 1_ik, FOS_BEAK_SCAN_POINTS
            sum_f = 0.0_rk
            do k = 1_ik, tables%k_max
                if (.not. active(k)) cycle
                sum_f = sum_f + a_even(k) * tables%bk_cos_even(i, k) &
                        + a_odd(k) * tables%bk_sin_odd(i, k)
            end do
            u = tables%u_beak(i)
            f_val = 1.0_rk - u**2 - sum_f
            if (f_val < f_min) f_min = f_val
        end do

        beak_ok = f_min > F_MIN_THRESHOLD

    end subroutine beak_scan_f_min_s

    !===========================================================================
    ! CYLINDRICAL GRID
    !===========================================================================

    !> Scales a tabled f-grid into the cylindrical rho(z) grid for elongation c.
    !!
    !! z(i) = c*u(i) + z_shift_intrinsic, rho = sqrt(f/c), drho/dz =
    !! f'/(2c*sqrt(c*f)). The tips (|u| = 1 up to roundoff) are rho = 0,
    !! drho/dz = 0 by convention: f(±1) = 0 analytically, and the roundoff
    !! residue would otherwise be amplified into a ~1e7 slope.
    !!
    !! Preconditions (caller's, unchecked): c > C_MIN, `tables` initialized, and
    !! every array at least tables%n_points long.
    !!
    !! @param[in]  tables             Initialized trig tables
    !! @param[in]  c                  Elongation, params(1)
    !! @param[in]  z_shift_intrinsic  COM shift baked into z
    !! @param[in]  f_grid             f at tables%u
    !! @param[in]  fp_grid            df/du at tables%u
    !! @param[out] z                  Axial coordinate, shifted
    !! @param[out] rho                Cylindrical radius
    !! @param[out] drho_dz            Slope
    !! @param[out] rho_max            Largest rho on the grid (Newton bracket input)
    !! @param[out] rho_positive       .false. iff an interior node has rho <= RHO_TOLERANCE
    subroutine scale_rho_grid_s(tables, c, z_shift_intrinsic, f_grid, fp_grid, &
            z, rho, drho_dz, rho_max, rho_positive)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: c
        real(kind = rk), intent(in) :: z_shift_intrinsic
        real(kind = rk), intent(in) :: f_grid(:)
        real(kind = rk), intent(in) :: fp_grid(:)
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: rho_max
        logical, intent(out) :: rho_positive

        real(kind = rk) :: c_inv, u, rho_sq, sqrt_cf
        integer(kind = ik) :: i

        c_inv = 1.0_rk / c
        rho_max = 0.0_rk
        rho_positive = .true.

        do i = 1_ik, tables%n_points
            u = tables%u(i)
            z(i) = c * u + z_shift_intrinsic

            if (abs(u) >= 1.0_rk - U_TIP_TOL) then
                rho(i) = 0.0_rk
                drho_dz(i) = 0.0_rk
                cycle
            end if

            rho_sq = f_grid(i) * c_inv

            if (rho_sq > 0.0_rk) then
                rho(i) = sqrt(rho_sq)
                sqrt_cf = sqrt(c * f_grid(i))
                drho_dz(i) = fp_grid(i) / (2.0_rk * c * sqrt_cf)
            else
                rho(i) = 0.0_rk
                drho_dz(i) = 0.0_rk
            end if

            ! The positivity verdict covers exactly the 1.x interior range: the
            ! tips are the only nodes taking the branch above, and they are
            ! always i = 1 and i = n_points on the symmetric u-grid.
            if (rho(i) > rho_max) rho_max = rho(i)
            if (rho(i) <= RHO_TOLERANCE) rho_positive = .false.
        end do

    end subroutine scale_rho_grid_s

    !===========================================================================
    ! NEWTON RADIUS CORE
    !===========================================================================

    !> Elemental R(theta) and dR/dtheta at x = cos(theta) for a resolved shape.
    !!
    !! Solves F(r) = r*sin(theta) - rho(r*cos(theta)) = 0 with a
    !! bisection-safeguarded Newton-Raphson. For a star-convex shape the ray
    !! crosses the surface exactly once, F < 0 inside and F > 0 outside, so the
    !! bracket [R_LO_FLOOR, bundle%r_hi_bound] always contains the root — the
    !! bound is analytic, computed by the caller from the resolved rho(z) grid.
    !! Newton steps leaving the bracket are replaced by bisection, which keeps
    !! steep polar lobes out of a limit cycle.
    !!
    !! The derivative is implicit differentiation at the root:
    !!   dR/dtheta = -(r cos + drho_dz r sin) / (sin - drho_dz cos)
    !! Poles are analytic (converged, dR/dtheta = 0). Degenerate bundles are the
    !! caller's problem: validity is established before a bundle is built.
    !!
    !! Preconditions (caller's, unchecked): c > C_MIN, the shape star-convex
    !! about the shifted origin, and `r_hi_bound` outside the surface at every
    !! theta — a zero bound collapses the bracket and every radius with it.
    !!
    !! @param[in]  bundle     Resolved shape (params, z_shift, bracket bound)
    !! @param[in]  x          cos(theta), theta in [0, pi]
    !! @param[out] r          Radius at theta
    !! @param[out] dr_dtheta  dR/dtheta at theta
    !! @param[out] converged  .true. iff the final residual met NR_TOLERANCE
    elemental subroutine newton_radius_s(bundle, x, r, dr_dtheta, converged)

        type(fos_bundle_t), intent(in)  :: bundle
        real(kind = rk),    intent(in)  :: x
        real(kind = rk),    intent(out) :: r
        real(kind = rk),    intent(out) :: dr_dtheta
        logical,            intent(out) :: converged

        real(kind = rk) :: c, sin_theta, cos_theta
        real(kind = rk) :: z_max, z_min, r_north, r_south
        real(kind = rk) :: rho, drho_dz, z
        real(kind = rk) :: r_lo, r_hi, r_curr, r_new, delta_r, F_val, dF_dr
        integer(kind = ik) :: iter

        dr_dtheta = 0.0_rk

        c = bundle%params(1)

        cos_theta = x
        sin_theta = sqrt(max(1.0_rk - x**2, 0.0_rk))

        ! In the shifted frame the shape spans z in [-c + z_shift, c + z_shift]
        z_max = c + bundle%z_shift
        z_min = -c + bundle%z_shift
        r_north = z_max
        r_south = abs(z_min)

        ! Poles are analytic
        if (x > POLE_THRESH) then
            r = r_north
            converged = .true.
            return
        end if

        if (x < -POLE_THRESH) then
            r = r_south
            converged = .true.
            return
        end if

        ! Analytic bracket: F(r_lo) < 0 (origin inside the body),
        ! F(r_hi) > 0 (beyond the surface, by construction of r_hi_bound).
        r_lo = R_LO_FLOOR
        r_hi = bundle%r_hi_bound

        ! Initial guess
        r_curr = 0.5_rk * ((1.0_rk + x) * r_north + (1.0_rk - x) * r_south)
        r_curr = min(max(r_curr, 0.01_rk), r_hi)

        do iter = 1_ik, NR_MAX_ITER
            z = r_curr * cos_theta
            call bundle_rho_s(bundle, z, rho, drho_dz)

            F_val = r_curr * sin_theta - rho
            dF_dr = sin_theta - drho_dz * cos_theta

            ! Maintain the sign-change bracket
            if (F_val < 0.0_rk) then
                r_lo = r_curr
            else
                r_hi = r_curr
            end if

            ! Residual-based convergence: |F| is the geometric distance between
            ! the trial point and the surface, which is what callers care about.
            if (abs(F_val) < NR_TOLERANCE * max(1.0_rk, r_curr)) exit

            if (abs(dF_dr) > DF_DR_FLOOR) then
                delta_r = F_val / dF_dr
                r_new = r_curr - delta_r
            else
                r_new = r_lo - 1.0_rk  ! force bisection
            end if

            ! Newton step leaving the bracket -> bisect instead
            if (r_new <= r_lo .or. r_new >= r_hi) then
                r_new = 0.5_rk * (r_lo + r_hi)
            end if
            r_curr = r_new
        end do

        r = r_curr

        ! Recompute rho and drho_dz at the final r so the convergence verdict and
        ! the implicit-differentiation inputs both match the returned radius.
        ! This also covers the max-iterations exit, where the loop-carried values
        ! lag one iterate.
        z = r * cos_theta
        call bundle_rho_s(bundle, z, rho, drho_dz)

        F_val = r * sin_theta - rho
        converged = abs(F_val) < NR_TOLERANCE * max(1.0_rk, r)

        dF_dr = sin_theta - drho_dz * cos_theta
        if (abs(dF_dr) > DF_DR_FLOOR) then
            dr_dtheta = -(r * cos_theta + drho_dz * r * sin_theta) / dF_dr
        else
            ! Vertical tangent — excluded for star-convex shapes by the
            ! conversion margin; return 0 rather than a garbage slope.
            dr_dtheta = 0.0_rk
        end if

    end subroutine newton_radius_s

    !===========================================================================
    ! PRIVATE HELPERS
    !===========================================================================

    !> a2 from the volume constraint, without the length check.
    pure function a2_f(params) result(a2)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: a2

        integer(kind = ik) :: n, idx, n_params
        real(kind = rk) :: a_2n, sign_factor

        a2 = 0.0_rk
        n_params = size(params, kind = ik)

        do n = 2_ik, FOS_MAX_K
            idx = 2_ik * n - 1_ik
            if (idx > n_params) exit
            a_2n = params(idx)
            if (abs(a_2n) < COEFF_NEGLIGIBLE) cycle
            sign_factor = merge(1.0_rk, -1.0_rk, mod(n, 2_ik) == 0_ik)
            a2 = a2 + sign_factor * a_2n / real(2_ik * n - 1_ik, rk)
        end do

    end function a2_f

    !> Fourier coefficient a_k, with a2 supplied by the caller (it costs a full
    !! volume-constraint sum, so callers that need many coefficients hoist it).
    pure function coefficient_f(params, k, a2) result(a_k)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: k
        real(kind = rk), intent(in) :: a2
        real(kind = rk) :: a_k

        integer(kind = ik) :: idx

        if (k < 2_ik) then
            a_k = 0.0_rk
        else if (k == 2_ik) then
            a_k = a2
        else
            idx = k - 1_ik
            if (idx <= size(params, kind = ik)) then
                a_k = params(idx)
            else
                a_k = 0.0_rk
            end if
        end if

    end function coefficient_f

    !> Splits a parameter vector into the (even, odd) coefficient pair per
    !! Fourier order, flagging the pairs that contribute.
    !!
    !! `active` reproduces the 1.x skip rule verbatim — a pair with both
    !! coefficients below COEFF_NEGLIGIBLE is dropped from the sum, not added as
    !! zero, so tabled sums are bitwise-identical to the live evaluator's.
    pure subroutine pair_coefficients_s(params, k_max, a_even, a_odd, active)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: k_max
        real(kind = rk), intent(out) :: a_even(:)
        real(kind = rk), intent(out) :: a_odd(:)
        logical, intent(out) :: active(:)

        real(kind = rk) :: a2
        integer(kind = ik) :: k

        a2 = a2_f(params)

        do k = 1_ik, k_max
            a_even(k) = coefficient_f(params, 2_ik * k, a2)
            a_odd(k) = coefficient_f(params, 2_ik * k + 1_ik, a2)
            active(k) = abs(a_even(k)) >= COEFF_NEGLIGIBLE &
                    .or. abs(a_odd(k)) >= COEFF_NEGLIGIBLE
        end do

    end subroutine pair_coefficients_s

    !> Live (untabled) f(u) and f'(u) — the Newton path evaluates f at arbitrary
    !! u, off every grid node.
    pure subroutine eval_f_s(params, u, f, fp)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: u
        real(kind = rk), intent(out) :: f
        real(kind = rk), intent(out) :: fp

        integer(kind = ik) :: k, k_max, n_params
        real(kind = rk) :: a2, omega_k, psi_k, a_even, a_odd
        real(kind = rk) :: cos_even, sin_even, cos_odd, sin_odd
        real(kind = rk) :: sum_f, sum_fp

        a2 = a2_f(params)
        n_params = size(params, kind = ik)
        k_max = min((n_params + 2_ik) / 2_ik + 1_ik, FOS_MAX_K)

        sum_f = 0.0_rk
        sum_fp = 0.0_rk

        do k = 1_ik, k_max
            a_even = coefficient_f(params, 2_ik * k, a2)
            a_odd = coefficient_f(params, 2_ik * k + 1_ik, a2)
            if (abs(a_even) < COEFF_NEGLIGIBLE .and. abs(a_odd) < COEFF_NEGLIGIBLE) cycle

            omega_k = real(2_ik * k - 1_ik, rk) * PI_C / 2.0_rk
            psi_k = real(k, rk) * PI_C

            cos_even = cos(omega_k * u)
            sin_even = sin(omega_k * u)
            cos_odd = cos(psi_k * u)
            sin_odd = sin(psi_k * u)

            sum_f = sum_f + a_even * cos_even + a_odd * sin_odd
            sum_fp = sum_fp - a_even * omega_k * sin_even + a_odd * psi_k * cos_odd
        end do

        f = 1.0_rk - u**2 - sum_f
        fp = -2.0_rk * u - sum_fp

    end subroutine eval_f_s

    !> rho and drho/dz at an axial coordinate in the shifted frame.
    !!
    !! Same tip and degenerate-c conventions as the 1.x point evaluator: outside
    !! the shape, at a tip, or on a non-positive f, the surface is rho = 0 with
    !! zero slope.
    pure subroutine bundle_rho_s(bundle, z, rho, drho_dz)

        type(fos_bundle_t), intent(in) :: bundle
        real(kind = rk), intent(in) :: z
        real(kind = rk), intent(out) :: rho
        real(kind = rk), intent(out) :: drho_dz

        real(kind = rk) :: c, c_inv, u, f_val, fp_val, sqrt_cf

        rho = 0.0_rk
        drho_dz = 0.0_rk

        if (bundle%n_params < 1_ik) return
        c = bundle%params(1)
        if (c <= C_MIN) return

        c_inv = 1.0_rk / c
        u = (z - bundle%z_shift) * c_inv
        if (abs(u) >= 1.0_rk - U_TIP_TOL) return

        call eval_f_s(bundle%params(1:bundle%n_params), u, f_val, fp_val)

        if (f_val > 0.0_rk) then
            sqrt_cf = sqrt(c * f_val)
            rho = sqrt(f_val * c_inv)
            drho_dz = fp_val / (2.0_rk * c * sqrt_cf)
        end if

    end subroutine bundle_rho_s

end module fos_parameterization_workers_mod
