!> Module for Fourier-over-Spheroid (FoS) nuclear shape parameterization.
!!
!! ## Architecture
!!
!! This module is responsible for:
!! 1. Computing the FoS shape in cylindrical coordinates (ρ, z)
!! 2. Validating that ρ ≤ 0 only at the poles (valid shape)
!! 3. Validating f_min threshold to prevent "beak" singularities
!! 4. Computing the intrinsic z-shift and baking it into the grid
!! 5. Checking star-convexity and finding additional z-shift if needed
!! 6. Converting the ρ(z) representation to R(θ) for the radius grid
!!
!! Physics policy is NOT this module's job. Rejections here are mathematics
!! (c > 0, rho > 0) plus this module's own numerical constraints (beak f_min,
!! star-convexity margin, both required by the R(theta) conversion). Filtering
!! of mathematically valid but physically meaningless shapes (e.g. necked
!! shapes at low elongation) belongs to the consumer: see
!! radius_grid_mod's validate_physics.
!!
!! ## Defense in Depth: Beak Singularity Detection
!!
!! The FoS parameterization can develop "beak" cusps when f(u) approaches zero
!! in the interior. While mathematically valid (f > 0), these shapes cause:
!! - Derivative singularities: dρ/dz → ∞ as f(u) → 0
!! - Newton-Raphson convergence failures
!! - Catastrophic errors in surface/Coulomb integrals
!!
!! We detect and reject such shapes early by checking f_min > F_MIN_THRESHOLD.
!!
!! ## Coordinate Systems and Shifts
!!
!! The FoS parameterization defines shapes in terms of a parameter u ∈ [-1, 1]:
!!   - z = c × u (in reduced units where R0 = 1)
!!   - ρ² = f(u) / c
!!
!! The intrinsic z-shift (z_shift_intrinsic) places the center of mass at the origin.
!! For R(θ) computation, we need the shape to be star-convex from the origin.
!!
!! **Key convention:** The z_shift returned by this module represents the total
!! shift applied to the shape's z-coordinates. Positive z_shift means the shape
!! is shifted in the +z direction. When checking star-convexity for a shape
!! shifted by z_shift, we check: (z + z_shift) × dρ/dz - ρ ≤ 0
!!
!! ## Workflow
!!
!! ```
!! 1. compute_rho_z_grid_s()     - Create internal ρ(z) grid with n_grid points
!! 2. validate_rho_grid_s()      - Check ρ > 0 for interior points AND f_min threshold
!! 3. compute_fos_z_shift_f()    - Calculate intrinsic z_shift
!! 4. apply_z_shift_to_grid_s()  - Shift grid to place COM at origin
!! 5. check_star_convexity_s()   - Verify star-convexity, find additional shift if needed
!! 6. compute_radius_grid_s()    - Convert to R(θ) grid
!! ```
!!
module fos_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C

    implicit none

    private

    ! Main entry point
    public :: compute_fos_radius_grid_s

    ! Per-shape resolve step (validity + z_shift + analytic pole radii)
    public :: compute_fos_shape_s

    ! Batch R(theta) + dR/dtheta evaluation at caller-supplied thetas
    public :: compute_fos_radius_and_derivative_at_thetas_s

    ! Individual workflow components
    public :: compute_rho_z_grid_s
    public :: validate_rho_grid_s
    public :: check_star_convexity_s
    public :: compute_radius_at_theta_s

    ! Helper functions
    public :: compute_fos_a2_f
    public :: get_fos_coefficient_f
    public :: compute_fos_f_and_derivatives_s
    public :: compute_fos_z_shift_f
    public :: compute_rho_at_z_s
    public :: compute_radius_fos_with_zshift_s
    public :: compute_fos_neck_s

    ! Shape/evaluation split: scalar shape bundle + elemental R, dR/dtheta core
    public :: make_fos_shape_f
    public :: compute_fos_radius_and_derivative_s

    ! Internal rho(z) grid type — public so consumers and tests can drive the
    ! grid-level workflow (compute_rho_z_grid_s / validate_rho_grid_s /
    ! check_star_convexity_s) directly.
    public :: rho_z_grid_t

    ! Module constants
    integer(kind = ik), parameter :: MAX_K = 50_ik
    real(kind = rk), parameter :: NR_TOLERANCE = 1.0e-12_rk
    integer(kind = ik), parameter :: NR_MAX_ITER = 50_ik
    real(kind = rk), parameter :: POLE_THRESH = 1.0_rk - 1.0e-10_rk
    real(kind = rk), parameter, public :: C_MIN = 1.0e-10_rk
    real(kind = rk), parameter :: COEFF_NEGLIGIBLE = 1.0e-30_rk
    ! Tip detection tolerance: grid construction reaches u = ±1 only up to
    ! roundoff (~1 ulp), and f(±1) = 0 analytically, so drho/dz = f'/(2c*sqrt(cf))
    ! amplifies that residue to ~1e7 unless u this close to a tip is treated AS
    ! the tip (rho = 0, drho/dz = 0). No interior evaluation comes within 1e-15.
    real(kind = rk), parameter :: U_TIP_TOL = 4.0_rk * epsilon(1.0_rk)

    !---------------------------------------------------------------------------
    ! Defense in Depth: Beak Singularity Threshold
    !---------------------------------------------------------------------------
    !! Minimum f(u) threshold - shapes with f_min below this are rejected.
    !! This prevents numerical instabilities from "beak" singularities where
    !! the shape approaches but doesn't quite reach self-intersection.
    !!
    !! Physical interpretation: f_min < 1e-3 corresponds to shapes with
    !! extremely concentrated deformations that nuclear matter cannot support.
    !!
    !! Tuning guidance:
    !!   - 1e-2: Very conservative, rejects moderately deformed shapes
    !!   - 1e-3: Recommended default, balances safety and coverage
    !!   - 1e-4: Aggressive, allows shapes closer to the validity boundary
    real(kind = rk), parameter, public :: F_MIN_THRESHOLD = 1.0e-3_rk

    !---------------------------------------------------------------------------
    ! Number of sampling points for f_min search
    !---------------------------------------------------------------------------
    integer(kind = ik), parameter :: F_MIN_SAMPLE_POINTS = 1001_ik

    ! Validation error codes
    integer(kind = ik), parameter, public :: FOS_VALID = 0_ik
    integer(kind = ik), parameter, public :: FOS_ERROR_RHO_NEGATIVE = 1_ik
    integer(kind = ik), parameter, public :: FOS_ERROR_NOT_STAR_CONVEX = 2_ik
    integer(kind = ik), parameter, public :: FOS_ERROR_INVALID_C = 3_ik
    integer(kind = ik), parameter, public :: FOS_ERROR_BEAK_SINGULARITY = 4_ik

    !---------------------------------------------------------------------------
    ! Star-convexity safety margin
    !---------------------------------------------------------------------------
    !! The star-convexity condition (z + z_shift)·dρ/dz - ρ ≤ 0 admits shapes
    !! whose surface grazes a ray from the origin (condition ≈ 0). Such shapes
    !! produce a near-vertical wall in R(θ) (|dR/dθ| ~ 1/margin) that no
    !! practical θ grid can resolve, corrupting surface and Coulomb integrals.
    !! Requiring (z + z_shift)·dρ/dz - ρ ≤ -margin rejects shapes that the
    !! R(θ) representation cannot carry accurately. In reduced units (R0 = 1).
    real(kind = rk), parameter, public :: STAR_CONVEXITY_MARGIN = 1.0e-1_rk

    !> Scalar bundle of FoS parameters plus a resolved z-shift. Exists so the
    !! radius evaluator can be elemental: Fortran requires every elemental
    !! dummy to be scalar, which an assumed-shape params(:) can never be.
    type, public :: fos_shape_t
        integer(kind = ik) :: n_params = 0_ik
        real(kind = rk)    :: params(MAX_K) = 0.0_rk
        real(kind = rk)    :: z_shift = 0.0_rk
    end type fos_shape_t

    ! Internal ρ(z) grid type
    type :: rho_z_grid_t
        integer(kind = ik) :: n_points = 0_ik
        real(kind = rk), allocatable :: z(:)
        real(kind = rk), allocatable :: rho(:)
        real(kind = rk), allocatable :: drho_dz(:)
        real(kind = rk) :: z_min = 0.0_rk
        real(kind = rk) :: z_max = 0.0_rk
        logical :: initialized = .false.
    end type rho_z_grid_t

contains

    !===========================================================================
    ! MAIN ENTRY POINT
    !===========================================================================

    !> Validates FoS shape and computes radius grid for radius_grid_mod.
    subroutine compute_fos_radius_grid_s(params, n_grid, radii, z_shift, is_valid, message, &
            n_rho_grid, error_code)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_grid
        real(kind = rk), intent(out) :: radii(n_grid)
        real(kind = rk), intent(out) :: z_shift
        logical, intent(out) :: is_valid
        character(len = *), intent(out) :: message
        integer(kind = ik), intent(in), optional :: n_rho_grid
        integer(kind = ik), intent(out), optional :: error_code

        integer(kind = ik) :: n_internal
        real(kind = rk) :: r_north, r_south
        integer(kind = ik) :: err_local

        radii = 0.0_rk
        z_shift = 0.0_rk
        is_valid = .false.
        message = ''

        n_internal = n_grid
        if (present(n_rho_grid)) n_internal = max(100_ik, n_rho_grid)

        ! Steps 1-4 live in compute_fos_shape_s; codes and messages unchanged.
        call compute_fos_shape_s(params, n_internal, z_shift, r_north, r_south, &
                is_valid, message, err_local)
        if (.not. is_valid) then
            if (present(error_code)) error_code = err_local
            return
        end if

        ! Step 5: Compute R(θ) grid using the final z_shift
        call compute_radius_grid_internal_s(params, n_grid, z_shift, radii)

        is_valid = .true.
        message = ''
        if (present(error_code)) error_code = FOS_VALID

    end subroutine compute_fos_radius_grid_s

    !> Per-shape resolve step: validity + total z-shift + analytic pole radii,
    !! computed once on the internal rho(z) grid — independent of any theta
    !! node set. Feed the resulting z_shift to
    !! compute_fos_radius_and_derivative_at_thetas_s for any number of node sets.
    !!
    !! On any failure all outputs are zero-filled and is_valid is .false.
    !!
    !! @param[in]  params      FoS parameters: params(1) = c, params(k-1) = a_k for k >= 3
    !! @param[in]  n_rho_grid  Internal rho(z) grid size, used verbatim
    !! @param[out] z_shift     Total shift (intrinsic COM + star-convexity search)
    !! @param[out] r_north     R(0) = c + z_shift (pole extent in the shifted frame)
    !! @param[out] r_south     R(pi) = |-c + z_shift|
    !! @param[out] is_valid    .true. iff the shape passed all mathematical checks
    !! @param[out] message     Empty on success
    !! @param[out] error_code  FOS_VALID on success
    subroutine compute_fos_shape_s(params, n_rho_grid, z_shift, r_north, r_south, &
            is_valid, message, error_code)
        real(kind = rk),    intent(in)  :: params(:)
        integer(kind = ik), intent(in)  :: n_rho_grid
        real(kind = rk),    intent(out) :: z_shift
        real(kind = rk),    intent(out) :: r_north
        real(kind = rk),    intent(out) :: r_south
        logical,            intent(out) :: is_valid
        character(len = *), intent(out) :: message
        integer(kind = ik), intent(out) :: error_code

        type(rho_z_grid_t) :: rho_grid
        real(kind = rk)    :: z_shift_intrinsic
        logical            :: rho_valid, star_convex

        z_shift    = 0.0_rk
        r_north    = 0.0_rk
        r_south    = 0.0_rk
        is_valid   = .false.
        message    = ''
        error_code = FOS_VALID

        ! Step 1: internal rho(z) grid (unshifted, z in [-c, c])
        call compute_rho_z_grid_s(params, n_rho_grid, rho_grid, error_code, message)
        if (error_code /= FOS_VALID) then
            call deallocate_rho_grid(rho_grid)
            return
        end if

        ! Step 2: interior rho > 0 AND f_min above the beak threshold
        call validate_rho_grid_s(rho_grid, params, rho_valid, error_code, message)
        if (.not. rho_valid) then
            call deallocate_rho_grid(rho_grid)
            return
        end if

        ! Step 3 + 3b: intrinsic COM shift, baked into the grid
        z_shift_intrinsic = compute_fos_z_shift_f(params)
        call apply_z_shift_to_grid_s(rho_grid, z_shift_intrinsic)

        ! Step 4: star-convexity (may find an additional shift)
        call check_star_convexity_s(params, rho_grid, z_shift, star_convex, message)
        if (.not. star_convex) then
            error_code = FOS_ERROR_NOT_STAR_CONVEX
            z_shift    = 0.0_rk
            call deallocate_rho_grid(rho_grid)
            return
        end if

        z_shift = z_shift_intrinsic + z_shift

        ! Analytic pole radii in the shifted frame — the same expressions the
        ! Newton evaluator's pole branch uses.
        r_north = params(1) + z_shift
        r_south = abs(-params(1) + z_shift)

        is_valid   = .true.
        message    = ''
        error_code = FOS_VALID
        call deallocate_rho_grid(rho_grid)

    end subroutine compute_fos_shape_s

    !> Batch-evaluate R(theta) and dR/dtheta at caller-supplied thetas, for a
    !! shape whose validity and z_shift were already established by
    !! compute_fos_shape_s. One scalar shape build, one elemental sweep — the
    !! N Newton solves are independent and free to vectorize.
    !!
    !! Caller contract: size(radii) == size(dr_dthetas) == size(thetas).
    !! Degenerate params (empty, c <= C_MIN) yield the unit-sphere fallback
    !! r = 1, dr_dtheta = 0 — same library guarantee as the scalar evaluator.
    pure subroutine compute_fos_radius_and_derivative_at_thetas_s(params, thetas, &
            z_shift, radii, dr_dthetas)
        real(kind = rk), intent(in)  :: params(:)
        real(kind = rk), intent(in)  :: thetas(:)
        real(kind = rk), intent(in)  :: z_shift
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dthetas(:)

        type(fos_shape_t) :: shape

        shape = make_fos_shape_f(params, z_shift)
        call compute_fos_radius_and_derivative_s(shape, cos(thetas), radii, dr_dthetas)

    end subroutine compute_fos_radius_and_derivative_at_thetas_s

    !===========================================================================
    ! STEP 1: CREATE ρ(z) GRID
    !===========================================================================

    subroutine compute_rho_z_grid_s(params, n_points, grid, error_code, message)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        type(rho_z_grid_t), intent(out) :: grid
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        real(kind = rk) :: c, c_inv, u, dz, f_val, fp_val, sqrt_cf, rho_sq
        integer(kind = ik) :: i

        error_code = FOS_VALID
        message = ''
        grid%initialized = .false.

        if (size(params) < 1) then
            error_code = FOS_ERROR_INVALID_C
            message = 'Empty parameter array'
            return
        end if

        c = params(1)
        if (c <= C_MIN) then
            error_code = FOS_ERROR_INVALID_C
            write(message, '(A,ES12.4)') 'Elongation c must be positive, got: ', c
            return
        end if

        c_inv = 1.0_rk / c

        grid%n_points = n_points
        allocate(grid%z(n_points), grid%rho(n_points), grid%drho_dz(n_points))

        grid%z_min = -c
        grid%z_max = c
        dz = 2.0_rk * c / real(n_points - 1_ik, rk)

        do i = 1_ik, n_points
            grid%z(i) = -c + real(i - 1_ik, rk) * dz
            u = grid%z(i) * c_inv

            ! Same tip convention as compute_rho_at_z_s: at |u| = 1 (up to
            ! roundoff) the shape ends, rho = 0 and drho/dz = 0 by convention.
            if (abs(u) >= 1.0_rk - U_TIP_TOL) then
                grid%rho(i) = 0.0_rk
                grid%drho_dz(i) = 0.0_rk
                cycle
            end if

            call compute_fos_f_and_derivatives_s(params, u, f_val, fp_val)
            rho_sq = f_val * c_inv

            if (rho_sq > 0.0_rk) then
                grid%rho(i) = sqrt(rho_sq)
                sqrt_cf = sqrt(c * f_val)
                grid%drho_dz(i) = fp_val / (2.0_rk * c * sqrt_cf)
            else
                grid%rho(i) = 0.0_rk
                grid%drho_dz(i) = 0.0_rk
            end if
        end do

        grid%initialized = .true.

    end subroutine compute_rho_z_grid_s

    !> Applies z-shift to the ρ(z) grid coordinates.
    !!
    !! Shifts all z-coordinates: z_new = z_old + z_shift
    !! Places center of mass at origin when z_shift = z_shift_intrinsic.
    subroutine apply_z_shift_to_grid_s(grid, z_shift)
        type(rho_z_grid_t), intent(inout) :: grid
        real(kind = rk), intent(in) :: z_shift

        if (.not. grid%initialized) return

        grid%z = grid%z + z_shift
        grid%z_min = grid%z_min + z_shift
        grid%z_max = grid%z_max + z_shift

    end subroutine apply_z_shift_to_grid_s

    !===========================================================================
    ! STEP 2: VALIDATE ρ(z) GRID (Enhanced with beak detection)
    !===========================================================================

    !> Validates the ρ(z) grid for physical consistency.
    !!
    !! This enhanced validation includes:
    !! 1. Interior point check: ρ > 0 for all interior points
    !! 2. Beak singularity detection: f_min > F_MIN_THRESHOLD
    !!
    !! The beak detection prevents numerical catastrophes from shapes that are
    !! mathematically valid but have f(u) approaching zero, causing derivative
    !! singularities that corrupt surface and Coulomb integrals.
    subroutine validate_rho_grid_s(grid, params, is_valid, error_code, message)
        type(rho_z_grid_t), intent(in) :: grid
        real(kind = rk), intent(in) :: params(:)
        logical, intent(out) :: is_valid
        integer(kind = ik), intent(out) :: error_code
        character(len = *), intent(out) :: message

        integer(kind = ik) :: i
        real(kind = rk), parameter :: RHO_TOLERANCE = 1.0e-12_rk
        real(kind = rk) :: f_min, u, f_val, u_at_f_min

        is_valid = .true.
        error_code = FOS_VALID
        message = ''

        if (.not. grid%initialized) then
            is_valid = .false.
            error_code = FOS_ERROR_RHO_NEGATIVE
            message = 'rho(z) grid not initialized'
            return
        end if

        !-----------------------------------------------------------------------
        ! Check 1: Interior points must have ρ > 0
        !-----------------------------------------------------------------------
        do i = 2_ik, grid%n_points - 1_ik
            if (grid%rho(i) <= RHO_TOLERANCE) then
                is_valid = .false.
                error_code = FOS_ERROR_RHO_NEGATIVE
                write(message, '(A,ES12.4,A,ES12.4)') &
                        'Invalid shape: rho <= 0 at z = ', grid%z(i), ', rho = ', grid%rho(i)
                return
            end if
        end do

        !-----------------------------------------------------------------------
        ! Check 2: Beak singularity detection via f_min threshold
        !-----------------------------------------------------------------------
        ! Sample f(u) at many points to find the minimum value.
        ! Even if the grid rho values are positive, f(u) approaching zero
        ! causes derivative singularities that corrupt calculations.
        f_min = huge(1.0_rk)
        u_at_f_min = 0.0_rk

        do i = 1_ik, F_MIN_SAMPLE_POINTS
            u = -1.0_rk + 2.0_rk * real(i - 1_ik, rk) / real(F_MIN_SAMPLE_POINTS - 1_ik, rk)
            ! Avoid exact poles where f=0 is expected
            u = max(-0.999_rk, min(0.999_rk, u))

            call compute_fos_f_and_derivatives_s(params, u, f_val)

            if (f_val < f_min) then
                f_min = f_val
                u_at_f_min = u
            end if
        end do

        if (f_min < F_MIN_THRESHOLD) then
            is_valid = .false.
            error_code = FOS_ERROR_BEAK_SINGULARITY
            write(message, '(A,ES10.3,A,ES10.3,A,F6.3)') &
                    'Shape too close to beak singularity: f_min = ', f_min, &
                    ' < threshold = ', F_MIN_THRESHOLD, &
                    ' at u = ', u_at_f_min
            return
        end if

    end subroutine validate_rho_grid_s

    !===========================================================================
    ! STEP 4: CHECK STAR-CONVEXITY
    !===========================================================================

    !> Checks if shape is star-convex from origin and finds valid shift if not.
    !!
    !! A shape is star-convex from the origin if every surface point can be
    !! connected to the origin by a straight line that doesn't intersect the
    !! surface elsewhere. For a shape already shifted (grid z-coordinates
    !! shifted), we check with z_shift=0 first, then search for additional shift.
    !!
    !! Returns z_shift as the ADDITIONAL shift needed (beyond what's already
    !! baked into the grid). If star-convex at z_shift=0, returns z_shift=0.
    subroutine check_star_convexity_s(params, rho_grid, z_shift, is_star_convex, message)
        real(kind = rk), intent(in) :: params(:)
        type(rho_z_grid_t), intent(in) :: rho_grid
        real(kind = rk), intent(out) :: z_shift
        logical, intent(out) :: is_star_convex
        character(len = *), intent(out) :: message

        z_shift = 0.0_rk
        is_star_convex = .false.
        message = ''

        ! First try z_shift=0 (grid already shifted to place COM at origin)
        if (is_star_convex_from_grid_f(rho_grid, 0.0_rk)) then
            is_star_convex = .true.
            z_shift = 0.0_rk
            return
        end if

        ! If not star-convex at z_shift=0, search for an additional shift
        call find_star_convex_shift_from_grid_s(params, rho_grid, z_shift, is_star_convex)

        if (.not. is_star_convex) then
            message = 'Shape is not star-convex and no valid z-shift found'
            z_shift = 0.0_rk
        end if

    end subroutine check_star_convexity_s

    !> Checks star-convexity for a given z-shift.
    !!
    !! The condition for star-convexity from the origin when the shape is
    !! shifted by z_shift is: (z + z_shift) × dρ/dz - ρ ≤ 0
    !!
    !! **IMPORTANT FIX:** The z_shift is ADDED (not subtracted) because positive
    !! z_shift means the shape is shifted in the +z direction, moving surface
    !! points from z_old to z_new = z_old + z_shift.
    pure function is_star_convex_from_grid_f(grid, z_shift) result(is_convex)
        type(rho_z_grid_t), intent(in) :: grid
        real(kind = rk), intent(in) :: z_shift
        logical :: is_convex

        integer(kind = ik) :: i
        real(kind = rk) :: z_shifted, test_val

        is_convex = .true.
        if (.not. grid%initialized) then
            is_convex = .false.
            return
        end if

        ! Check star-convexity condition at each interior point
        ! For a shape shifted by z_shift: check (z + z_shift) × dρ/dz - ρ ≤ -margin
        ! The margin rejects grazing-ray shapes that R(θ) cannot represent.
        do i = 2_ik, grid%n_points - 1_ik
            z_shifted = grid%z(i) + z_shift  ! ADD z_shift, not subtract!
            test_val = z_shifted * grid%drho_dz(i) - grid%rho(i)
            if (test_val > -STAR_CONVEXITY_MARGIN) then
                is_convex = .false.
                return
            end if
        end do

    end function is_star_convex_from_grid_f

    !> Maximum star-convexity test value over interior points for a given shift.
    !!
    !! g(s) = max_i [ (z_i + s) * drho_dz_i - rho_i ]. This is a pointwise max of
    !! affine functions of s, hence convex in s. The shape is star-convex at shift
    !! s iff g(s) <= -STAR_CONVEXITY_MARGIN (see is_star_convex_from_grid_f).
    pure function max_star_convexity_value_f(grid, z_shift) result(g)
        type(rho_z_grid_t), intent(in) :: grid
        real(kind = rk), intent(in) :: z_shift
        real(kind = rk) :: g

        integer(kind = ik) :: i
        real(kind = rk) :: test_val

        if (.not. grid%initialized) then
            g = huge(1.0_rk)
            return
        end if

        g = -huge(1.0_rk)
        do i = 2_ik, grid%n_points - 1_ik
            test_val = (grid%z(i) + z_shift) * grid%drho_dz(i) - grid%rho(i)
            if (test_val > g) g = test_val
        end do
    end function max_star_convexity_value_f

    !> Searches for a z-shift that makes the shape star-convex.
    !!
    !! The grid is assumed to be already shifted by z_shift_intrinsic (COM at origin).
    !! This function searches for an ADDITIONAL shift beyond that.
    !! The z_shift=0 case is already tested in check_star_convexity_s.
    !!
    !! g(s) = max_i[(z_i + s) drho_dz_i - rho_i] is convex piecewise-linear in s
    !! (pointwise max of affine functions), so golden-section finds its exact
    !! global minimum with no local-minimum traps. This resolves the true optimum
    !! shift, unlike the former neck-centered guesses + 0.1c-step fallback, which
    !! stepped over narrow acceptance windows near the star-convexity margin.
    subroutine find_star_convex_shift_from_grid_s(params, grid, z_shift, is_convex)
        real(kind = rk), intent(in) :: params(:)
        type(rho_z_grid_t), intent(in) :: grid
        real(kind = rk), intent(out) :: z_shift
        logical, intent(out) :: is_convex

        ! (sqrt(5) - 1) / 2; SHIFT_TOL and bracket are in reduced units (R0 = 1).
        real(kind = rk), parameter :: GOLDEN = 0.6180339887498949_rk
        real(kind = rk), parameter :: SHIFT_TOL = 1.0e-6_rk
        integer(kind = ik), parameter :: MAX_ITER = 200_ik

        real(kind = rk) :: c, a, b, x1, x2, f1, f2
        integer(kind = ik) :: it

        z_shift = 0.0_rk
        is_convex = .false.
        c = params(1)

        ! Any origin that can be star-convex lies inside the body, which is
        ! contained in +/-2c on the COM-shifted grid, so this brackets the minimum.
        a = -2.0_rk * c
        b = 2.0_rk * c
        x1 = b - GOLDEN * (b - a)
        x2 = a + GOLDEN * (b - a)
        f1 = max_star_convexity_value_f(grid, x1)
        f2 = max_star_convexity_value_f(grid, x2)

        do it = 1_ik, MAX_ITER
            if (b - a <= SHIFT_TOL) exit
            if (f1 < f2) then
                b = x2
                x2 = x1
                f2 = f1
                x1 = b - GOLDEN * (b - a)
                f1 = max_star_convexity_value_f(grid, x1)
            else
                a = x1
                x1 = x2
                f1 = f2
                x2 = a + GOLDEN * (b - a)
                f2 = max_star_convexity_value_f(grid, x2)
            end if
        end do

        z_shift = 0.5_rk * (a + b)
        is_convex = is_star_convex_from_grid_f(grid, z_shift)
        if (.not. is_convex) z_shift = 0.0_rk
    end subroutine find_star_convex_shift_from_grid_s

    !> Finds the z-position of the neck (minimum ρ between two maxima).
    subroutine find_neck_from_grid_s(grid, z_neck, found)
        type(rho_z_grid_t), intent(in) :: grid
        real(kind = rk), intent(out) :: z_neck
        logical, intent(out) :: found

        integer(kind = ik) :: i, n_maxima, max1_idx, max2_idx
        integer(kind = ik) :: left_idx, right_idx, neck_idx, j
        real(kind = rk) :: max1_rho, max2_rho, min_rho

        z_neck = 0.0_rk
        found = .false.
        if (.not. grid%initialized) return

        n_maxima = 0_ik
        max1_idx = 0_ik
        max2_idx = 0_ik
        max1_rho = -1.0_rk
        max2_rho = -1.0_rk

        do i = 2_ik, grid%n_points - 1_ik
            if (grid%rho(i) > grid%rho(i - 1) .and. grid%rho(i) > grid%rho(i + 1)) then
                n_maxima = n_maxima + 1_ik
                if (grid%rho(i) > max1_rho) then
                    max2_rho = max1_rho
                    max2_idx = max1_idx
                    max1_rho = grid%rho(i)
                    max1_idx = i
                else if (grid%rho(i) > max2_rho) then
                    max2_rho = grid%rho(i)
                    max2_idx = i
                end if
            end if
        end do

        if (n_maxima < 2_ik .or. max1_idx == 0_ik .or. max2_idx == 0_ik) return

        if (max1_idx < max2_idx) then
            left_idx = max1_idx
            right_idx = max2_idx
        else
            left_idx = max2_idx
            right_idx = max1_idx
        end if

        min_rho = huge(1.0_rk)
        neck_idx = left_idx

        do j = left_idx, right_idx
            if (grid%rho(j) < min_rho) then
                min_rho = grid%rho(j)
                neck_idx = j
            end if
        end do

        z_neck = grid%z(neck_idx)
        found = .true.

    end subroutine find_neck_from_grid_s

    !> Finds the neck (interior minimum of ρ between two maxima) of a FoS shape.
    !!
    !! A coarse grid scan brackets the neck, then Newton iteration on f'(u) = 0
    !! refines it to machine precision (the neck is a minimum of f, so f'' > 0
    !! there). Used to define the scission line: scission when rho_neck → 0.
    !!
    !! @param[in]  params    FoS parameters [c, a3, a4, ...]
    !! @param[out] z_neck    Neck z-position in the intrinsic frame (COM at origin)
    !! @param[out] rho_neck  Neck radius in reduced units (R0 = 1)
    !! @param[out] found     .false. if the shape has no neck or the grid is invalid
    !! @param[in]  n_grid    Optional scan resolution (default 1001)
    subroutine compute_fos_neck_s(params, z_neck, rho_neck, found, n_grid)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(in), optional :: n_grid

        integer(kind = ik), parameter :: NECK_NEWTON_MAX_ITER = 50_ik
        real(kind = rk), parameter :: NECK_NEWTON_TOL = 1.0e-14_rk

        type(rho_z_grid_t) :: grid
        integer(kind = ik) :: n_points, err, iter
        character(len = 256) :: message
        real(kind = rk) :: c, u, u_lo, u_hi, du, step
        real(kind = rk) :: f_val, fp_val, fpp_val, z_neck_grid

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        n_points = 1001_ik
        if (present(n_grid)) n_points = max(101_ik, n_grid)

        call compute_rho_z_grid_s(params, n_points, grid, err, message)
        if (err /= FOS_VALID) then
            call deallocate_rho_grid(grid)
            return
        end if

        ! Grid scan on the unshifted grid (z ∈ [-c, c])
        call find_neck_from_grid_s(grid, z_neck_grid, found)
        call deallocate_rho_grid(grid)
        if (.not. found) return

        ! Newton refinement of the f-minimum, bracketed to one grid spacing
        ! around the scan result so it cannot escape to a different extremum.
        c = params(1)
        du = 2.0_rk / real(n_points - 1_ik, rk)
        u = z_neck_grid / c
        u_lo = max(-1.0_rk, u - du)
        u_hi = min(1.0_rk, u + du)

        do iter = 1_ik, NECK_NEWTON_MAX_ITER
            call compute_fos_f_and_derivatives_s(params, u, f_val, fp_val, fpp_val)
            if (fpp_val <= 0.0_rk) exit
            step = fp_val / fpp_val
            u = min(u_hi, max(u_lo, u - step))
            if (abs(step) < NECK_NEWTON_TOL) exit
        end do

        call compute_fos_f_and_derivatives_s(params, u, f_val)
        rho_neck = sqrt(max(f_val, 0.0_rk) / c)
        z_neck = c * u + compute_fos_z_shift_f(params)

    end subroutine compute_fos_neck_s

    !===========================================================================
    ! STEP 5: COMPUTE R(θ) GRID
    !===========================================================================

    subroutine compute_radius_grid_internal_s(params, n_grid, z_shift, radii)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_grid
        real(kind = rk), intent(in) :: z_shift
        real(kind = rk), intent(out) :: radii(n_grid)

        integer(kind = ik) :: i
        real(kind = rk) :: h, theta, x, r

        h = PI_C / real(n_grid - 1_ik, rk)

        do i = 1_ik, n_grid
            theta = real(i - 1_ik, rk) * h
            x = cos(theta)
            call compute_radius_fos_with_zshift_s(params, x, z_shift, r)
            radii(i) = r
        end do

    end subroutine compute_radius_grid_internal_s

    subroutine compute_radius_at_theta_s(params, theta, z_shift, r)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: theta
        real(kind = rk), intent(in) :: z_shift
        real(kind = rk), intent(out) :: r

        call compute_radius_fos_with_zshift_s(params, cos(theta), z_shift, r)

    end subroutine compute_radius_at_theta_s

    !===========================================================================
    ! HELPER FUNCTIONS
    !===========================================================================

    subroutine deallocate_rho_grid(grid)
        type(rho_z_grid_t), intent(inout) :: grid
        if (allocated(grid%z)) deallocate(grid%z)
        if (allocated(grid%rho)) deallocate(grid%rho)
        if (allocated(grid%drho_dz)) deallocate(grid%drho_dz)
        grid%initialized = .false.
        grid%n_points = 0_ik
    end subroutine deallocate_rho_grid

    !> Computes ρ and optionally dρ/dz at a given z-coordinate.
    !!
    !! The z argument is in the shifted frame (where origin is for R(θ)).
    !! To get the FoS parameter u: u = (z - z_shift) / c
    pure subroutine compute_rho_at_z_s(params, z, z_shift, rho, drho_dz)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: z
        real(kind = rk), intent(in) :: z_shift
        real(kind = rk), intent(out) :: rho
        real(kind = rk), intent(out), optional :: drho_dz

        real(kind = rk) :: c, c_inv, u, f_val, fp_val, sqrt_cf

        rho = 0.0_rk
        if (present(drho_dz)) drho_dz = 0.0_rk

        if (size(params) < 1) return
        c = params(1)
        if (c <= C_MIN) return

        c_inv = 1.0_rk / c
        u = (z - z_shift) * c_inv
        if (abs(u) >= 1.0_rk - U_TIP_TOL) return

        if (present(drho_dz)) then
            call compute_fos_f_and_derivatives_s(params, u, f_val, fp_val)
        else
            call compute_fos_f_and_derivatives_s(params, u, f_val)
        end if

        if (f_val > 0.0_rk) then
            sqrt_cf = sqrt(c * f_val)
            rho = sqrt(f_val * c_inv)
            if (present(drho_dz)) drho_dz = fp_val / (2.0_rk * c * sqrt_cf)
        end if

    end subroutine compute_rho_at_z_s

    !===========================================================================
    ! CORE FoS FUNCTIONS
    !===========================================================================

    !> Computes a2 from the volume constraint formula.
    pure function compute_fos_a2_f(params) result(a2)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: a2
        integer(kind = ik) :: n, idx, n_params
        real(kind = rk) :: a_2n, sign_factor

        a2 = 0.0_rk
        n_params = size(params, kind = ik)

        do n = 2_ik, MAX_K
            idx = 2_ik * n - 1_ik
            if (idx > n_params) exit
            a_2n = params(idx)
            if (abs(a_2n) < COEFF_NEGLIGIBLE) cycle
            sign_factor = merge(1.0_rk, -1.0_rk, mod(n, 2_ik) == 0_ik)
            a2 = a2 + sign_factor * a_2n / real(2_ik * n - 1_ik, rk)
        end do

    end function compute_fos_a2_f

    !> Gets FoS coefficient a_k from parameter array.
    pure function get_fos_coefficient_f(params, k) result(a_k)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: k
        real(kind = rk) :: a_k
        integer(kind = ik) :: idx

        if (k < 2_ik) then
            a_k = 0.0_rk
        else if (k == 2_ik) then
            a_k = compute_fos_a2_f(params)
        else
            idx = k - 1_ik
            if (idx <= size(params, kind = ik)) then
                a_k = params(idx)
            else
                a_k = 0.0_rk
            end if
        end if

    end function get_fos_coefficient_f

    !> Computes f(u) and its derivatives.
    pure subroutine compute_fos_f_and_derivatives_s(params, u, f, fp, fpp)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: u
        real(kind = rk), intent(out) :: f
        real(kind = rk), intent(out), optional :: fp, fpp

        integer(kind = ik) :: k, k_max, n_params
        real(kind = rk) :: omega_k, psi_k, a_even, a_odd
        real(kind = rk) :: cos_even, sin_even, cos_odd, sin_odd
        real(kind = rk) :: sum_f, sum_fp, sum_fpp
        logical :: need_fp, need_fpp

        need_fp = present(fp)
        need_fpp = present(fpp)
        n_params = size(params, kind = ik)
        k_max = min((n_params + 2_ik) / 2_ik + 1_ik, MAX_K)

        sum_f = 0.0_rk
        sum_fp = 0.0_rk
        sum_fpp = 0.0_rk

        do k = 1_ik, k_max
            a_even = get_fos_coefficient_f(params, 2_ik * k)
            a_odd = get_fos_coefficient_f(params, 2_ik * k + 1_ik)
            if (abs(a_even) < COEFF_NEGLIGIBLE .and. abs(a_odd) < COEFF_NEGLIGIBLE) cycle

            omega_k = real(2_ik * k - 1_ik, rk) * PI_C / 2.0_rk
            psi_k = real(k, rk) * PI_C

            cos_even = cos(omega_k * u)
            sin_even = sin(omega_k * u)
            cos_odd = cos(psi_k * u)
            sin_odd = sin(psi_k * u)

            sum_f = sum_f + a_even * cos_even + a_odd * sin_odd

            if (need_fp .or. need_fpp) then
                sum_fp = sum_fp - a_even * omega_k * sin_even + a_odd * psi_k * cos_odd
            end if

            if (need_fpp) then
                sum_fpp = sum_fpp - a_even * omega_k**2 * cos_even - a_odd * psi_k**2 * sin_odd
            end if
        end do

        f = 1.0_rk - u**2 - sum_f
        if (need_fp) fp = -2.0_rk * u - sum_fp
        if (need_fpp) fpp = -2.0_rk - sum_fpp

    end subroutine compute_fos_f_and_derivatives_s

    !> Computes the intrinsic z-shift to place COM at origin.
    pure function compute_fos_z_shift_f(params) result(z_sh)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk) :: z_sh
        real(kind = rk) :: c, sum_term, a_odd, sign_factor
        integer(kind = ik) :: n

        z_sh = 0.0_rk
        if (size(params) < 1) return
        c = params(1)
        if (c <= C_MIN) return

        sum_term = 0.0_rk
        do n = 1_ik, MAX_K
            a_odd = get_fos_coefficient_f(params, 2_ik * n + 1_ik)
            if (abs(a_odd) < COEFF_NEGLIGIBLE) cycle
            sign_factor = merge(-1.0_rk, 1.0_rk, mod(n, 2_ik) == 0_ik)
            sum_term = sum_term + sign_factor * a_odd / real(n, rk)
        end do

        z_sh = (3.0_rk / (2.0_rk * PI_C)) * c * sum_term

    end function compute_fos_z_shift_f

    !> Computes R(θ) at x = cos(θ). Thin wrapper over the elemental core
    !! compute_fos_radius_and_derivative_s — same Newton path, derivative
    !! discarded. Kept for its established callers.
    pure subroutine compute_radius_fos_with_zshift_s(params, x, z_shift, r)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in), value :: x
        real(kind = rk), intent(in), value :: z_shift
        real(kind = rk), intent(out) :: r

        type(fos_shape_t) :: shape
        real(kind = rk)   :: dr_unused

        shape = make_fos_shape_f(params, z_shift)
        call compute_fos_radius_and_derivative_s(shape, x, r, dr_unused)

    end subroutine compute_radius_fos_with_zshift_s

    !> Bundle a params vector and resolved z_shift into the scalar shape type.
    !!
    !! Entries beyond MAX_K are dropped: radius evaluation reads coefficients
    !! only up to k = MAX_K (compute_fos_f_and_derivatives_s caps there, so
    !! params indices above MAX_K - 1 are never touched), so truncation cannot
    !! change any evaluated radius. z_shift arrives already resolved — it may
    !! depend on higher entries, which is why compute_fos_shape_s takes the
    !! full vector, not this type.
    pure function make_fos_shape_f(params, z_shift) result(shape)
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: z_shift
        type(fos_shape_t) :: shape

        integer(kind = ik) :: n

        n = min(size(params, kind = ik), MAX_K)
        shape%n_params    = n
        shape%params(:)   = 0.0_rk
        shape%params(1:n) = params(1:n)
        shape%z_shift     = z_shift
    end function make_fos_shape_f

    !> Elemental core: R(theta) and dR/dtheta at x = cos(theta), theta in [0, pi],
    !! for a resolved shape.
    !!
    !! Solves F(r) = r×sin(θ) - ρ(r×cos(θ)) = 0 with bisection-safeguarded
    !! Newton-Raphson. For a star-convex shape the ray crosses the surface
    !! exactly once, with F < 0 inside the body and F > 0 outside, so a
    !! sign-change bracket [r_lo, r_hi] always exists. Newton steps that would
    !! leave the bracket are replaced by bisection, which guarantees convergence
    !! even for steep polar lobes where plain Newton enters a limit cycle
    !! (stepping outside the body collapses the derivative to sin(θ) and the
    !! resulting jump overshoots).
    !!
    !! The derivative is implicit differentiation at the root:
    !!   dR/dtheta = -(r cos + drho_dz r sin) / (sin - drho_dz cos)
    !! Pole branches and the degenerate unit-sphere fallback return dr_dtheta = 0
    !! (smooth axisymmetric surface has zero slope in R(theta) at the poles).
    elemental subroutine compute_fos_radius_and_derivative_s(shape, x, r, dr_dtheta)
        type(fos_shape_t), intent(in)  :: shape
        real(kind = rk),   intent(in)  :: x
        real(kind = rk),   intent(out) :: r
        real(kind = rk),   intent(out) :: dr_dtheta

        real(kind = rk) :: c, sin_theta, cos_theta
        real(kind = rk) :: z_max, z_min, r_north, r_south
        real(kind = rk) :: rho, drho_dz, z
        real(kind = rk) :: r_lo, r_hi, r_curr, r_new, delta_r, F_val, dF_dr
        integer(kind = ik) :: iter

        dr_dtheta = 0.0_rk

        if (shape%n_params < 1_ik) then
            r = 1.0_rk
            return
        end if

        c = shape%params(1)
        if (c <= C_MIN) then
            r = 1.0_rk
            return
        end if

        cos_theta = x
        sin_theta = sqrt(max(1.0_rk - x**2, 0.0_rk))

        ! In shifted frame, shape spans z ∈ [-c + z_shift, c + z_shift]
        z_max = c + shape%z_shift
        z_min = -c + shape%z_shift
        r_north = z_max
        r_south = abs(z_min)

        ! Handle poles analytically
        if (x > POLE_THRESH) then
            r = r_north
            return
        end if

        if (x < -POLE_THRESH) then
            r = r_south
            return
        end if

        ! Bracket the root: F(r_lo) < 0 (origin inside the body),
        ! F(r_hi) > 0 (beyond the surface). Expand r_hi if needed (very
        ! oblate shapes have equatorial radii exceeding the polar extents).
        r_lo = 1.0e-10_rk
        r_hi = 2.0_rk * max(r_north, r_south)
        do iter = 1_ik, 8_ik
            z = r_hi * cos_theta
            call compute_rho_at_z_s(shape%params(1:shape%n_params), z, shape%z_shift, rho)
            if (r_hi * sin_theta - rho > 0.0_rk) exit
            r_hi = 2.0_rk * r_hi
        end do

        ! Initial guess
        r_curr = 0.5_rk * ((1.0_rk + x) * r_north + (1.0_rk - x) * r_south)
        r_curr = min(max(r_curr, 0.01_rk), r_hi)

        ! Safeguarded Newton: solve F(r) = r×sin(θ) - ρ(r×cos(θ)) = 0
        do iter = 1_ik, NR_MAX_ITER
            z = r_curr * cos_theta
            call compute_rho_at_z_s(shape%params(1:shape%n_params), z, shape%z_shift, rho, drho_dz)

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

            if (abs(dF_dr) > 1.0e-14_rk) then
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

        ! Recompute rho and drho_dz at the final r so the implicit-differentiation
        ! inputs match the returned radius exactly. This also covers the
        ! max-iterations exit, where the loop-carried values lag one iterate.
        z = r * cos_theta
        call compute_rho_at_z_s(shape%params(1:shape%n_params), z, shape%z_shift, &
                rho, drho_dz)
        dF_dr = sin_theta - drho_dz * cos_theta
        if (abs(dF_dr) > 1.0e-14_rk) then
            dr_dtheta = -(r * cos_theta + drho_dz * r * sin_theta) / dF_dr
        else
            ! Vertical tangent — excluded for star-convex shapes by the
            ! conversion margin; return 0 rather than a garbage slope.
            dr_dtheta = 0.0_rk
        end if

    end subroutine compute_fos_radius_and_derivative_s

end module fos_parameterization_mod