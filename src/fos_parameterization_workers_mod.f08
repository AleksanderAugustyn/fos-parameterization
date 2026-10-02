!> Arithmetic of the Fourier-over-Spheroid (FoS) parameterization.
!!
!! This module holds every floating-point kernel of the library and nothing
!! else: the trig tables (`tables_t`) and their build, the stage kernels (a2,
!! intrinsic z-shift, f grid, beak scan, rho grid, origin resolve), the Newton
!! bundle, the R(theta) solve loop and the neck refinement. No tier logic, no
!! `iso_c_binding`, no derived-type I/O.
!!
!! The tiers live in `fos_parameterization_mod`: `cache_t`, the usage checks,
!! the parameter trimming, the per-call scratch, the gates and both the cached
!! and the one-shot routines. That module calls the kernels here and performs
!! no floating-point arithmetic on an output value.
!!
!! ## Why the split is a correctness rule (normative)
!!
!! One-shot == cached, short == zero-padded and repeat == first call are
!! BITWISE claims. Two tier paths that produce the same output call the same
!! entry points of this module with the same data. The library's sources are
!! compiled without LTO, so a procedure of this module cannot be inlined into
!! the tier module: both paths execute the same machine code, on every
!! compiler. Inside one file the compiler still inlines and specializes —
!! 2.0.0 measured a one-ulp split in dR/dtheta from two inlined copies of the
!! solve loop — which is why no tier routine may live here, and why the tier
!! module must not do arithmetic of its own.
!!
!! `tables_t` is internal: the public surface wraps it in `cache_t` and does
!! not export it. Its components stay public so the kernel-level tests can
!! read them.
module fos_parameterization_workers_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_WRONG_PARAM_COUNT

    implicit none

    private

    ! Types
    public :: tables_t
    public :: fos_bundle_t

    ! Tables
    public :: build_tables_s
    public :: tables_free_s

    ! Parameter-vector helpers
    public :: active_length_f
    public :: table_order_f
    public :: shifted_origin_f

    ! Stage kernels, in pipeline order
    public :: compute_a2_s
    public :: compute_z_shift_s
    public :: compute_f_grid_s
    public :: beak_scan_f_min_s
    public :: scale_rho_grid_s
    public :: resolve_origin_s

    ! R(theta) conversion and neck
    public :: fos_bundle_f
    public :: newton_radius_s
    public :: solve_thetas_s
    public :: refine_neck_s

    !---------------------------------------------------------------------------
    ! Table-construction parameters
    !---------------------------------------------------------------------------
    !! Minimum u-grid resolution. Below this the rho(z) grid is too coarse for
    !! the neck/star-convexity scans that read it.
    integer(kind = ik), parameter, public :: FOS_N_POINTS_FLOOR = 100_ik

    !! Sample count of the beak (f_min) scan grid.
    integer(kind = ik), parameter, public :: FOS_BEAK_SCAN_POINTS = 1001_ik

    !! Beak-scan clamp. f(±1) = 0 analytically, so an unclamped scan would find
    !! f_min = 0 for every shape — including the sphere — and reject it all.
    real(kind = rk), parameter :: U_BEAK_CLAMP = 0.999_rk

    !---------------------------------------------------------------------------
    ! Kernel constants
    !---------------------------------------------------------------------------
    ! This module is the single owner of every kernel constant. The 1.x
    ! duplicates in fos_parameterization_mod were deleted with the 1.x surface;
    ! the two constants the raw evaluators there still need
    ! (FOS_COEFF_NEGLIGIBLE, FOS_U_TIP_TOL) are public and imported from here.

    !> Largest Fourier coefficient index the evaluator reads.
    integer(kind = ik), parameter, public :: FOS_MAX_K = 50_ik

    !> N_max: the longest parameter vector either tier accepts. Equal to
    !! FOS_MAX_K = 50 by construction — `fos_bundle_t` carries exactly that many
    !! slots, so a longer vector could not be evaluated even if it were
    !! accepted. A cache is built for any `max_params` up to this limit.
    integer(kind = ik), parameter, public :: FOS_MAX_PARAMS = FOS_MAX_K

    !> Coefficient pairs below this magnitude contribute nothing and are
    !! skipped — skipping is load-bearing for bitwise parity between the tabled
    !! kernels and the live evaluator. Public: the raw evaluators in
    !! fos_parameterization_mod share it.
    real(kind = rk), parameter, public :: FOS_COEFF_NEGLIGIBLE = 1.0e-30_rk

    !> Smallest elongation the evaluator treats as a shape.
    real(kind = rk), parameter, public :: C_MIN = 1.0e-10_rk

    !---------------------------------------------------------------------------
    ! Beak-singularity threshold
    !---------------------------------------------------------------------------
    !> A shape is accepted iff f_min > threshold, where f_min is the minimum of
    !! f(u) over the 1001-point scan clamped to |u| <= 0.999. f(+-1) = 0 for
    !! every shape, so near a tip f ~ |f'(+-1)| * 1e-3: a boundary-clamped
    !! f_min is a tip-slope proxy (sphere: f' = -2 -> baseline 1.999e-3, the
    !! hard upper bound for any threshold), and f' -> 0 makes the tip a beak.
    !! Two failure branches, both diverging as f_min^(-1/2): boundary-clamped
    !! minima (polar branch, drho/dz ~ f'/(2c sqrt(c f_min)) -> R' blow-up
    !! near a pole, prefactor c^(-3/2)) and strict-interior minima (neck
    !! branch, slope bounded but d2rho/dz2 ~ f''/(2c^2 sqrt(c f_min)) -> R''
    !! blow-up at the waist, prefactor c^(-5/2)).
    !!
    !! Value retuned 2026-08-11 (1e-3 -> 5e-4) from the geometry-sweep
    !! representability probe (fos_param_geometry_sweep_test --probe; f_min
    !! binned 6/decade, branch-split, GL-4096 V/S/round-trip through the
    !! beak-ungated diagnostic conversion) on the validated box
    !! c in [0.75, 3.50], a3 in [0, 0.60], a4 in [-0.20, 0.75],
    !! a5 in [-0.15, 0.15], a6 in [-0.10, 0.10], coarse (0.05^5) plus fine
    !! (a5/a6 at 0.01, 9.5M points). Measured structure, fine grid: below
    !! ~2.2e-4 accuracy is catastrophic (dS to 4e-2, round-trip breaches);
    !! the band [3.162e-4, 4.642e-4) carries a quadrature tail to dS 2.4e-4
    !! (high-c, a3 = 0.60 shapes whose resolved origin sits ~0.04 off the
    !! south tip); every accepted bin above 5e-4 meets dV/V <= 8.5e-8,
    !! dS/S <= 4.3e-5, round-trip <= 6.3e-12 over 7.5M shapes with zero
    !! conversion failures. 5e-4 (half the legacy 1e-3) clears both dirty
    !! bands and stays 4x below the sphere clamp bound; the residual
    !! in-domain worst (dS ~ 4.3e-5 at f_min ~ 3e-3) is a family tail no
    !! legal threshold can evict. The a6 bound is part of the statement: at
    !! |a6| = 0.15 an a4 = 0.75, a3 = 0.55, c <= 0.85 family exists whose
    !! broad interior near-neck (R'' ~ 5e4) degrades GL-4096 V/S at
    !! f_min ~ 3.9e-3 — healthy by any legal f_min threshold; catching it
    !! (and the high-c tails) needs a future curvature/origin-aware
    !! criterion. The scan clamp is load-bearing: f_min scales linearly with
    !! the clamp distance, so threshold and clamp (0.999) move together. In
    !! reduced units (R0 = 1).
    real(kind = rk), parameter, public :: F_MIN_THRESHOLD = 5.0e-4_rk

    !> Tip detection tolerance: f(±1) = 0 analytically, so u within roundoff of
    !! a tip is treated AS the tip (rho = 0, drho/dz = 0). Public: the raw
    !! evaluators in fos_parameterization_mod share it.
    real(kind = rk), parameter, public :: FOS_U_TIP_TOL = 4.0_rk * epsilon(1.0_rk)

    !> Interior nodes at or below this rho are a pinched (invalid) shape.
    real(kind = rk), parameter :: RHO_TOLERANCE = 1.0e-12_rk

    !> Newton residual tolerance, scaled by max(1, r).
    real(kind = rk), parameter :: NR_TOLERANCE = 1.0e-12_rk

    !> Newton-phase length and total iteration cap. Iterations 1 to
    !! NR_NEWTON_PHASE run the historical safeguarded-Newton update untouched
    !! (50 was the old hard cap, so every node that converged before this
    !! change still walks a bit-identical iterate path). A node still open
    !! after that is treated as Newton-resistant — e.g. an attracting
    !! period-2 orbit strictly inside the bracket, invisible to the
    !! leaves-the-bracket safeguard — and the remaining iterations bisect
    !! unconditionally: ~70 halvings collapse any bracket (r_hi_bound <= ~9)
    !! below the NR_BRACKET_ULPS acceptance width.
    integer(kind = ik), parameter :: NR_NEWTON_PHASE = 50_ik
    integer(kind = ik), parameter :: NR_MAX_ITER = 120_ik

    !> Bracket width, in ulps of the larger endpoint, at or below which the
    !! root is pinned to machine precision and the node is accepted even if
    !! the residual test is unsatisfiable (pole boundary layer: dF/dr ~ 1e4
    !! makes one ulp of r move F by more than NR_TOLERANCE).
    real(kind = rk), parameter :: NR_BRACKET_ULPS = 2.0_rk

    !> |cos(theta)| above which the radius is taken from the analytic pole.
    real(kind = rk), parameter :: POLE_THRESH = 1.0_rk - 1.0e-10_rk

    !> Newton bracket floor and the derivative magnitude below which a Newton
    !! step is abandoned for bisection.
    real(kind = rk), parameter :: R_LO_FLOOR = 1.0e-10_rk
    real(kind = rk), parameter :: DF_DR_FLOOR = 1.0e-14_rk

    !---------------------------------------------------------------------------
    ! Star-convexity safety margin
    !---------------------------------------------------------------------------
    !> A shape is accepted iff, at its best origin s* = argmin_s g(s),
    !! g(s*) = max_i[(z_i+s*)*drho/dz_i - rho_i] <= -margin. The mathematical
    !! single-valued (representable) boundary is g(s*) = 0, where R(theta)
    !! develops a vertical wall (|dR/dtheta| -> infinity); the margin holds
    !! shapes off that singularity.
    !!
    !! Value empirically retuned 2026-07-05 (0.1 -> 0.01) via the RewriteProject
    !! representability probe. At spectral GL-4096 density the FoS -> R(theta)
    !! conversion reproduces V/S to well within tolerance for EVERY
    !! single-valued shape (round-trip ~3e-12, dS <= 1.2e-6 even in the last bin
    !! before g=0), so the intrinsic representability limit is g=0 and the
    !! margin is only a safety buffer off that singularity. At 0.01 the
    !! production sweep keeps ~400x (volume) / 1e6x (surface) headroom under the
    !! tolerances; 0.01 is one bin off the g=0 singularity. Lowering it required
    !! fixing origin selection (see ORIGIN_CONDITION_MARGIN): the former s=0 fast
    !! path evaluated marginal-at-COM shapes from a steep origin, corrupting the
    !! volume integral. The energy model's neck-radius cutoff selects physical
    !! shapes downstream. In reduced units (R0 = 1).
    real(kind = rk), parameter, public :: STAR_CONVEXITY_MARGIN = 1.0e-2_rk

    !---------------------------------------------------------------------------
    ! COM-origin conditioning threshold for star-convexity origin selection
    !---------------------------------------------------------------------------
    !> The R(theta) origin is the COM (additional shift 0) when the shape is
    !! already well-conditioned there — g(0) = max_i[z_i*drho/dz_i - rho_i] <=
    !! -threshold — and otherwise the star-convexity optimum s* = argmin_s g(s),
    !! which is always at least as well-conditioned (g(s*) <= g(0)). 0.1 is the
    !! empirically-validated COM conditioning bound: at g(0) <= -0.1 the GL
    !! volume/surface integrals converge to ~1e-12 from the COM origin. Held
    !! fixed (independent of the acceptance margin above) so lowering the margin
    !! preserves the origin — and thus the z_shift — of every previously-accepted
    !! shape, while routing the newly-admitted marginal shapes to their
    !! best-conditioned origin.
    real(kind = rk), parameter :: ORIGIN_CONDITION_MARGIN = 1.0e-1_rk

    !> Golden-section constants for the origin search: (sqrt(5) - 1) / 2, the
    !! bracket half-width factor (the bracket is [-2c, 2c]), the stopping width,
    !! and the iteration cap. Values carried over from the 1.x minimizer, which
    !! this module replaced.
    real(kind = rk), parameter :: GOLDEN = 0.6180339887498949_rk
    real(kind = rk), parameter :: SHIFT_TOL = 1.0e-6_rk
    integer(kind = ik), parameter :: SHIFT_MAX_ITER = 200_ik

    !> Newton refinement of the neck: iteration cap and the step size below
    !! which the minimum of f is taken as located. Carried over from the 1.x
    !! neck scanner, which this module replaced.
    integer(kind = ik), parameter :: NECK_NEWTON_MAX_ITER = 50_ik
    real(kind = rk), parameter :: NECK_NEWTON_TOL = 1.0e-14_rk

    !> Degenerate or missing elongation. Owned here and re-exported by
    !! fos_parameterization_mod, which is where consumers meet it.
    integer(kind = ik), parameter, public :: FOS_ERROR_INVALID_C = 102_ik

    !> Interior node with rho <= 0. Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_RHO_NEGATIVE = 100_ik

    !> A Newton radius solve that did not meet NR_TOLERANCE at some node.
    !! Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_CONVERGENCE = 104_ik

    !> Output array whose size differs from the grid it must hold.
    !! Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_BUFFER_MISMATCH = 105_ik

    !> Shape not star-convex from any origin. Owned here, re-exported by the
    !! main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_NOT_STAR_CONVEX = 101_ik

    !> f_min below the beak threshold. Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_BEAK_SINGULARITY = 103_ik

    !> Parameter-independent trigonometric tables for FoS shape evaluation.
    !!
    !! Immutable after `build_tables_s`, therefore shareable across threads.
    !! Internal to the library: the public `cache_t` of
    !! `fos_parameterization_mod` holds one by value.
    !!
    !! The theta nodes are stored raw, WITHOUT a precomputed cosine: the bitwise
    !! contract between the fixed-grid and at-thetas forms requires ONE cos
    !! evaluation site (see `solve_thetas_s`), so both pass `thetas` and let the
    !! solver take the cosine.
    type, public :: tables_t
        integer(kind = ik) :: n_points = 0_ik, n_theta = 0_ik, k_max = 0_ik
        real(kind = rk), allocatable :: u(:)                              !! n_points
        real(kind = rk), allocatable :: cos_even(:, :), sin_even(:, :)    !! (n_points, k_max)
        real(kind = rk), allocatable :: cos_odd(:, :), sin_odd(:, :)      !! (n_points, k_max)
        real(kind = rk), allocatable :: bk_cos_even(:, :), bk_sin_odd(:, :) !! (FOS_BEAK_SCAN_POINTS, k_max)
        real(kind = rk), allocatable :: u_beak(:)                         !! FOS_BEAK_SCAN_POINTS
        real(kind = rk), allocatable :: thetas(:)                         !! n_theta
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

    !> Releases every table array and clears the initialized flag. Infallible.
    !!
    !! @param[in,out] tables  Tables to reset; safe on an already-free instance
    pure subroutine tables_free_s(tables)

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
        if (allocated(tables%omega)) deallocate(tables%omega)
        if (allocated(tables%psi)) deallocate(tables%psi)

        tables%n_points = 0_ik
        tables%n_theta = 0_ik
        tables%k_max = 0_ik
        tables%initialized = .false.

    end subroutine tables_free_s

    !> Table-construction worker: u-grid, Fourier bases, beak grid, theta nodes.
    !!
    !! Accepts an empty theta set — the theta-less one-shot forms (shape,
    !! rho_z_grid, neck, optimum) need no theta nodes, and then `thetas` stays
    !! unallocated with `n_theta = 0`. The public `cache_init_s` applies the
    !! one-node floor before it calls this.
    !!
    !! A grid the process cannot allocate is SHAPE_ERROR_INVALID_GRID, not an
    !! abort: `n_points` and `k_max` come from the caller, so a request the heap
    !! cannot satisfy is a rejected grid like any other. Every allocation is
    !! checked, and a failure releases whatever was already taken, leaving
    !! `tables` in the freed (uninitialized) state.
    !!
    !! @param[out] tables    Freshly built tables (reset on any rejection)
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  thetas    Polar angles in [0, pi]; may be empty
    !! @param[in]  k_max     Number of Fourier orders to tabulate, >= 1
    !! @param[out] status    SHAPE_VALID, or SHAPE_ERROR_INVALID_GRID
    pure subroutine build_tables_s(tables, n_points, thetas, k_max, status)

        type(tables_t), intent(out) :: tables
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: k_max
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i, k, n_theta
        integer :: alloc_stat
        real(kind = rk) :: u_raw

        if (n_points < FOS_N_POINTS_FLOOR) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        ! Zero Fourier orders would tabulate nothing and silently return the
        ! sphere for every parameter vector.
        if (k_max < 1_ik) then
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

        allocate(tables%omega(k_max), tables%psi(k_max), stat = alloc_stat)
        if (alloc_stat /= 0) then
            call reject_tables_s(tables, status)
            return
        end if
        do k = 1_ik, k_max
            tables%omega(k) = real(2_ik * k - 1_ik, rk) * PI_C / 2.0_rk
            tables%psi(k) = real(k, rk) * PI_C
        end do

        ! u-grid and its Fourier basis
        allocate(tables%u(n_points), stat = alloc_stat)
        if (alloc_stat /= 0) then
            call reject_tables_s(tables, status)
            return
        end if
        do i = 1_ik, n_points
            tables%u(i) = -1.0_rk + 2.0_rk * real(i - 1_ik, rk) / real(n_points - 1_ik, rk)
        end do

        allocate(tables%cos_even(n_points, k_max), tables%sin_even(n_points, k_max), &
                tables%cos_odd(n_points, k_max), tables%sin_odd(n_points, k_max), &
                stat = alloc_stat)
        if (alloc_stat /= 0) then
            call reject_tables_s(tables, status)
            return
        end if
        do k = 1_ik, k_max
            do i = 1_ik, n_points
                tables%cos_even(i, k) = cos(tables%omega(k) * tables%u(i))
                tables%sin_even(i, k) = sin(tables%omega(k) * tables%u(i))
                tables%cos_odd(i, k) = cos(tables%psi(k) * tables%u(i))
                tables%sin_odd(i, k) = sin(tables%psi(k) * tables%u(i))
            end do
        end do

        ! Beak-scan grid: same construction, clamped away from the poles
        allocate(tables%u_beak(FOS_BEAK_SCAN_POINTS), stat = alloc_stat)
        if (alloc_stat /= 0) then
            call reject_tables_s(tables, status)
            return
        end if
        do i = 1_ik, FOS_BEAK_SCAN_POINTS
            u_raw = -1.0_rk + 2.0_rk * real(i - 1_ik, rk) &
                    / real(FOS_BEAK_SCAN_POINTS - 1_ik, rk)
            tables%u_beak(i) = max(-U_BEAK_CLAMP, min(U_BEAK_CLAMP, u_raw))
        end do

        allocate(tables%bk_cos_even(FOS_BEAK_SCAN_POINTS, k_max), &
                tables%bk_sin_odd(FOS_BEAK_SCAN_POINTS, k_max), stat = alloc_stat)
        if (alloc_stat /= 0) then
            call reject_tables_s(tables, status)
            return
        end if
        do k = 1_ik, k_max
            do i = 1_ik, FOS_BEAK_SCAN_POINTS
                tables%bk_cos_even(i, k) = cos(tables%omega(k) * tables%u_beak(i))
                tables%bk_sin_odd(i, k) = sin(tables%psi(k) * tables%u_beak(i))
            end do
        end do

        ! Theta nodes (absent for the theta-less one-shot forms)
        if (n_theta > 0_ik) then
            allocate(tables%thetas(n_theta), stat = alloc_stat)
            if (alloc_stat /= 0) then
                call reject_tables_s(tables, status)
                return
            end if
            do i = 1_ik, n_theta
                tables%thetas(i) = thetas(i)
            end do
        end if

        tables%initialized = .true.
        status = SHAPE_VALID

    end subroutine build_tables_s

    !> Releases a half-built `tables_t` and reports the grid as unsatisfiable.
    !!
    !! The single exit `build_tables_s` uses when an allocation fails: whatever
    !! was already taken is given back, so the caller sees the same freed state a
    !! never-built tables object has.
    !!
    !! @param[in,out] tables  Partially built tables, freed here
    !! @param[out]    status  Always SHAPE_ERROR_INVALID_GRID
    pure subroutine reject_tables_s(tables, status)

        type(tables_t), intent(inout) :: tables
        integer(kind = ik), intent(out) :: status

        call tables_free_s(tables)
        status = SHAPE_ERROR_INVALID_GRID

    end subroutine reject_tables_s

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

        ! Trimmed, so a short vector and its zero-padded form run the same loop
        ! on the same data.
        a2 = a2_f(params(1:active_length_f(params)))
        status = SHAPE_VALID

    end subroutine compute_a2_s

    !> Computes the intrinsic z-shift that places the COM at the origin.
    !!
    !! @param[in]  params   FoS parameters; params(1) = c is mandatory
    !! @param[out] z_shift  Intrinsic shift (0 on rejection)
    !! @param[out] status   SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS,
    !!                      SHAPE_ERROR_WRONG_PARAM_COUNT (empty vector), or
    !!                      FOS_ERROR_INVALID_C (c <= C_MIN)
    pure subroutine compute_z_shift_s(params, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        real(kind = rk) :: c, sum_term, a_odd, sign_factor
        integer(kind = ik) :: n, n_active

        z_shift = 0.0_rk

        if (size(params, kind = ik) > FOS_MAX_K) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        ! An empty vector has no c: the contract's code for it is 4, as in every
        ! routine that needs params(1).
        if (size(params, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if

        c = params(1)
        if (c <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        ! Trimmed like `compute_a2_s`; c > C_MIN makes n_active >= 1.
        n_active = active_length_f(params)

        sum_term = 0.0_rk
        do n = 1_ik, FOS_MAX_K
            a_odd = coefficient_f(params(1:n_active), 2_ik * n + 1_ik, 0.0_rk)
            if (abs(a_odd) < FOS_COEFF_NEGLIGIBLE) cycle
            sign_factor = merge(-1.0_rk, 1.0_rk, mod(n, 2_ik) == 0_ik)
            sum_term = sum_term + sign_factor * a_odd / real(n, rk)
        end do

        z_shift = (3.0_rk / (2.0_rk * PI_C)) * c * sum_term
        status = SHAPE_VALID

    end subroutine compute_z_shift_s

    !===========================================================================
    ! PARAMETER-VECTOR HELPERS
    !===========================================================================

    !> Index of the last nonzero parameter; 0 for an empty or all-zero vector.
    !!
    !! Both tiers trim a vector to this length before any kernel runs, so a
    !! short vector and its zero-padded form reach the kernels as the same data
    !! with the same trip counts. The test is `abs(x) > 0`, so a trailing -0.0
    !! counts as zero. Interior zeros are kept.
    !!
    !! @param[in] params  Parameter vector, any length
    !! @return            Largest k with params(k) /= 0, or 0
    pure function active_length_f(params) result(n_active)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik) :: n_active

        integer(kind = ik) :: k

        n_active = 0_ik
        do k = size(params, kind = ik), 1_ik, -1_ik
            if (abs(params(k)) > 0.0_rk) then
                n_active = k
                exit
            end if
        end do

    end function active_length_f

    !> Fourier orders a vector of `n_params` parameters needs.
    !!
    !! (n_params + 2)/2 + 1, capped at FOS_MAX_K: the formula of the live
    !! evaluator `eval_f_s`. It sizes a cache's tables from `max_params` and
    !! bounds the kernel loops from the vector's own length.
    !!
    !! @param[in] n_params  Parameter count, >= 0
    !! @return              Number of Fourier orders
    pure function table_order_f(n_params) result(k_max)

        integer(kind = ik), intent(in) :: n_params
        integer(kind = ik) :: k_max

        k_max = min((n_params + 2_ik) / 2_ik + 1_ik, FOS_MAX_K)

    end function table_order_f

    !> Total z-shift for an origin moved by `s` from the COM frame.
    !!
    !! One addition, kept in this module so the tier code performs no arithmetic
    !! on an output value.
    !!
    !! @param[in] z_shift_intrinsic  COM shift
    !! @param[in] s                  Additional origin shift
    !! @return                       z_shift_intrinsic + s
    pure function shifted_origin_f(z_shift_intrinsic, s) result(z_shift)

        real(kind = rk), intent(in) :: z_shift_intrinsic
        real(kind = rk), intent(in) :: s
        real(kind = rk) :: z_shift

        z_shift = z_shift_intrinsic + s

    end function shifted_origin_f

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
    !! covers the vector, k_max >= (size(params) + 2)/2 + 1. The sum runs over
    !! the vector's orders, `table_order_f(size(params))`, never further.
    !!
    !! @param[in]  tables   Initialized trig tables
    !! @param[in]  params   FoS parameters (params(1) = c is not read)
    !! @param[out] f_grid   f at tables%u
    !! @param[out] fp_grid  df/du at tables%u
    pure subroutine compute_f_grid_s(tables, params, f_grid, fp_grid)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: f_grid(:)
        real(kind = rk), intent(out) :: fp_grid(:)

        real(kind = rk) :: a_even(tables%k_max), a_odd(tables%k_max)
        logical :: active(tables%k_max)
        real(kind = rk) :: sum_f, sum_fp, u
        integer(kind = ik) :: i, k, k_hi

        ! The loop runs over the orders the VECTOR needs, not the orders the
        ! table holds: tables of different order then execute the same trip
        ! count on the same data, which is what makes a larger cache
        ! bit-identical to a smaller one.
        k_hi = min(table_order_f(size(params, kind = ik)), tables%k_max)

        call pair_coefficients_s(params, k_hi, a_even, a_odd, active)

        do i = 1_ik, tables%n_points
            sum_f = 0.0_rk
            sum_fp = 0.0_rk
            do k = 1_ik, k_hi
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
    !! @param[in]  tables        Initialized trig tables
    !! @param[in]  params        FoS parameters (params(1) = c is not read)
    !! @param[out] f_min         Smallest f over the scan grid
    !! @param[out] beak_ok       .true. iff f_min > F_MIN_THRESHOLD
    !! @param[out] u_at_min      Optional: scan-grid u of the minimum
    !! @param[out] interior_min  Optional: .true. iff the minimum is a strict
    !!                           interior scan point (not at the clamp) — a
    !!                           minimum away from the clamp IS a local minimum
    !!                           of the scan, so no neighbor test is needed
    pure subroutine beak_scan_f_min_s(tables, params, f_min, beak_ok, u_at_min, &
            interior_min)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: f_min
        logical, intent(out) :: beak_ok
        real(kind = rk), intent(out), optional :: u_at_min
        logical, intent(out), optional :: interior_min

        real(kind = rk) :: a_even(tables%k_max), a_odd(tables%k_max)
        logical :: active(tables%k_max)
        real(kind = rk) :: sum_f, f_val, u
        integer(kind = ik) :: i, k, i_min, k_hi

        ! Same bound as `compute_f_grid_s`: the vector's orders, not the table's.
        k_hi = min(table_order_f(size(params, kind = ik)), tables%k_max)

        call pair_coefficients_s(params, k_hi, a_even, a_odd, active)

        f_min = huge(1.0_rk)
        i_min = 1_ik
        do i = 1_ik, FOS_BEAK_SCAN_POINTS
            sum_f = 0.0_rk
            do k = 1_ik, k_hi
                if (.not. active(k)) cycle
                sum_f = sum_f + a_even(k) * tables%bk_cos_even(i, k) &
                        + a_odd(k) * tables%bk_sin_odd(i, k)
            end do
            u = tables%u_beak(i)
            f_val = 1.0_rk - u**2 - sum_f
            if (f_val < f_min) then
                f_min = f_val
                i_min = i
            end if
        end do

        beak_ok = f_min > F_MIN_THRESHOLD
        if (present(u_at_min)) u_at_min = tables%u_beak(i_min)
        if (present(interior_min)) interior_min = &
                i_min > 1_ik .and. i_min < FOS_BEAK_SCAN_POINTS

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
    pure subroutine scale_rho_grid_s(tables, c, z_shift_intrinsic, f_grid, fp_grid, &
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

            if (abs(u) >= 1.0_rk - FOS_U_TIP_TOL) then
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
    ! ORIGIN RESOLVE (STAR-CONVEXITY)
    !===========================================================================

    !> Resolves the R(theta) origin of a COM-shifted rho(z) grid.
    !!
    !! g(s) = max_i[(z_i + s) drho_dz_i - rho_i] over the interior nodes is the
    !! star-convexity test value at additional shift s: the shape is
    !! single-valued in R(theta) about that origin iff g(s) <= 0, and this
    !! library accepts it iff g(s) <= -STAR_CONVEXITY_MARGIN. g is a pointwise
    !! max of affine functions of s, hence convex, so the golden-section search
    !! over [-2c, 2c] returns the global minimum.
    !!
    !! The verdict and the origin are two separate decisions, exactly as in 1.x:
    !!   - accept iff g(s*) <= -STAR_CONVEXITY_MARGIN — the deepest achievable
    !!     value is the true representability test;
    !!   - then keep the COM origin (additional shift 0) iff it is already
    !!     well-conditioned, g(0) <= -ORIGIN_CONDITION_MARGIN, else move to s*.
    !!     Evaluating a marginal-at-COM shape from its steep COM origin is what
    !!     corrupts the volume integral.
    !! There is no lazy path: g(0) is computed AND the minimizer always runs.
    !!
    !! Infallible: the caller gates on `star_ok`. `z` arrives already COM-shifted
    !! (scale_rho_grid_s baked z_shift_intrinsic in), so s is measured from the
    !! COM frame and the total shift is z_shift_intrinsic + s.
    !!
    !! Preconditions (caller's, unchecked): the three arrays are the same length
    !! and hold a resolved rho(z) grid; c > C_MIN.
    !!
    !! @param[in]  z                  Axial coordinate, COM-shifted
    !! @param[in]  rho                Cylindrical radius
    !! @param[in]  drho_dz            Slope
    !! @param[in]  c                  Elongation, params(1)
    !! @param[in]  z_shift_intrinsic  COM shift already baked into z
    !! @param[out] g0                 g at the COM origin
    !! @param[out] s_opt              argmin_s g(s)
    !! @param[out] g_opt              g(s_opt)
    !! @param[out] star_ok            .true. iff g_opt <= -STAR_CONVEXITY_MARGIN
    !! @param[out] z_shift_total      Intrinsic shift plus the chosen origin
    !! @param[out] r_north            R(0) = c + z_shift_total
    !! @param[out] r_south            R(pi) = |-c + z_shift_total|
    pure subroutine resolve_origin_s(z, rho, drho_dz, c, z_shift_intrinsic, &
            g0, s_opt, g_opt, star_ok, z_shift_total, r_north, r_south)

        real(kind = rk), intent(in) :: z(:)
        real(kind = rk), intent(in) :: rho(:)
        real(kind = rk), intent(in) :: drho_dz(:)
        real(kind = rk), intent(in) :: c
        real(kind = rk), intent(in) :: z_shift_intrinsic
        real(kind = rk), intent(out) :: g0
        real(kind = rk), intent(out) :: s_opt
        real(kind = rk), intent(out) :: g_opt
        logical, intent(out) :: star_ok
        real(kind = rk), intent(out) :: z_shift_total
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south

        real(kind = rk) :: s_origin

        g0 = star_convexity_max_f(z, rho, drho_dz, 0.0_rk)
        call minimize_star_convexity_s(z, rho, drho_dz, c, s_opt, g_opt)

        star_ok = g_opt <= -STAR_CONVEXITY_MARGIN

        if (g0 <= -ORIGIN_CONDITION_MARGIN) then
            s_origin = 0.0_rk
        else
            s_origin = s_opt
        end if

        z_shift_total = z_shift_intrinsic + s_origin
        r_north = c + z_shift_total
        r_south = abs(-c + z_shift_total)

    end subroutine resolve_origin_s

    !> g(s) = max_i[(z_i + s) drho_dz_i - rho_i] over the INTERIOR nodes.
    !!
    !! The tips carry the rho = 0, drho/dz = 0 convention and would otherwise
    !! contribute a spurious 0; the 1.x scan skips them the same way.
    !!
    !! Precondition (caller's, unchecked): size(z) >= 3. A shorter grid has no
    !! interior node, and the empty max returns -huge — which every caller would
    !! read as a wildly star-convex shape. Unreachable in this library: every
    !! grid comes from `build_tables_s`, which floors n_points at
    !! FOS_N_POINTS_FLOOR = 100.
    pure function star_convexity_max_f(z, rho, drho_dz, s) result(g)

        real(kind = rk), intent(in) :: z(:)
        real(kind = rk), intent(in) :: rho(:)
        real(kind = rk), intent(in) :: drho_dz(:)
        real(kind = rk), intent(in) :: s
        real(kind = rk) :: g

        integer(kind = ik) :: i, n
        real(kind = rk) :: test_val

        n = size(z, kind = ik)

        g = -huge(1.0_rk)
        do i = 2_ik, n - 1_ik
            test_val = (z(i) + s) * drho_dz(i) - rho(i)
            if (test_val > g) g = test_val
        end do

    end function star_convexity_max_f

    !> Golden-section minimization of g(s) over the bracket [-2c, 2c], which
    !! contains any origin that can be star-convex. g is convex piecewise-linear
    !! in s, so the search returns the exact global minimum.
    pure subroutine minimize_star_convexity_s(z, rho, drho_dz, c, s_opt, g_opt)

        real(kind = rk), intent(in) :: z(:)
        real(kind = rk), intent(in) :: rho(:)
        real(kind = rk), intent(in) :: drho_dz(:)
        real(kind = rk), intent(in) :: c
        real(kind = rk), intent(out) :: s_opt
        real(kind = rk), intent(out) :: g_opt

        real(kind = rk) :: a, b, x1, x2, f1, f2
        integer(kind = ik) :: it

        a = -2.0_rk * c
        b = 2.0_rk * c
        x1 = b - GOLDEN * (b - a)
        x2 = a + GOLDEN * (b - a)
        f1 = star_convexity_max_f(z, rho, drho_dz, x1)
        f2 = star_convexity_max_f(z, rho, drho_dz, x2)

        do it = 1_ik, SHIFT_MAX_ITER
            if (b - a <= SHIFT_TOL) exit
            if (f1 < f2) then
                b = x2
                x2 = x1
                f2 = f1
                x1 = b - GOLDEN * (b - a)
                f1 = star_convexity_max_f(z, rho, drho_dz, x1)
            else
                a = x1
                x1 = x2
                f1 = f2
                x2 = a + GOLDEN * (b - a)
                f2 = star_convexity_max_f(z, rho, drho_dz, x2)
            end if
        end do

        s_opt = 0.5_rk * (a + b)
        g_opt = star_convexity_max_f(z, rho, drho_dz, s_opt)

    end subroutine minimize_star_convexity_s

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
    !! Two further safeguards (2026-08-11 geometry sweep, GL-4096 nodes):
    !!   - Bracket-collapse acceptance: within ~ulp(r_south)/r of a pole the
    !!     surface slope makes dF/dr * ulp(r) exceed NR_TOLERANCE, so NO
    !!     representable r passes the residual test. Once the sign-change
    !!     bracket is NR_BRACKET_ULPS wide the root is pinned to machine
    !!     precision and the node is accepted on that evidence instead.
    !!   - Bisection fallback: Newton can enter an attracting period-2 orbit
    !!     strictly inside the bracket (the sweep found one), where the
    !!     leaves-the-bracket safeguard is blind. Past NR_NEWTON_PHASE
    !!     iterations every step bisects, which converges unconditionally on
    !!     the sign-change bracket. Nodes that converge never see either
    !!     safeguard: the residual exit fires within the Newton phase, whose
    !!     update rule is unchanged.
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
        logical :: bracket_collapsed
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

        bracket_collapsed = .false.

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

            ! Bracket collapsed to machine precision: the root is pinned even
            ! though the residual floor dF/dr * ulp(r) sits above NR_TOLERANCE
            ! (pole boundary layer). Accept on the bracket evidence.
            if (r_hi - r_lo <= NR_BRACKET_ULPS &
                    * spacing(max(abs(r_lo), abs(r_hi)))) then
                bracket_collapsed = .true.
                exit
            end if

            if (iter <= NR_NEWTON_PHASE) then
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
            else
                ! Newton-resistant node (e.g. attracting 2-cycle inside the
                ! bracket): bisect unconditionally from here on.
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
        converged = abs(F_val) < NR_TOLERANCE * max(1.0_rk, r) &
                .or. bracket_collapsed

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
    ! R(THETA) SOLVE, NECK, AND PRIVATE HELPERS
    !===========================================================================

    !> The evaluator bundle for a resolved shape.
    !!
    !! `r_hi_bound` is the analytic Newton bracket: twice the distance from the
    !! shifted origin to the corner of the shape's bounding box
    !! (rho_max by the polar extents), which is outside the surface at every
    !! theta. It replaces the 1.x doubling loop, whose eight doublings silently
    !! returned the initial guess for extreme-oblate shapes.
    !!
    !! Preconditions (caller's, unchecked): 1 <= size(params) <= FOS_MAX_K,
    !! params(1) > C_MIN, and `rho_max` / `z_shift_total` from an up-to-date
    !! resolve. The tier module guarantees the length bound: it accepts at most
    !! FOS_MAX_PARAMS parameters and passes the trimmed vector.
    !!
    !! @param[in] params         Trimmed parameter vector
    !! @param[in] z_shift_total  Intrinsic COM shift plus the chosen origin
    !! @param[in] rho_max        Largest rho on the resolved grid
    !! @return                   Bundle for `newton_radius_s`
    pure function fos_bundle_f(params, z_shift_total, rho_max) result(bundle)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: z_shift_total
        real(kind = rk), intent(in) :: rho_max
        type(fos_bundle_t) :: bundle

        real(kind = rk) :: z_max, z_min
        integer(kind = ik) :: i, n

        n = size(params, kind = ik)
        bundle%n_params = n
        do i = 1_ik, n
            bundle%params(i) = params(i)
        end do
        bundle%z_shift = z_shift_total

        z_max = params(1) + z_shift_total
        z_min = -params(1) + z_shift_total
        bundle%r_hi_bound = 2.0_rk * sqrt(rho_max**2 + max(z_max, abs(z_min))**2)

    end function fos_bundle_f

    !> The ONE R(theta) solve loop in the library: cos(theta) at every node, the
    !! Newton root, and the convergence verdict folded into a single status.
    !!
    !! Every R(theta) output of the tier module — the fixed-grid forms and the
    !! at-thetas form alike — calls this routine, and that is a correctness
    !! requirement, not tidiness. The contract is that identical thetas in give
    !! bit-identical radii and derivatives out, in every build configuration.
    !! Two mechanisms enforce it, and both are needed:
    !!
    !!   - ONE cos evaluation site. The fixed-grid forms pass `tables%thetas`,
    !!     not a precomputed cosine table, because under `-ffast-math` the
    !!     compiler may lower a vectorizable cos loop (the one in
    !!     `build_tables_s`) to libmvec and a scalar one to the libm call, whose
    !!     results differ by an ulp — and one ulp on cos(theta) moves the Newton
    !!     iterate.
    !!   - ONE machine-code copy. Every caller is in the tier module, a
    !!     different file compiled without LTO, so this body cannot be inlined
    !!     into any of them. (Through 2.1.0 the callers shared this file and a
    !!     `noinline` attribute did the job on GCC >= 12 only; the manylinux
    !!     wheel's GCC 10 ignored it.)
    !!
    !! Preconditions (caller's, unchecked): `bundle` from a resolved,
    !! representable shape; `thetas` in [0, pi]; `radii` and `dr_dtheta` at
    !! least size(thetas) long.
    !!
    !! @param[in]  bundle     Resolved shape (params, z_shift, bracket bound)
    !! @param[in]  thetas     Polar angles in [0, pi]
    !! @param[out] radii      R(theta) at those angles
    !! @param[out] dr_dtheta  dR/dtheta at those angles
    !! @param[out] status     SHAPE_VALID, or FOS_ERROR_CONVERGENCE if ANY node
    !!                        missed NR_TOLERANCE
    pure subroutine solve_thetas_s(bundle, thetas, radii, dr_dtheta, status)

        type(fos_bundle_t), intent(in) :: bundle
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i
        logical :: converged

        status = SHAPE_VALID

        do i = 1_ik, size(thetas, kind = ik)
            call newton_radius_s(bundle, cos(thetas(i)), radii(i), dr_dtheta(i), &
                    converged)
            if (.not. converged) status = FOS_ERROR_CONVERGENCE
        end do

    end subroutine solve_thetas_s

    !> Grid index of the neck: the smallest rho between the two largest interior
    !! rho maxima. Reproduces the 1.x scan node for node.
    !!
    !! Fewer than two maxima means no neck — a shape property, not a failure.
    pure subroutine find_neck_index_s(rho, neck_idx, found)

        real(kind = rk), intent(in) :: rho(:)
        integer(kind = ik), intent(out) :: neck_idx
        logical, intent(out) :: found

        integer(kind = ik) :: i, j, n, n_maxima, max1_idx, max2_idx
        integer(kind = ik) :: left_idx, right_idx
        real(kind = rk) :: max1_rho, max2_rho, min_rho

        neck_idx = 0_ik
        found = .false.

        n = size(rho, kind = ik)

        n_maxima = 0_ik
        max1_idx = 0_ik
        max2_idx = 0_ik
        max1_rho = -1.0_rk
        max2_rho = -1.0_rk

        do i = 2_ik, n - 1_ik
            if (rho(i) > rho(i - 1_ik) .and. rho(i) > rho(i + 1_ik)) then
                n_maxima = n_maxima + 1_ik
                if (rho(i) > max1_rho) then
                    max2_rho = max1_rho
                    max2_idx = max1_idx
                    max1_rho = rho(i)
                    max1_idx = i
                else if (rho(i) > max2_rho) then
                    max2_rho = rho(i)
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
            if (rho(j) < min_rho) then
                min_rho = rho(j)
                neck_idx = j
            end if
        end do

        found = .true.

    end subroutine find_neck_index_s

    !> Neck of a resolved rho(z) grid: the interior rho minimum between the two
    !! largest rho maxima, refined to machine precision.
    !!
    !! A coarse scan over the grid brackets the neck, then Newton iteration on
    !! f'(u) = 0, bracketed to one grid spacing around the scan result so it
    !! cannot escape to a different extremum, refines it. `found` is .false.,
    !! with both outputs zero, for a profile with fewer than two rho maxima:
    !! having no neck is an answer, not an error.
    !!
    !! Preconditions (caller's, unchecked): `tables` initialized, `rho` the
    !! tables%n_points-long grid of this vector, params(1) > C_MIN.
    !!
    !! @param[in]  tables             Initialized trig tables (for the u nodes)
    !! @param[in]  params             FoS parameters
    !! @param[in]  rho                Cylindrical radius on the u grid
    !! @param[in]  z_shift_intrinsic  COM shift of the profile
    !! @param[out] z_neck             Neck z-position in the COM frame
    !! @param[out] rho_neck           Neck radius
    !! @param[out] found              .true. iff the profile has a neck
    pure subroutine refine_neck_s(tables, params, rho, z_shift_intrinsic, &
            z_neck, rho_neck, found)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: rho(:)
        real(kind = rk), intent(in) :: z_shift_intrinsic
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found

        integer(kind = ik) :: neck_idx, iter
        real(kind = rk) :: c, u, u_lo, u_hi, du, step
        real(kind = rk) :: f_val, fp_val, fpp_val

        z_neck = 0.0_rk
        rho_neck = 0.0_rk

        call find_neck_index_s(rho, neck_idx, found)
        if (.not. found) return

        c = params(1)
        du = 2.0_rk / real(tables%n_points - 1_ik, rk)
        u = tables%u(neck_idx)
        u_lo = max(-1.0_rk, u - du)
        u_hi = min(1.0_rk, u + du)

        do iter = 1_ik, NECK_NEWTON_MAX_ITER
            call eval_f_s(params, u, f_val, fp_val, fpp_val)
            if (fpp_val <= 0.0_rk) exit
            step = fp_val / fpp_val
            u = min(u_hi, max(u_lo, u - step))
            if (abs(step) < NECK_NEWTON_TOL) exit
        end do

        call eval_f_s(params, u, f_val, fp_val)
        rho_neck = sqrt(max(f_val, 0.0_rk) / c)
        z_neck = c * u + z_shift_intrinsic

    end subroutine refine_neck_s

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
            if (abs(a_2n) < FOS_COEFF_NEGLIGIBLE) cycle
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
    !! coefficients below FOS_COEFF_NEGLIGIBLE is dropped from the sum, not added as
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
            active(k) = abs(a_even(k)) >= FOS_COEFF_NEGLIGIBLE &
                    .or. abs(a_odd(k)) >= FOS_COEFF_NEGLIGIBLE
        end do

    end subroutine pair_coefficients_s

    !> Live (untabled) f(u), f'(u) and — on request — f''(u). The Newton paths
    !! evaluate f at arbitrary u, off every grid node; only the neck refinement
    !! needs the curvature, so it is optional rather than always summed.
    pure subroutine eval_f_s(params, u, f, fp, fpp)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: u
        real(kind = rk), intent(out) :: f
        real(kind = rk), intent(out) :: fp
        real(kind = rk), intent(out), optional :: fpp

        integer(kind = ik) :: k, k_max, n_params
        real(kind = rk) :: a2, omega_k, psi_k, a_even, a_odd
        real(kind = rk) :: cos_even, sin_even, cos_odd, sin_odd
        real(kind = rk) :: sum_f, sum_fp, sum_fpp
        logical :: need_fpp

        need_fpp = present(fpp)

        a2 = a2_f(params)
        n_params = size(params, kind = ik)
        k_max = min((n_params + 2_ik) / 2_ik + 1_ik, FOS_MAX_K)

        sum_f = 0.0_rk
        sum_fp = 0.0_rk
        sum_fpp = 0.0_rk

        do k = 1_ik, k_max
            a_even = coefficient_f(params, 2_ik * k, a2)
            a_odd = coefficient_f(params, 2_ik * k + 1_ik, a2)
            if (abs(a_even) < FOS_COEFF_NEGLIGIBLE .and. abs(a_odd) < FOS_COEFF_NEGLIGIBLE) cycle

            omega_k = real(2_ik * k - 1_ik, rk) * PI_C / 2.0_rk
            psi_k = real(k, rk) * PI_C

            cos_even = cos(omega_k * u)
            sin_even = sin(omega_k * u)
            cos_odd = cos(psi_k * u)
            sin_odd = sin(psi_k * u)

            sum_f = sum_f + a_even * cos_even + a_odd * sin_odd
            sum_fp = sum_fp - a_even * omega_k * sin_even + a_odd * psi_k * cos_odd
            if (need_fpp) then
                sum_fpp = sum_fpp - a_even * omega_k**2 * cos_even &
                        - a_odd * psi_k**2 * sin_odd
            end if
        end do

        f = 1.0_rk - u**2 - sum_f
        fp = -2.0_rk * u - sum_fp
        if (need_fpp) fpp = -2.0_rk - sum_fpp

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
        if (abs(u) >= 1.0_rk - FOS_U_TIP_TOL) return

        call eval_f_s(bundle%params(1:bundle%n_params), u, f_val, fp_val)

        if (f_val > 0.0_rk) then
            sqrt_cf = sqrt(c * f_val)
            rho = sqrt(f_val * c_inv)
            drho_dz = fp_val / (2.0_rk * c * sqrt_cf)
        end if

    end subroutine bundle_rho_s

end module fos_parameterization_workers_mod
