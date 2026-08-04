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
!!
!! `cache_t` is the per-shape working level: one owner, one thread. It embeds a
!! `shape_engine_t` recompute tracker and every per-shape buffer, and it is
!! mutated by every compute call. See the ownership and lifetime notes on the
!! type itself.
!!
!! ## Cached intermediates (normative dependency map)
!!
!! Parameter slots: p1 = c, p2 = a3, p3 = a4, p4 = a5, p5 = a6, p6 = a7,
!! p7 = a8, p8 = a9.
!!
!! | # | intermediate                          | depends on      |
!! |---|---------------------------------------|-----------------|
!! | 1 | a2 (volume constraint)                | p3, p5, p7      |
!! | 2 | z_shift_intrinsic                     | p1, p2, p4, p6, p8 |
!! | 3 | f_grid: f, f' at the u nodes          | p2-p8 (NOT c)   |
!! | 4 | beak verdict: f_min over the scan grid| p2-p8 (NOT c)   |
!! | 5 | rho_grid: z, rho, drho/dz, rho_max    | p1-p8           |
!! | 6 | resolve: g(0), optimum, z_shift_total | p1-p8           |
!! | 7 | radius grid at the init thetas        | p1-p8           |
!! | 8 | radius and derivative at those thetas | p1-p8           |
!!
!! Producer order is 1 -> 8. Every cached compute presents its parameters to
!! the engine once (`ensure_list_s`), then recomputes only the intermediates
!! its own output needs, in that order.
module fos_parameterization_workers_mod

    use precision_utilities_mod, only: ik, ikl, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_TOO_MANY_PARAMS, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED, SHAPE_CACHE_MAX_PARAMS, &
            shape_engine_t, shape_engine_init_s, shape_engine_begin_s, &
            shape_engine_needs_f, shape_engine_note_computed_s, &
            shape_engine_invalidate_all_s, shape_engine_recompute_count_f

    implicit none

    private

    ! Types
    public :: tables_t
    public :: cache_t
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
    public :: resolve_origin_s
    public :: newton_radius_s

    ! Cache lifecycle, cached outputs, and engine introspection
    public :: fos_masks_f
    public :: cache_init_s
    public :: cache_init_shared_s
    public :: cache_free_s
    public :: cache_shape_s
    public :: cache_rho_z_grid_s
    public :: cache_radius_grid_s
    public :: cache_radius_and_derivative_s
    public :: cache_radius_and_derivative_at_thetas_s
    public :: cache_neck_s
    public :: cache_star_convexity_optimum_s
    public :: cache_recompute_count_f

    ! Standalone (tier-1) forms: one call in, one answer out, nothing retained
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s
    public :: compute_shape_standalone_s
    public :: compute_rho_z_grid_standalone_s
    public :: compute_neck_standalone_s
    public :: compute_star_convexity_optimum_standalone_s

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
    ! This module is the single owner of every kernel constant. The 1.x
    ! duplicates in fos_parameterization_mod were deleted with the 1.x surface;
    ! the two constants the raw evaluators there still need
    ! (FOS_COEFF_NEGLIGIBLE, FOS_U_TIP_TOL) are public and imported from here.

    !> Largest Fourier coefficient index the evaluator reads.
    integer(kind = ik), parameter, public :: FOS_MAX_K = 50_ik

    !> N_max: the longest parameter vector the standalone (tier-1) forms accept.
    !! Equal to FOS_MAX_K = 50 by construction — the evaluator reads no
    !! coefficient past a_50 and `fos_bundle_t` carries exactly that many slots,
    !! so a longer vector could not be evaluated even if it were accepted. The
    !! cached tier stops much earlier, at SHAPE_CACHE_MAX_PARAMS, because its
    !! tables are sized once at init; tier-1 sizes its tables per call and so
    !! goes all the way to N_max.
    integer(kind = ik), parameter, public :: FOS_MAX_PARAMS = FOS_MAX_K

    !> Coefficient pairs below this magnitude contribute nothing and are
    !! skipped — skipping is load-bearing for bitwise parity between the tabled
    !! kernels and the live evaluator. Public: the raw evaluators in
    !! fos_parameterization_mod share it.
    real(kind = rk), parameter, public :: FOS_COEFF_NEGLIGIBLE = 1.0e-30_rk

    !> Smallest elongation the evaluator treats as a shape.
    real(kind = rk), parameter, public :: C_MIN = 1.0e-10_rk

    !> Beak-singularity threshold: shapes whose f_min falls below it are rejected.
    real(kind = rk), parameter, public :: F_MIN_THRESHOLD = 1.0e-3_rk

    !> Tip detection tolerance: f(±1) = 0 analytically, so u within roundoff of
    !! a tip is treated AS the tip (rho = 0, drho/dz = 0). Public: the raw
    !! evaluators in fos_parameterization_mod share it.
    real(kind = rk), parameter, public :: FOS_U_TIP_TOL = 4.0_rk * epsilon(1.0_rk)

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

    !> Output array whose size differs from the cache's init-time resolution.
    !! Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_BUFFER_MISMATCH = 105_ik

    !> Shape not star-convex from any origin. Owned here, re-exported by the
    !! main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_NOT_STAR_CONVEX = 101_ik

    !> f_min below the beak threshold. Owned here, re-exported by the main module.
    integer(kind = ik), parameter, public :: FOS_ERROR_BEAK_SINGULARITY = 103_ik

    !---------------------------------------------------------------------------
    ! Cached intermediates
    !---------------------------------------------------------------------------
    !> Number of intermediates the engine tracks (see the dependency map in the
    !! module header). Must equal the extent of `fos_masks_f`'s result.
    integer(kind = ik), parameter, public :: FOS_N_INTERMEDIATES = 8_ik

    !> Intermediate indices, in producer order.
    integer(kind = ik), parameter, public :: I_A2 = 1_ik
    integer(kind = ik), parameter, public :: I_Z_SHIFT = 2_ik
    integer(kind = ik), parameter, public :: I_F_GRID = 3_ik
    integer(kind = ik), parameter, public :: I_BEAK = 4_ik
    integer(kind = ik), parameter, public :: I_RHO_GRID = 5_ik
    integer(kind = ik), parameter, public :: I_RESOLVE = 6_ik
    integer(kind = ik), parameter, public :: I_RADIUS_GRID = 7_ik
    integer(kind = ik), parameter, public :: I_RADIUS_DERIV = 8_ik

    !> Full-width dependency masks, one per intermediate: bit j-1 set means the
    !! intermediate depends on parameter j. `fos_masks_f` truncates them to the
    !! declared parameter count.
    integer(kind = ik), parameter :: FOS_MASKS_FULL(FOS_N_INTERMEDIATES) = &
            [84_ik, 171_ik, 254_ik, 254_ik, 255_ik, 255_ik, 255_ik, 255_ik]

    !> Needs-list of the cylindrical output: a2, z_shift, f-grid, rho-grid.
    !! The beak scan (#4) is deliberately absent — it gates the R(theta)
    !! conversion, not the cylindrical profile — so it is neither computed nor
    !! stamped by this path.
    integer(kind = ik), parameter :: NEEDS_RHO_Z_GRID(4) = &
            [I_A2, I_Z_SHIFT, I_F_GRID, I_RHO_GRID]

    !> Needs-list of the resolved shape: everything up to and including the
    !! origin search. The beak verdict (#4) IS part of it — the resolve exists
    !! only to serve the R(theta) conversion, which a beak singularity breaks.
    integer(kind = ik), parameter :: NEEDS_SHAPE(6) = &
            [I_A2, I_Z_SHIFT, I_F_GRID, I_BEAK, I_RHO_GRID, I_RESOLVE]

    !> Needs-lists of the two fixed-grid R(theta) outputs: the resolve, plus the
    !! one radius intermediate each consumes. The at-thetas form shares the
    !! resolve prefix but stamps neither #7 nor #8 — its nodes are the caller's,
    !! not the tables', so nothing it computes can be reused.
    integer(kind = ik), parameter :: NEEDS_RADIUS_GRID(7) = &
            [I_A2, I_Z_SHIFT, I_F_GRID, I_BEAK, I_RHO_GRID, I_RESOLVE, I_RADIUS_GRID]
    integer(kind = ik), parameter :: NEEDS_RADIUS_DERIV(7) = &
            [I_A2, I_Z_SHIFT, I_F_GRID, I_BEAK, I_RHO_GRID, I_RESOLVE, I_RADIUS_DERIV]

    !> Parameter-independent trigonometric tables for FoS shape evaluation.
    !!
    !! Immutable after `tables_init_s`, therefore shareable across threads. A
    !! tables object bound to a cache with `cache_init_shared_s` MUST be
    !! declared `target` and MUST outlive every cache bound to it — the cache
    !! stores a pointer, not a copy.
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

    !> Per-shape working level: the embedded recompute tracker, the resolved
    !! scalars, and every per-shape buffer.
    !!
    !! ## Ownership
    !!
    !! `tables` points at the trig tables. Two modes:
    !!   - `cache_init_s` (private): the cache heap-allocates its own `tables_t`
    !!     through `tables` and sets `owns_tables = .true.`; `cache_free_s`
    !!     frees it.
    !!   - `cache_init_shared_s`: `tables` points at a caller-owned `tables_t`
    !!     and `owns_tables = .false.`. **The caller MUST declare that
    !!     `tables_t` with the `target` attribute**, keep it alive, and not free
    !!     it for the whole lifetime of the cache.
    !!
    !! ## `cache_free_s` is mandatory — there is no finalizer
    !!
    !! Both init routines take `intent(out) :: cache`, which default-initializes
    !! on entry and so nulls the pointer before the previous target can be
    !! released. A private-mode cache therefore leaks its heap `tables_t` if you
    !! re-initialize a live cache or let it go out of scope unfreed.
    !!
    !! ## Never copy-assign a `cache_t`
    !!
    !! Intrinsic assignment copies `tables` and `owns_tables` shallowly,
    !! producing a dangling pointer or a double free. One cache has one owner
    !! and one thread; share the `tables_t` instead.
    type :: cache_t
        type(shape_engine_t) :: engine                    !! Recompute tracker

        type(tables_t), pointer :: tables => null()       !! Trig tables (see Ownership)
        logical :: owns_tables = .false.                  !! .true. => cache_free_s frees them

        integer(kind = ik) :: n_params = 0_ik             !! Fixed at init

        real(kind = rk) :: a2 = 0.0_rk                    !! #1
        real(kind = rk) :: z_shift_intrinsic = 0.0_rk     !! #2

        real(kind = rk), allocatable :: f_grid(:)         !! #3, n_points
        real(kind = rk), allocatable :: fp_grid(:)        !! #3, n_points

        real(kind = rk) :: f_min = 0.0_rk                 !! #4
        logical :: beak_ok = .false.                      !! #4 verdict

        real(kind = rk), allocatable :: z(:)              !! #5, n_points
        real(kind = rk), allocatable :: rho(:)            !! #5, n_points
        real(kind = rk), allocatable :: drho_dz(:)        !! #5, n_points
        real(kind = rk) :: rho_max = 0.0_rk               !! #5, Newton bracket input
        logical :: rho_positive = .false.                 !! #5 verdict

        real(kind = rk) :: g0 = 0.0_rk                    !! #6, g at the COM origin
        real(kind = rk) :: s_opt = 0.0_rk                 !! #6, best origin
        real(kind = rk) :: g_opt = 0.0_rk                 !! #6, g at s_opt
        real(kind = rk) :: z_shift_total = 0.0_rk         !! #6
        real(kind = rk) :: r_north = 0.0_rk               !! #6
        real(kind = rk) :: r_south = 0.0_rk               !! #6
        logical :: star_ok = .false.                      !! #6 verdict

        real(kind = rk), allocatable :: radii(:)          !! #7, n_theta
        real(kind = rk), allocatable :: dr_scratch(:)     !! #7 sink, n_theta
        real(kind = rk), allocatable :: radii_d(:)        !! #8, n_theta
        real(kind = rk), allocatable :: dr_dtheta(:)      !! #8, n_theta
    end type cache_t

    !> Per-call working state of the standalone (tier-1) tier: local tables plus
    !! the same intermediates `cache_t` stores, minus the engine — nothing here
    !! survives the call that built it, so there is nothing to track.
    !!
    !! Never public, never returned: the calling form copies out what it needs
    !! and the state dies with it. Every allocatable component (the tables'
    !! included) is released automatically when that form returns, which is why
    !! the standalone tier needs no free routine.
    type :: sa_state_t
        type(tables_t) :: tables

        real(kind = rk) :: a2 = 0.0_rk                    !! #1
        real(kind = rk) :: z_shift_intrinsic = 0.0_rk     !! #2

        real(kind = rk), allocatable :: f_grid(:)         !! #3, n_points
        real(kind = rk), allocatable :: fp_grid(:)        !! #3, n_points

        real(kind = rk) :: f_min = 0.0_rk                 !! #4
        logical :: beak_ok = .false.                      !! #4 verdict

        real(kind = rk), allocatable :: z(:)              !! #5, n_points
        real(kind = rk), allocatable :: rho(:)            !! #5, n_points
        real(kind = rk), allocatable :: drho_dz(:)        !! #5, n_points
        real(kind = rk) :: rho_max = 0.0_rk               !! #5
        logical :: rho_positive = .false.                 !! #5 verdict

        real(kind = rk) :: g0 = 0.0_rk                    !! #6
        real(kind = rk) :: s_opt = 0.0_rk                 !! #6
        real(kind = rk) :: g_opt = 0.0_rk                 !! #6
        real(kind = rk) :: z_shift_total = 0.0_rk         !! #6
        real(kind = rk) :: r_north = 0.0_rk               !! #6
        real(kind = rk) :: r_south = 0.0_rk               !! #6
        logical :: star_ok = .false.                      !! #6 verdict
    end type sa_state_t

    !> The empty theta set of the four theta-less standalone forms.
    real(kind = rk), parameter :: NO_THETAS(0) = [real(kind = rk) ::]

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
    !! then `thetas` stays unallocated with `n_theta = 0`.
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
            allocate(tables%thetas(n_theta))
            do i = 1_ik, n_theta
                tables%thetas(i) = thetas(i)
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
            if (abs(a_odd) < FOS_COEFF_NEGLIGIBLE) cycle
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
    !! read as a wildly star-convex shape. Unreachable in this library: both
    !! tiers route their grid through `build_tables_s`, which floors n_points at
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
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Dependency masks for the embedded engine, truncated to `n_params` bits.
    !!
    !! An out-of-range `n_params` yields all-zero masks: `maskr` is undefined
    !! outside [0, bit_size], and the engine rejects the count before the mask
    !! values can matter.
    !!
    !! @param[in] n_params  Requested parameter count (may be out of range)
    !! @return              One mask per tracked intermediate
    pure function fos_masks_f(n_params) result(masks)

        integer(kind = ik), intent(in) :: n_params
        integer(kind = ik) :: masks(FOS_N_INTERMEDIATES)

        integer(kind = ik) :: j

        masks = 0_ik
        if (n_params < 1_ik .or. n_params > SHAPE_CACHE_MAX_PARAMS) return

        do j = 1_ik, FOS_N_INTERMEDIATES
            masks(j) = iand(FOS_MASKS_FULL(j), maskr(n_params, ik))
        end do

    end function fos_masks_f

    !> Initializes a cache that owns its trig tables (private mode).
    !!
    !! The tables are heap-allocated through `cache%tables` and released by
    !! `cache_free_s`. Use `cache_init_shared_s` when many caches should share
    !! one `tables_t`.
    !!
    !! **Call `cache_free_s` before re-initializing a live cache and before it
    !! goes out of scope** — there is no finalizer, and `intent(out)` nulls the
    !! pointer on entry, so a re-init without a free leaks the previous tables.
    !!
    !! @param[out] cache     Ready on success; default (uninitialized) state otherwise
    !! @param[in]  n_params  1 <= n_params <= SHAPE_CACHE_MAX_PARAMS
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  thetas    Polar angles in [0, pi]; may be empty
    !! @param[out] status    SHAPE_VALID, SHAPE_ERROR_INVALID_INIT,
    !!                       SHAPE_ERROR_TOO_MANY_PARAMS or SHAPE_ERROR_INVALID_GRID
    subroutine cache_init_s(cache, n_params, n_points, thetas, status)

        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(in) :: n_params
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: masks(FOS_N_INTERMEDIATES)

        ! c is mandatory, so an empty parameter vector has no shape to describe.
        if (n_params < 1_ik) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        allocate(cache%tables)
        call build_tables_s(cache%tables, n_points, thetas, FOS_TABLES_K_MAX, status)
        if (status /= SHAPE_VALID) then
            deallocate(cache%tables)
            nullify(cache%tables)
            return
        end if
        cache%owns_tables = .true.

        ! The engine owns the upper end of the n_params contract: n_params > 8
        ! is rejected HERE, which is what keeps every cached path inside the
        ! FOS_TABLES_K_MAX = (8+2)/2+1 orders the tables above were built with.
        masks = fos_masks_f(n_params)
        call shape_engine_init_s(cache%engine, n_params, masks, status)
        if (status /= SHAPE_VALID) then
            call cache_free_s(cache)
            return
        end if

        call cache_alloc_buffers_s(cache)
        cache%n_params = n_params

    end subroutine cache_init_s

    !> Initializes a cache over caller-owned shared tables.
    !!
    !! The cache stores a pointer to `tables` and never frees it. **The caller
    !! MUST declare `tables` with the `target` attribute**; without it the
    !! association is undefined once this routine returns. `tables` must also
    !! outlive the cache and must not be freed while the cache is in use.
    !!
    !! @param[out] cache     Ready on success; default (uninitialized) state otherwise
    !! @param[in]  tables    Initialized shared tables; declared `target` by the caller
    !! @param[in]  n_params  1 <= n_params <= SHAPE_CACHE_MAX_PARAMS
    !! @param[out] status    SHAPE_VALID, SHAPE_ERROR_TABLES_NOT_INITIALIZED,
    !!                       SHAPE_ERROR_INVALID_INIT or SHAPE_ERROR_TOO_MANY_PARAMS
    subroutine cache_init_shared_s(cache, tables, n_params, status)

        type(cache_t), intent(out) :: cache
        type(tables_t), intent(in), target :: tables
        integer(kind = ik), intent(in) :: n_params
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: masks(FOS_N_INTERMEDIATES)

        if (.not. tables%initialized) then
            status = SHAPE_ERROR_TABLES_NOT_INITIALIZED
            return
        end if

        ! The engine owns the whole n_params contract here: < 1 is
        ! SHAPE_ERROR_INVALID_INIT, > SHAPE_CACHE_MAX_PARAMS is
        ! SHAPE_ERROR_TOO_MANY_PARAMS. Nothing is allocated before it agrees.
        masks = fos_masks_f(n_params)
        call shape_engine_init_s(cache%engine, n_params, masks, status)
        if (status /= SHAPE_VALID) return

        cache%tables => tables
        cache%owns_tables = .false.

        call cache_alloc_buffers_s(cache)
        cache%n_params = n_params

    end subroutine cache_init_shared_s

    !> Releases the cache. Infallible, and safe on an already-free instance.
    !!
    !! Frees the tables only in private mode; shared tables are left for their
    !! owner. The engine is reset to its pristine default so that a compute call
    !! on a freed cache reports SHAPE_ERROR_CACHE_NOT_INITIALIZED instead of
    !! dereferencing the dead tables pointer.
    !!
    !! @param[in,out] cache  Reset to the default (uninitialized) state
    subroutine cache_free_s(cache)

        type(cache_t), intent(inout) :: cache

        ! shape_engine_t has private components, so a structure constructor is
        ! not available: assign a default-initialized instance instead.
        type(shape_engine_t) :: blank

        if (cache%owns_tables .and. associated(cache%tables)) then
            call tables_free_s(cache%tables)
            deallocate(cache%tables)
        end if
        nullify(cache%tables)
        cache%owns_tables = .false.

        cache%engine = blank

        if (allocated(cache%f_grid)) deallocate(cache%f_grid)
        if (allocated(cache%fp_grid)) deallocate(cache%fp_grid)
        if (allocated(cache%z)) deallocate(cache%z)
        if (allocated(cache%rho)) deallocate(cache%rho)
        if (allocated(cache%drho_dz)) deallocate(cache%drho_dz)
        if (allocated(cache%radii)) deallocate(cache%radii)
        if (allocated(cache%dr_scratch)) deallocate(cache%dr_scratch)
        if (allocated(cache%radii_d)) deallocate(cache%radii_d)
        if (allocated(cache%dr_dtheta)) deallocate(cache%dr_dtheta)

        cache%n_params = 0_ik
        cache%a2 = 0.0_rk
        cache%z_shift_intrinsic = 0.0_rk
        cache%f_min = 0.0_rk
        cache%beak_ok = .false.
        cache%rho_max = 0.0_rk
        cache%rho_positive = .false.
        cache%g0 = 0.0_rk
        cache%s_opt = 0.0_rk
        cache%g_opt = 0.0_rk
        cache%z_shift_total = 0.0_rk
        cache%r_north = 0.0_rk
        cache%r_south = 0.0_rk
        cache%star_ok = .false.

    end subroutine cache_free_s

    !> Times an intermediate was recomputed since init; 0 on a freed or
    !! out-of-range query. Forwards the engine's counter for the contract's
    !! mandated minimality tests.
    !!
    !! @param[in] cache         Cache to query
    !! @param[in] intermediate  1-based intermediate index
    !! @return                  Recompute count
    pure function cache_recompute_count_f(cache, intermediate) result(count)

        type(cache_t), intent(in) :: cache
        integer(kind = ik), intent(in) :: intermediate
        integer(kind = ikl) :: count

        count = shape_engine_recompute_count_f(cache%engine, intermediate)

    end function cache_recompute_count_f

    !===========================================================================
    ! CACHED COMPUTES
    !===========================================================================

    !> Resolved shape: total z-shift and the analytic pole radii.
    !!
    !! Runs the whole resolve pipeline (#1-#6) and applies the three shape
    !! gates, in producer order: beak singularity (103), interior rho <= 0
    !! (100), star-convexity margin (101). Wrong parameter count (4) and
    !! degenerate c (102) are rejected first, inside `ensure_list_s`.
    !!
    !! NOTE on 103 vs 100: an interior rho <= 0 means f <= 0 there, which the
    !! denser beak scan also sees, so a shape that 1.x rejected with 100 is
    !! rejected here with 103. The cylindrical output, which never runs #4,
    !! still reports 100 for it.
    !!
    !! Every failure zero-fills all three outputs and returns the engine to
    !! cold, so the next call runs from scratch.
    !!
    !! @param[in,out] cache    Initialized cache
    !! @param[in]     params   Parameter vector, length == the cache's n_params
    !! @param[out]    z_shift  Total shift: intrinsic COM shift + chosen origin
    !! @param[out]    r_north  R(0) = c + z_shift
    !! @param[out]    r_south  R(pi) = |-c + z_shift|
    !! @param[out]    status   SHAPE_VALID on success, else the rejecting code
    subroutine cache_shape_s(cache, params, z_shift, r_north, r_south, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south
        integer(kind = ik), intent(out) :: status

        ! Zeroed here once: every rejection below returns without touching them.
        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk

        call ensure_resolved_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        z_shift = cache%z_shift_total
        r_north = cache%r_north
        r_south = cache%r_south

    end subroutine cache_shape_s

    !> Cylindrical rho(z) profile in the COM frame.
    !!
    !! Gates: wrong parameter count (4), degenerate c (102), interior rho <= 0
    !! (100), output size mismatch (105). Beak and star-convexity do NOT gate
    !! here — both exist only to guarantee the R(theta) conversion, and a
    !! beak-marginal shape still has a well-defined cylindrical profile.
    !!
    !! Every failure zero-fills the ENTIRE actual extent of every output and
    !! returns the engine to cold.
    !!
    !! @param[in,out] cache    Initialized cache
    !! @param[in]     params   Parameter vector, length == the cache's n_params
    !! @param[out]    z        Axial coordinate, n_points long
    !! @param[out]    rho      Cylindrical radius, n_points long
    !! @param[out]    drho_dz  Slope, n_points long
    !! @param[out]    z_shift  Intrinsic COM shift baked into z
    !! @param[out]    status   SHAPE_VALID on success, else the rejecting code
    subroutine cache_rho_z_grid_s(cache, params, z, rho, drho_dz, z_shift, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i, n

        z_shift = 0.0_rk

        call ensure_list_s(cache, params, NEEDS_RHO_Z_GRID, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        ! The buffer check follows the engine so that a wrong parameter count
        ! always outranks a buffer mismatch, and so that cache%tables is known
        ! to be associated by the time its resolution is read.
        n = cache%tables%n_points
        if (size(z, kind = ik) /= n .or. size(rho, kind = ik) /= n &
                .or. size(drho_dz, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call shape_engine_invalidate_all_s(cache%engine)
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        if (.not. cache%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            call shape_engine_invalidate_all_s(cache%engine)
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        do i = 1_ik, n
            z(i) = cache%z(i)
            rho(i) = cache%rho(i)
            drho_dz(i) = cache%drho_dz(i)
        end do
        z_shift = cache%z_shift_intrinsic

    end subroutine cache_rho_z_grid_s

    !> R(theta) at the cache's own theta nodes.
    !!
    !! Gate order: the parameter-vector gates of `ensure_resolved_s` (4, 102,
    !! then 103/100/101), then the buffer size (105), then the Newton solve —
    !! a node that misses NR_TOLERANCE rejects the whole grid with 104.
    !!
    !! The result lives in `cache%radii` (#7), so a repeat call on the same
    !! vector copies it out without solving again.
    !!
    !! Every failure zero-fills the ENTIRE actual extent of `radii` and returns
    !! the engine to cold.
    !!
    !! @param[in,out] cache   Initialized cache
    !! @param[in]     params  Parameter vector, length == the cache's n_params
    !! @param[out]    radii   R(theta) at the init thetas, n_theta long
    !! @param[out]    status  SHAPE_VALID on success, else the rejecting code
    subroutine cache_radius_grid_s(cache, params, radii, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i, n

        call ensure_resolved_s(cache, params, status)
        if (status /= SHAPE_VALID) then
            radii = 0.0_rk
            return
        end if

        n = cache%tables%n_theta
        if (size(radii, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call shape_engine_invalidate_all_s(cache%engine)
            radii = 0.0_rk
            return
        end if

        call ensure_list_s(cache, params, NEEDS_RADIUS_GRID, status)
        if (status /= SHAPE_VALID) then
            radii = 0.0_rk
            return
        end if

        do i = 1_ik, n
            radii(i) = cache%radii(i)
        end do

    end subroutine cache_radius_grid_s

    !> R(theta) and dR/dtheta at the cache's own theta nodes.
    !!
    !! Same gates and the same failure semantics as `cache_radius_grid_s`; the
    !! result lives in `cache%radii_d` / `cache%dr_dtheta` (#8). #7 and #8 are
    !! separate intermediates: an output that needs the derivative never reads
    !! the radius-only grid, so the two never share a solve.
    !!
    !! @param[in,out] cache      Initialized cache
    !! @param[in]     params     Parameter vector, length == the cache's n_params
    !! @param[out]    radii      R(theta) at the init thetas, n_theta long
    !! @param[out]    dr_dtheta  dR/dtheta at those thetas, n_theta long
    !! @param[out]    status     SHAPE_VALID on success, else the rejecting code
    subroutine cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i, n

        call ensure_resolved_s(cache, params, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        n = cache%tables%n_theta
        if (size(radii, kind = ik) /= n .or. size(dr_dtheta, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call shape_engine_invalidate_all_s(cache%engine)
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        call ensure_list_s(cache, params, NEEDS_RADIUS_DERIV, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        do i = 1_ik, n
            radii(i) = cache%radii_d(i)
            dr_dtheta(i) = cache%dr_dtheta(i)
        end do

    end subroutine cache_radius_and_derivative_s

    !> R(theta) and dR/dtheta at caller-supplied theta nodes.
    !!
    !! Shares the resolve prefix with the fixed-grid forms but stamps NO
    !! intermediate: the nodes are the caller's, so nothing computed here can be
    !! reused by a later call. Feeding it the cache's own theta nodes reproduces
    !! `cache_radius_and_derivative_s` bit for bit — same bundle, same cosines,
    !! same solver.
    !!
    !! Gate order: the gates of `ensure_resolved_s` (4, 102, 103/100/101), then
    !! the buffer sizes (105), then the theta range (invalid grid), then the
    !! Newton solve (104). Sizes precede the range check for the same reason
    !! they precede everything else in the cached tier: a buffer the caller
    !! cannot fill is a harder error than a node value it can.
    !!
    !! Every failure zero-fills the ENTIRE actual extent of both outputs and
    !! returns the engine to cold.
    !!
    !! @param[in,out] cache      Initialized cache
    !! @param[in]     params     Parameter vector, length == the cache's n_params
    !! @param[in]     thetas     Polar angles in [0, pi]
    !! @param[out]    radii      R(theta), size(thetas) long
    !! @param[out]    dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out]    status     SHAPE_VALID on success, else the rejecting code
    subroutine cache_radius_and_derivative_at_thetas_s(cache, params, thetas, &
            radii, dr_dtheta, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(fos_bundle_t) :: bundle
        integer(kind = ik) :: i, n

        call ensure_resolved_s(cache, params, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        n = size(thetas, kind = ik)
        if (size(radii, kind = ik) /= n .or. size(dr_dtheta, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call shape_engine_invalidate_all_s(cache%engine)
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        do i = 1_ik, n
            if (thetas(i) < 0.0_rk .or. thetas(i) > PI_C) then
                status = SHAPE_ERROR_INVALID_GRID
                call shape_engine_invalidate_all_s(cache%engine)
                call zero_fill_pair_s(radii, dr_dtheta)
                return
            end if
        end do

        bundle = fos_bundle_f(params, cache%z_shift_total, cache%rho_max)
        call solve_thetas_s(bundle, thetas, radii, dr_dtheta, status)

        if (status /= SHAPE_VALID) then
            call shape_engine_invalidate_all_s(cache%engine)
            call zero_fill_pair_s(radii, dr_dtheta)
        end if

    end subroutine cache_radius_and_derivative_at_thetas_s

    !> Neck of the cylindrical profile: the interior rho minimum between the two
    !! largest rho maxima, refined to machine precision.
    !!
    !! Gates: wrong parameter count (4), degenerate c (102), interior rho <= 0
    !! (100). Beak (103) and star-convexity (101) do NOT gate here — the neck is
    !! a property of the rho(z) profile, which a beak-marginal or
    !! non-star-convex shape still has. That asymmetry is the point of this
    !! entry: `rho_neck -> 0` defines the scission line, and scission-adjacent
    !! shapes are exactly the ones the R(theta) gates reject.
    !!
    !! A coarse scan over the cached grid (#5) brackets the neck, then Newton
    !! iteration on f'(u) = 0, bracketed to one grid spacing around the scan
    !! result, refines it. `found` is .false. — with status SHAPE_VALID — for a
    !! shape with fewer than two rho maxima: having no neck is an answer, not an
    !! error.
    !!
    !! Every failure zero-fills both outputs, clears `found`, and returns the
    !! engine to cold.
    !!
    !! @param[in,out] cache     Initialized cache
    !! @param[in]     params    Parameter vector, length == the cache's n_params
    !! @param[out]    z_neck    Neck z-position in the COM frame
    !! @param[out]    rho_neck  Neck radius, reduced units (R0 = 1)
    !! @param[out]    found     .true. iff the profile has a neck
    !! @param[out]    status    SHAPE_VALID on success, else the rejecting code
    subroutine cache_neck_s(cache, params, z_neck, rho_neck, found, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: neck_idx, iter, n
        real(kind = rk) :: c, u, u_lo, u_hi, du, step
        real(kind = rk) :: f_val, fp_val, fpp_val

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call ensure_list_s(cache, params, NEEDS_RHO_Z_GRID, status)
        if (status /= SHAPE_VALID) return

        if (.not. cache%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            call shape_engine_invalidate_all_s(cache%engine)
            return
        end if

        call find_neck_index_s(cache%rho, neck_idx, found)
        if (.not. found) return

        ! Newton refinement of the f-minimum, bracketed to one grid spacing
        ! around the scan result so it cannot escape to a different extremum.
        c = params(1)
        n = cache%tables%n_points
        du = 2.0_rk / real(n - 1_ik, rk)
        u = cache%tables%u(neck_idx)
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
        z_neck = c * u + cache%z_shift_intrinsic

    end subroutine cache_neck_s

    !> Diagnostic: the raw, UNGATED star-convexity optimum of a shape.
    !!
    !! Returns g(s*) = min_s max_i[(z_i + s) drho_dz_i - rho_i] and the total
    !! shift z_shift_intrinsic + s*, WITHOUT applying STAR_CONVEXITY_MARGIN. The
    !! mathematical single-valued (representable) boundary is g(s*) = 0; the
    !! production gate accepts g(s*) <= -STAR_CONVEXITY_MARGIN. This entry
    !! exposes the quantity the margin thresholds, so a shape the R(theta) path
    !! rejects with 101 is still reported on here — that is the whole purpose,
    !! and the reason star-convexity is the one resolve gate this path omits.
    !!
    !! The returned shift is ALWAYS the optimum, intrinsic + s*, never the
    !! branch-selected R(theta) origin (which keeps the COM when the COM is
    !! already well-conditioned).
    !!
    !! Gates: wrong parameter count (4), degenerate c (102), beak singularity
    !! (103), interior rho <= 0 (100) — the resolve gate set minus 101, in the
    !! same order, so a shape that 1.x rejected with 100 is rejected here with
    !! 103 exactly as in `cache_shape_s`.
    !!
    !! Reads the STORED resolve (#6): after a successful `cache_shape_s` on the
    !! same vector this call recomputes nothing.
    !!
    !! Every failure zero-fills both outputs and returns the engine to cold.
    !!
    !! @param[in,out] cache          Initialized cache
    !! @param[in]     params         Parameter vector, length == the cache's n_params
    !! @param[out]    z_shift_total  Intrinsic COM shift plus the optimum s*
    !! @param[out]    g_opt          g(s*), the star-convexity test value there
    !! @param[out]    status         SHAPE_VALID on success, else the rejecting code
    subroutine cache_star_convexity_optimum_s(cache, params, z_shift_total, g_opt, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift_total
        real(kind = rk), intent(out) :: g_opt
        integer(kind = ik), intent(out) :: status

        z_shift_total = 0.0_rk
        g_opt = 0.0_rk

        call ensure_gated_s(cache, params, .false., status)
        if (status /= SHAPE_VALID) return

        z_shift_total = cache%z_shift_intrinsic + cache%s_opt
        g_opt = cache%g_opt

    end subroutine cache_star_convexity_optimum_s

    !===========================================================================
    ! STANDALONE (TIER-1) FORMS
    !===========================================================================
    ! One call in, one answer out. No engine, no cache, no tables handed to the
    ! caller: `standalone_run_s` builds a local `tables_t`, drives the SAME
    ! kernels in the SAME cold order as the cached tier's ensure path, and
    ! everything is released when the form returns.
    !
    ! That shared order is the whole mechanism behind the normative bitwise
    ! rule: at matching (n_points, thetas) a standalone call and a COLD cached
    ! call must agree to the last bit, because they are the same arithmetic in
    ! the same sequence.
    !===========================================================================

    !> R(theta) at caller-supplied theta nodes, standalone.
    !!
    !! Gate order matches `cache_radius_grid_s`: the parameter-vector gates (1,
    !! 102), the grid (invalid grid, covering both `n_points` and the theta
    !! range), the shape gates (103, 100, 101), the buffer size (105), then the
    !! Newton solve (104).
    !!
    !! Every failure zero-fills the ENTIRE actual extent of `radii`.
    !!
    !! @param[in]  params    FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  thetas    Polar angles in [0, pi]
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] radii     R(theta), size(thetas) long
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    subroutine compute_radius_grid_standalone_s(params, thetas, n_points, radii, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st
        type(fos_bundle_t) :: bundle
        real(kind = rk), allocatable :: dr_scratch(:)
        integer(kind = ik) :: n

        call standalone_run_s(params, thetas, n_points, .true., st, status)
        if (status /= SHAPE_VALID) then
            radii = 0.0_rk
            return
        end if

        call sa_gate_s(st, .true., status)
        if (status /= SHAPE_VALID) then
            radii = 0.0_rk
            return
        end if

        n = st%tables%n_theta
        if (size(radii, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            radii = 0.0_rk
            return
        end if

        ! An empty theta set is a well-formed (empty) answer, and `tables%thetas`
        ! is unallocated then — there is nothing to pass to the solver.
        if (n < 1_ik) return

        ! The derivative lands in a write-only sink: the shared solver always
        ! produces one, and this form does not report it.
        allocate(dr_scratch(n))
        bundle = fos_bundle_f(params, st%z_shift_total, st%rho_max)
        call solve_thetas_s(bundle, st%tables%thetas, radii, dr_scratch, status)

        if (status /= SHAPE_VALID) radii = 0.0_rk

    end subroutine compute_radius_grid_standalone_s

    !> R(theta) and dR/dtheta at caller-supplied theta nodes, standalone.
    !!
    !! Same gates and the same failure semantics as
    !! `compute_radius_grid_standalone_s`, and the same solve — so the radii of
    !! the two forms agree bitwise on identical input.
    !!
    !! @param[in]  params     FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  thetas     Polar angles in [0, pi]
    !! @param[in]  n_points   u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] radii      R(theta), size(thetas) long
    !! @param[out] dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out] status     SHAPE_VALID on success, else the rejecting code
    subroutine compute_radius_and_derivative_standalone_s(params, thetas, &
            n_points, radii, dr_dtheta, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st
        type(fos_bundle_t) :: bundle
        integer(kind = ik) :: n

        call standalone_run_s(params, thetas, n_points, .true., st, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        call sa_gate_s(st, .true., status)
        if (status /= SHAPE_VALID) then
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        n = st%tables%n_theta
        if (size(radii, kind = ik) /= n .or. size(dr_dtheta, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call zero_fill_pair_s(radii, dr_dtheta)
            return
        end if

        if (n < 1_ik) return

        bundle = fos_bundle_f(params, st%z_shift_total, st%rho_max)
        call solve_thetas_s(bundle, st%tables%thetas, radii, dr_dtheta, status)

        if (status /= SHAPE_VALID) call zero_fill_pair_s(radii, dr_dtheta)

    end subroutine compute_radius_and_derivative_standalone_s

    !> Resolved shape — total z-shift and the analytic pole radii — standalone.
    !!
    !! Gates: the parameter-vector gates (1, 102), the grid, then the three
    !! shape gates in producer order (103, 100, 101). Needs no theta nodes.
    !!
    !! @param[in]  params    FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_shift   Total shift: intrinsic COM shift + chosen origin
    !! @param[out] r_north   R(0) = c + z_shift
    !! @param[out] r_south   R(pi) = |-c + z_shift|
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    subroutine compute_shape_standalone_s(params, n_points, z_shift, r_north, &
            r_south, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_shift
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st

        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk

        call standalone_run_s(params, NO_THETAS, n_points, .true., st, status)
        if (status /= SHAPE_VALID) return

        call sa_gate_s(st, .true., status)
        if (status /= SHAPE_VALID) return

        z_shift = st%z_shift_total
        r_north = st%r_north
        r_south = st%r_south

    end subroutine compute_shape_standalone_s

    !> Cylindrical rho(z) profile in the COM frame, standalone.
    !!
    !! Gates: the parameter-vector gates (1, 102), the grid, the buffer sizes
    !! (105), then interior rho <= 0 (100). Beak and star-convexity do NOT gate
    !! here, exactly as in `cache_rho_z_grid_s` — the beak scan is not even run.
    !!
    !! Every failure zero-fills the ENTIRE actual extent of every output.
    !!
    !! @param[in]  params    FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z         Axial coordinate, n_points long
    !! @param[out] rho       Cylindrical radius, n_points long
    !! @param[out] drho_dz   Slope, n_points long
    !! @param[out] z_shift   Intrinsic COM shift baked into z
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    subroutine compute_rho_z_grid_standalone_s(params, n_points, z, rho, &
            drho_dz, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st
        integer(kind = ik) :: i, n

        z_shift = 0.0_rk

        call standalone_run_s(params, NO_THETAS, n_points, .false., st, status)
        if (status /= SHAPE_VALID) then
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        n = st%tables%n_points
        if (size(z, kind = ik) /= n .or. size(rho, kind = ik) /= n &
                .or. size(drho_dz, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        if (.not. st%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            call zero_fill_grid_s(z, rho, drho_dz)
            return
        end if

        do i = 1_ik, n
            z(i) = st%z(i)
            rho(i) = st%rho(i)
            drho_dz(i) = st%drho_dz(i)
        end do
        z_shift = st%z_shift_intrinsic

    end subroutine compute_rho_z_grid_standalone_s

    !> Neck of the cylindrical profile, standalone.
    !!
    !! Gates: the parameter-vector gates (1, 102), the grid, then interior
    !! rho <= 0 (100). Beak (103) and star-convexity (101) do NOT gate — the
    !! neck is a property of the rho(z) profile, and scission-adjacent shapes are
    !! exactly the ones the R(theta) gates reject. Same scan-then-Newton
    !! refinement as `cache_neck_s`.
    !!
    !! `found` is .false. with status SHAPE_VALID for a profile with fewer than
    !! two rho maxima: having no neck is an answer, not an error.
    !!
    !! @param[in]  params    FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_neck    Neck z-position in the COM frame
    !! @param[out] rho_neck  Neck radius, reduced units (R0 = 1)
    !! @param[out] found     .true. iff the profile has a neck
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    subroutine compute_neck_standalone_s(params, n_points, z_neck, rho_neck, &
            found, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st
        integer(kind = ik) :: neck_idx, iter, n
        real(kind = rk) :: c, u, u_lo, u_hi, du, step
        real(kind = rk) :: f_val, fp_val, fpp_val

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call standalone_run_s(params, NO_THETAS, n_points, .false., st, status)
        if (status /= SHAPE_VALID) return

        if (.not. st%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            return
        end if

        call find_neck_index_s(st%rho, neck_idx, found)
        if (.not. found) return

        c = params(1)
        n = st%tables%n_points
        du = 2.0_rk / real(n - 1_ik, rk)
        u = st%tables%u(neck_idx)
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
        z_neck = c * u + st%z_shift_intrinsic

    end subroutine compute_neck_standalone_s

    !> Ungated star-convexity optimum, standalone.
    !!
    !! Returns g(s*) and the total shift z_shift_intrinsic + s* WITHOUT applying
    !! STAR_CONVEXITY_MARGIN, so a shape the R(theta) path rejects with 101 is
    !! still reported on. Gates: the parameter-vector gates (1, 102), the grid,
    !! beak (103), interior rho <= 0 (100) — the resolve gate set minus 101, in
    !! the same order as `cache_star_convexity_optimum_s`.
    !!
    !! @param[in]  params         FoS parameters, 1 to FOS_MAX_PARAMS entries
    !! @param[in]  n_points       u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_shift_total  Intrinsic COM shift plus the optimum s*
    !! @param[out] g_opt          g(s*), the star-convexity test value there
    !! @param[out] status         SHAPE_VALID on success, else the rejecting code
    subroutine compute_star_convexity_optimum_standalone_s(params, n_points, &
            z_shift_total, g_opt, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_shift_total
        real(kind = rk), intent(out) :: g_opt
        integer(kind = ik), intent(out) :: status

        type(sa_state_t) :: st

        z_shift_total = 0.0_rk
        g_opt = 0.0_rk

        call standalone_run_s(params, NO_THETAS, n_points, .true., st, status)
        if (status /= SHAPE_VALID) return

        call sa_gate_s(st, .false., status)
        if (status /= SHAPE_VALID) return

        z_shift_total = st%z_shift_intrinsic + st%s_opt
        g_opt = st%g_opt

    end subroutine compute_star_convexity_optimum_standalone_s

    !> The standalone tier's single worker: validate, build local tables, run
    !! the kernels in cold order.
    !!
    !! ## k_max is sized per call — that is the point of tier-1
    !!
    !! `k_max = min((size(params) + 2)/2 + 1, FOS_MAX_K)` is the same formula the
    !! live evaluator `eval_f_s` uses, so a 50-parameter vector gets the 27
    !! Fourier orders it needs instead of the cached tier's fixed 6. This is what
    !! lifts the tier-1 cap from SHAPE_CACHE_MAX_PARAMS to N_max = 50.
    !!
    !! It does NOT break the bitwise rule at the overlap. The formula is
    !! deliberately conservative — a vector of n parameters carries coefficients
    !! only up to a_(n+1), i.e. through order floor((n+1)/2) — and every order
    !! above that has both coefficients zero, so `pair_coefficients_s` marks it
    !! inactive and both `compute_f_grid_s` and `beak_scan_f_min_s` skip it under
    !! the FOS_COEFF_NEGLIGIBLE rule. Two tables differing only in trailing inactive
    !! orders therefore sum the SAME terms in the SAME order. At n <= 8, where a
    !! cold cache exists to compare against, this k_max (2 to 6) and the cache's
    !! fixed FOS_TABLES_K_MAX = 6 are interchangeable bit for bit.
    !!
    !! `size(params) >= 1` also makes k_max >= 2, so `build_tables_s` — which
    !! does not validate k_max — is never handed a degenerate order count.
    !!
    !! @param[in]  params    Incoming parameter vector
    !! @param[in]  thetas    Polar angles in [0, pi]; empty for the theta-less forms
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  resolve   .true. to also run the beak scan (#4) and the origin
    !!                       resolve (#6); .false. stops at the rho(z) grid (#5)
    !! @param[out] st        Per-call state, valid only when status is SHAPE_VALID
    !! @param[out] status    SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS,
    !!                       FOS_ERROR_INVALID_C or SHAPE_ERROR_INVALID_GRID
    subroutine standalone_run_s(params, thetas, n_points, resolve, st, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        logical, intent(in) :: resolve
        type(sa_state_t), intent(out) :: st
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: n, k_max

        n = size(params, kind = ik)

        ! Gate order mirrors `ensure_list_s`: the vector's length is judged
        ! before any code reads params(1).
        if (n > FOS_MAX_PARAMS) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        if (n < 1_ik) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        if (params(1) <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        k_max = min((n + 2_ik) / 2_ik + 1_ik, FOS_MAX_K)

        call build_tables_s(st%tables, n_points, thetas, k_max, status)
        if (status /= SHAPE_VALID) return

        allocate(st%f_grid(n_points), st%fp_grid(n_points))
        allocate(st%z(n_points), st%rho(n_points), st%drho_dz(n_points))

        ! The same kernels the cached tier's ensure path drives, in the same
        ! cold order: #1, #2, #3, [#4], #5, [#6]. The two optional steps are
        ! exactly the ones NEEDS_RHO_Z_GRID omits.
        call compute_a2_s(params, st%a2, status)
        if (status /= SHAPE_VALID) return

        call compute_z_shift_s(params, st%z_shift_intrinsic, status)
        if (status /= SHAPE_VALID) return

        call compute_f_grid_s(st%tables, params, st%f_grid, st%fp_grid)

        if (resolve) then
            call beak_scan_f_min_s(st%tables, params, st%f_min, st%beak_ok)
        end if

        call scale_rho_grid_s(st%tables, params(1), st%z_shift_intrinsic, &
                st%f_grid, st%fp_grid, st%z, st%rho, st%drho_dz, &
                st%rho_max, st%rho_positive)

        if (resolve) then
            call resolve_origin_s(st%z, st%rho, st%drho_dz, params(1), &
                    st%z_shift_intrinsic, st%g0, st%s_opt, st%g_opt, &
                    st%star_ok, st%z_shift_total, st%r_north, st%r_south)
        end if

        status = SHAPE_VALID

    end subroutine standalone_run_s

    !> The standalone twin of `ensure_gated_s`: the shape gates in the normative
    !! producer order, with star-convexity optional.
    !!
    !! Precondition (caller's, unchecked): `st` comes from a `standalone_run_s`
    !! call that returned SHAPE_VALID with `resolve = .true.` — a state built
    !! without the beak scan and the resolve carries no verdicts to gate on.
    !!
    !! @param[in]  st         Resolved per-call state
    !! @param[in]  gate_star  .true. to reject on the star-convexity margin
    !! @param[out] status     SHAPE_VALID, or the rejecting code
    pure subroutine sa_gate_s(st, gate_star, status)

        type(sa_state_t), intent(in) :: st
        logical, intent(in) :: gate_star
        integer(kind = ik), intent(out) :: status

        status = SHAPE_VALID

        if (.not. st%beak_ok) then
            status = FOS_ERROR_BEAK_SINGULARITY
            return
        end if

        if (.not. st%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            return
        end if

        if (gate_star .and. .not. st%star_ok) status = FOS_ERROR_NOT_STAR_CONVEX

    end subroutine sa_gate_s

    !===========================================================================
    ! ENGINE DRIVER (the module's single engine touchpoint)
    !===========================================================================

    !> Presents `params` to the engine and brings the listed intermediates up to
    !! date, in list order.
    !!
    !! Gate order is normative: `shape_engine_begin_s` runs FIRST, so a vector
    !! of the wrong length (including length 0) is rejected before any code
    !! reads `params(1)`. Only then is c checked. Any failure leaves the engine
    !! cold; zero-filling the outputs is the caller's half of the contract.
    !!
    !! `needs` is an explicit index list, not a range: an output that does not
    !! consume an intermediate must leave it neither computed nor stamped.
    !!
    !! @param[in,out] cache   Cache to bring up to date
    !! @param[in]     params  Incoming parameter vector
    !! @param[in]     needs   Intermediate indices, in producer order
    !! @param[out]    status  SHAPE_VALID, or the first rejecting code
    subroutine ensure_list_s(cache, params, needs, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: needs(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: j, id

        call shape_engine_begin_s(cache%engine, params, status)
        if (status /= SHAPE_VALID) return

        ! begin_s accepted the vector, so size(params) == n_params >= 1.
        if (params(1) <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            call shape_engine_invalidate_all_s(cache%engine)
            return
        end if

        do j = 1_ik, size(needs, kind = ik)
            id = needs(j)
            if (.not. shape_engine_needs_f(cache%engine, id)) cycle
            call compute_intermediate_s(cache, params, id, status)
            if (status /= SHAPE_VALID) then
                call shape_engine_invalidate_all_s(cache%engine)
                return
            end if
            call shape_engine_note_computed_s(cache%engine, id)
        end do

    end subroutine ensure_list_s

    !> Brings intermediates #1-#6 up to date and applies all three shape gates.
    !!
    !! The entry point of every output that needs a resolved, R(theta)-
    !! representable shape. A thin name for `ensure_gated_s` with the
    !! star-convexity gate ON — that routine holds the gate order and the
    !! invalidation policy, and is the only place either is written down.
    !!
    !! @param[in,out] cache   Cache to resolve
    !! @param[in]     params  Incoming parameter vector
    !! @param[out]    status  SHAPE_VALID, or the rejecting code
    subroutine ensure_resolved_s(cache, params, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        call ensure_gated_s(cache, params, .true., status)

    end subroutine ensure_resolved_s

    !> Brings intermediates #1-#6 up to date and applies the shape gates, with
    !! the star-convexity gate optional.
    !!
    !! `gate_star = .false.` is the ungated-diagnostic path
    !! (`cache_star_convexity_optimum_s`): it must report g(s*) for exactly the
    !! shapes the margin rejects, so 101 cannot fire. Every other gate, and
    !! their order, is shared with the resolved path — one code path, one
    !! normative order.
    !!
    !! Gate order is normative: beak singularity (103), interior rho <= 0 (100),
    !! star-convexity margin (101), after the wrong-parameter-count (4) and
    !! degenerate-c (102) gates inside `ensure_list_s`. Any rejection leaves the
    !! engine cold; zero-filling the outputs is the caller's half of the
    !! contract.
    !!
    !! @param[in,out] cache       Cache to resolve
    !! @param[in]     params      Incoming parameter vector
    !! @param[in]     gate_star   .true. to reject on the star-convexity margin
    !! @param[out]    status      SHAPE_VALID, or the rejecting code
    subroutine ensure_gated_s(cache, params, gate_star, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: gate_star
        integer(kind = ik), intent(out) :: status

        call ensure_list_s(cache, params, NEEDS_SHAPE, status)
        if (status /= SHAPE_VALID) return

        if (.not. cache%beak_ok) then
            status = FOS_ERROR_BEAK_SINGULARITY
            call shape_engine_invalidate_all_s(cache%engine)
            return
        end if

        if (.not. cache%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            call shape_engine_invalidate_all_s(cache%engine)
            return
        end if

        if (gate_star .and. .not. cache%star_ok) then
            status = FOS_ERROR_NOT_STAR_CONVEX
            call shape_engine_invalidate_all_s(cache%engine)
            return
        end if

    end subroutine ensure_gated_s

    !> Recomputes one intermediate from its producers.
    !!
    !! Preconditions (`ensure_list_s`'s): the cache is initialized, c > C_MIN,
    !! and every producer of `id` is already up to date.
    !!
    !! @param[in,out] cache   Cache holding the intermediate
    !! @param[in]     params  Parameter vector the engine accepted
    !! @param[in]     id      Intermediate index
    !! @param[out]    status  SHAPE_VALID, or the rejecting code
    subroutine compute_intermediate_s(cache, params, id, status)

        type(cache_t), intent(inout) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: id
        integer(kind = ik), intent(out) :: status

        type(fos_bundle_t) :: bundle

        status = SHAPE_VALID

        select case (id)
        case (I_A2)
            call compute_a2_s(params, cache%a2, status)
        case (I_Z_SHIFT)
            call compute_z_shift_s(params, cache%z_shift_intrinsic, status)
        case (I_F_GRID)
            call compute_f_grid_s(cache%tables, params, cache%f_grid, cache%fp_grid)
        case (I_BEAK)
            ! Stores the verdict only; the consuming output owns the 103 gate.
            call beak_scan_f_min_s(cache%tables, params, cache%f_min, cache%beak_ok)
        case (I_RHO_GRID)
            call scale_rho_grid_s(cache%tables, params(1), cache%z_shift_intrinsic, &
                    cache%f_grid, cache%fp_grid, cache%z, cache%rho, cache%drho_dz, &
                    cache%rho_max, cache%rho_positive)
        case (I_RESOLVE)
            ! Infallible by construction; the consuming output owns the 101 gate.
            call resolve_origin_s(cache%z, cache%rho, cache%drho_dz, params(1), &
                    cache%z_shift_intrinsic, cache%g0, cache%s_opt, cache%g_opt, &
                    cache%star_ok, cache%z_shift_total, cache%r_north, cache%r_south)
        case (I_RADIUS_GRID)
            ! Precondition (ensure_resolved_s's): the shape passed 103/100/101,
            ! so #5's rho_max and #6's z_shift_total describe a representable
            ! surface and the analytic bracket encloses every root.
            !
            ! The derivative lands in the write-only `dr_scratch` sink: the
            ! shared solver always produces one, and #7 must not write #8's
            ! buffers.
            if (cache%tables%n_theta > 0_ik) then
                bundle = fos_bundle_f(params, cache%z_shift_total, cache%rho_max)
                call solve_thetas_s(bundle, cache%tables%thetas, cache%radii, &
                        cache%dr_scratch, status)
            end if
        case (I_RADIUS_DERIV)
            if (cache%tables%n_theta > 0_ik) then
                bundle = fos_bundle_f(params, cache%z_shift_total, cache%rho_max)
                call solve_thetas_s(bundle, cache%tables%thetas, cache%radii_d, &
                        cache%dr_dtheta, status)
            end if
        case default
            ! Unreachable today: every needs-list in this module is a named
            ! literal, and none of them contains an index without an arm above.
            ! Failing loudly is deliberate — a silent success here would stamp
            ! an UNCOMPUTED intermediate valid and hand out stale data.
            status = SHAPE_ERROR_INVALID_INIT
        end select

    end subroutine compute_intermediate_s

    !===========================================================================
    ! PRIVATE HELPERS
    !===========================================================================

    !> Allocates the per-shape buffers to the tables' resolution.
    subroutine cache_alloc_buffers_s(cache)

        type(cache_t), intent(inout) :: cache

        integer(kind = ik) :: n_points, n_theta

        n_points = cache%tables%n_points
        n_theta = cache%tables%n_theta

        allocate(cache%f_grid(n_points), cache%fp_grid(n_points))
        allocate(cache%z(n_points), cache%rho(n_points), cache%drho_dz(n_points))

        ! Zero-length for the theta-less tier-1 forms: the R(theta) outputs then
        ! accept only a zero-length buffer, which is the correct contract.
        allocate(cache%radii(n_theta), cache%dr_scratch(n_theta))
        allocate(cache%radii_d(n_theta), cache%dr_dtheta(n_theta))

        cache%f_grid = 0.0_rk
        cache%fp_grid = 0.0_rk
        cache%z = 0.0_rk
        cache%rho = 0.0_rk
        cache%drho_dz = 0.0_rk
        cache%radii = 0.0_rk
        cache%dr_scratch = 0.0_rk
        cache%radii_d = 0.0_rk
        cache%dr_dtheta = 0.0_rk

    end subroutine cache_alloc_buffers_s

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
    !! resolve. The cached tier satisfies the length bound through the engine,
    !! which caps n_params at SHAPE_CACHE_MAX_PARAMS.
    !!
    !! @param[in] params         Parameter vector the engine accepted
    !! @param[in] z_shift_total  Intrinsic COM shift plus the chosen origin (#6)
    !! @param[in] rho_max        Largest rho on the resolved grid (#5)
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
    !! Every R(theta) path routes through here — the two cached intermediates
    !! (#7, #8) and the at-thetas form alike — and that is a correctness
    !! requirement, not tidiness. The contract is that identical thetas in give
    !! bit-identical radii and derivatives out, in every build configuration.
    !! Two mechanisms enforce it, and both are needed:
    !!
    !!   - ONE cos evaluation site. The fixed-grid arms pass `tables%thetas`,
    !!     not a precomputed cosine table, because under `-ffast-math`
    !!     the compiler may lower a vectorizable cos loop (the one in
    !!     `build_tables_s`) to libmvec and a scalar one to the libm call, whose
    !!     results differ by an ulp — and one ulp on cos(theta) moves the Newton
    !!     iterate.
    !!   - ONE machine-code copy, via the `noinline` attribute below.
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
    subroutine solve_thetas_s(bundle, thetas, radii, dr_dtheta, status)

        ! One source loop is not enough: at -O3 with -flto and -finline-limit=1000
        ! this body is inlined into BOTH call sites, and -ffast-math then
        ! contracts the two copies differently — measured as a 1-ulp split in
        ! dR/dtheta while the radii still agreed. `noinline` collapses them back
        ! to one machine-code copy, which is what the bitwise contract actually
        ! needs. The cost is one call per R(theta) batch, not per node.
        !GCC$ ATTRIBUTES noinline :: solve_thetas_s

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

    !> Zeroes the full actual extent of a radius/derivative output pair.
    pure subroutine zero_fill_pair_s(a, b)

        real(kind = rk), intent(out) :: a(:)
        real(kind = rk), intent(out) :: b(:)

        a = 0.0_rk
        b = 0.0_rk

    end subroutine zero_fill_pair_s

    !> Zeroes the full actual extent of the cylindrical outputs, tail included:
    !! an oversized buffer must not keep stale data after a rejection.
    pure subroutine zero_fill_grid_s(z, rho, drho_dz)

        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)

        z = 0.0_rk
        rho = 0.0_rk
        drho_dz = 0.0_rk

    end subroutine zero_fill_grid_s

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
