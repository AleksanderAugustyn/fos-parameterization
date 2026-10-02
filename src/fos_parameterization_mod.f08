!> Public surface of the Fourier-over-Spheroid (FoS) nuclear shape library.
!!
!! This module is the library's front door and owns the two tiers: `cache_t`,
!! the usage checks, the parameter trimming, the per-call scratch, the shape
!! gates, and every cached and one-shot routine. It owns no arithmetic. Every
!! floating-point kernel lives in `fos_parameterization_workers_mod`; the code
!! here checks, trims, allocates, calls kernels, applies gates and copies
!! results. The only arithmetic written in this file belongs to the raw
!! evaluators (see "Outside the tier rules" below).
!!
!! ## Two tiers
!!
!! 1. **One-shot** — `compute_*_standalone_s`: one call in, one answer out.
!!    The routine builds a local cache, calls the cached routine and lets the
!!    cache go out of scope. Nothing is retained; there is nothing to free.
!! 2. **Read-only cache** — `cache_*_s` on a `cache_t`: the caller builds a
!!    cache once with `cache_init_s` and requests any number of shapes against
!!    it. The cache holds only what depends on the resolution: the u grid, the
!!    Fourier basis on it, the beak-scan basis and the primary thetas. Nothing
!!    derived from a shape parameter is ever stored.
!!
!! One-shot and cached calls return the same status and bitwise-identical
!! outputs for any cache with `max_params >= size(params)` at the same
!! resolution. There is one pipeline; the one-shot forms are wrappers.
!!
!! Physics policy is NOT this library's job. Rejections are mathematics
!! (c > 0, rho > 0) plus the numerical constraints the R(theta) conversion
!! itself needs (beak `F_MIN_THRESHOLD`, `STAR_CONVEXITY_MARGIN`). Filtering
!! mathematically valid but physically meaningless shapes belongs to the
!! consumer.
!!
!! ## Thread and lifetime rules (normative)
!!
!! - **A cache is immutable after `cache_init_s`.** Every compute takes it
!!   `intent(in)`, and every compute is `pure`. Any number of threads may
!!   compute on one cache concurrently.
!! - **`cache_init_s` and `cache_free_s` on a cache must not race with any
!!   other call on it.**
!! - **No hidden state.** The library holds no mutable module-level state.
!! - A `cache_t` has allocatable components only: intrinsic assignment is a
!!   deep copy, and a cache that goes out of scope releases its tables.
!!   `cache_free_s` is available but not mandatory.
!!
!! ## Parameters
!!
!! `params = [c, a3, a4, a5, ...]`; a2 is fixed by the volume constraint.
!! A cache accepts any vector of `1 .. max_params` entries; a one-shot call
!! accepts `1 .. FOS_MAX_PARAMS`. Missing trailing parameters are zero, and
!! trailing zeros are trimmed before any kernel runs, so a short vector and
!! its zero-padded form give bitwise-identical outputs. Inputs must be finite.
!!
!! ## Checks, in order (normative)
!!
!! The first failing check sets the status. Usage codes depend on sizes,
!! handles and grids only, never on parameter values.
!!
!!   cache_init_s   max_params < 1 -> 5; > FOS_MAX_PARAMS -> 1; n_points below
!!                  FOS_N_POINTS_FLOOR -> 3; empty theta set -> 3; a theta
!!                  outside [0, pi] -> 3; unallocatable table -> 3
!!   cached call    uninitialized cache -> 2; size(params) outside
!!                  1..max_params -> 4; output buffer size -> 105; at-thetas
!!                  only: empty theta set or a theta outside [0, pi] -> 3;
!!                  unallocatable scratch -> 3; then the value codes
!!   one-shot call  empty params -> 4; more than FOS_MAX_PARAMS -> 1; output
!!                  buffer size -> 105; grid and thetas -> 3; then the value
!!                  codes
!!
!! Value codes, in pipeline order: degenerate c (102), then the gates each
!! output applies — beak (103), interior rho <= 0 (100), star-convexity (101)
!! — then the Newton solve (104).
!!
!!   radius_grid, radius_and_derivative, at-thetas   103, 100, 101, 104
!!   shape                                           103, 100, 101
!!   star_convexity_optimum                          103, 100
!!   rho_z_grid, neck                                100
!!   rho_z_grid_unchecked                            none
!!
!! On any nonzero status every output argument is zero-filled over its entire
!! actual extent.
!!
!! ## Coordinate systems and shifts
!!
!! FoS defines the shape on u = z/c in [-1, 1]: z = c*u in reduced units
!! (R0 = 1), and rho^2 = f(u)/c. The intrinsic z-shift places the centre of
!! mass at the origin. The R(theta) conversion additionally needs a star-convex
!! origin, so the resolve step keeps the COM origin when it is already
!! well-conditioned and otherwise moves to the star-convexity optimum s*. The
!! R(theta) outputs and `shape` report the TOTAL shift; `rho_z_grid` and `neck`
!! report in the COM frame (intrinsic shift only). Positive means shifted
!! toward +z.
!!
!! ## Outside the tier rules
!!
!! - **Raw evaluators** — `compute_fos_f_and_derivatives_s`,
!!   `get_fos_coefficient_f` and `compute_rho_at_z_s`: unvalidated, no status,
!!   no trimming; degenerate input gives the documented fallback value.
!! - **Diagnostics** — `compute_f_min_standalone_s` and
!!   `compute_conversion_diagnostic_standalone_s`: one-shot only, probe use.
!!   They follow the one-shot check order, but the conversion diagnostic keeps
!!   `z_shift` and `g_opt` when the solve returns 104.
module fos_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_MAX_PARAMS
    use fos_parameterization_workers_mod, only: tables_t, build_tables_s, &
            tables_free_s, fos_bundle_t, fos_bundle_f, solve_thetas_s, &
            refine_neck_s, active_length_f, table_order_f, shifted_origin_f, &
            compute_a2_s, compute_z_shift_s, compute_f_grid_s, &
            beak_scan_f_min_s, scale_rho_grid_s, resolve_origin_s, &
            FOS_MAX_PARAMS, FOS_N_POINTS_FLOOR, FOS_MAX_K, FOS_COEFF_NEGLIGIBLE, &
            FOS_U_TIP_TOL, C_MIN, F_MIN_THRESHOLD, STAR_CONVEXITY_MARGIN, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE, FOS_ERROR_BUFFER_MISMATCH

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Read-only cache — build once, share across threads
    !---------------------------------------------------------------------------
    public :: cache_t
    public :: cache_init_s
    public :: cache_free_s
    public :: cache_max_params_f
    public :: cache_n_points_f
    public :: cache_n_thetas_f
    public :: cache_is_initialized_f
    public :: FOS_N_POINTS_FLOOR

    !---------------------------------------------------------------------------
    ! Cached tier
    !---------------------------------------------------------------------------
    public :: cache_radius_grid_s
    public :: cache_radius_and_derivative_s
    public :: cache_radius_and_derivative_at_thetas_s
    public :: cache_shape_s
    public :: cache_rho_z_grid_s
    public :: cache_rho_z_grid_unchecked_s
    public :: cache_neck_s
    public :: cache_star_convexity_optimum_s

    !---------------------------------------------------------------------------
    ! One-shot tier
    !---------------------------------------------------------------------------
    public :: FOS_MAX_PARAMS
    public :: SHAPE_MAX_PARAMS
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s
    public :: compute_shape_standalone_s
    public :: compute_rho_z_grid_standalone_s
    public :: compute_rho_z_grid_unchecked_standalone_s
    public :: compute_neck_standalone_s
    public :: compute_star_convexity_optimum_standalone_s

    !---------------------------------------------------------------------------
    ! Status-reporting scalar helpers
    !---------------------------------------------------------------------------
    public :: compute_a2_s
    public :: compute_z_shift_s

    !---------------------------------------------------------------------------
    ! Outside the tier rules — diagnostics and raw evaluators (module header)
    !---------------------------------------------------------------------------
    public :: compute_f_min_standalone_s
    public :: compute_conversion_diagnostic_standalone_s
    public :: get_fos_coefficient_f
    public :: compute_fos_f_and_derivatives_s
    public :: compute_rho_at_z_s

    !---------------------------------------------------------------------------
    ! Numerical constants a consumer may need to reason about a rejection
    !---------------------------------------------------------------------------
    !! Owned by the workers module (the single owner of every kernel constant);
    !! re-exported here so the public surface is self-contained.
    public :: C_MIN
    public :: F_MIN_THRESHOLD
    public :: STAR_CONVEXITY_MARGIN

    !---------------------------------------------------------------------------
    ! Status codes
    !---------------------------------------------------------------------------
    !! Two disjoint ranges share one integer space. 0-99 are the shared
    !! shape_core codes, re-exported here so a caller never has to use
    !! shape_core_mod directly; 100+ are FoS-specific and owned by the workers
    !! module. Code 6 is retired and never reused. FOS_VALID is an alias of
    !! SHAPE_VALID, so success compares equal across both libraries.
    public :: SHAPE_VALID
    public :: SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED
    public :: SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT
    public :: SHAPE_ERROR_INVALID_INIT

    public :: FOS_ERROR_RHO_NEGATIVE
    public :: FOS_ERROR_NOT_STAR_CONVEX
    public :: FOS_ERROR_INVALID_C
    public :: FOS_ERROR_BEAK_SINGULARITY
    public :: FOS_ERROR_CONVERGENCE
    public :: FOS_ERROR_BUFFER_MISMATCH

    integer(kind = ik), parameter, public :: FOS_VALID = SHAPE_VALID

    !> Length of every string returned by status_message.
    integer(kind = ik), parameter, public :: STATUS_MESSAGE_LEN = 64_ik
    public :: status_message

    !> L: the longest parameter vector either tier accepts.
    integer(kind = ik), parameter :: FOS_PARAM_LIMIT = min(SHAPE_MAX_PARAMS, FOS_MAX_PARAMS)

    !> The empty theta set of the theta-less one-shot forms.
    real(kind = rk), parameter :: NO_THETAS(0) = [real(kind = rk) ::]

    !> Read-only cache: everything determined by `max_params`, the u-grid
    !! resolution and the primary theta set.
    !!
    !! Opaque: the components are private and may grow. Immutable after
    !! `cache_init_s`, shared freely across threads. Allocatable components
    !! only, so assignment is a deep copy and scope exit releases the tables.
    type :: cache_t
        private
        integer(kind = ik) :: max_params = 0_ik   !! Longest accepted vector
        type(tables_t) :: tables                  !! Trig tables, order from max_params
    end type cache_t

    !> Per-call scratch: everything the pipeline derives from one parameter
    !! vector. Lives for one compute call and is never stored.
    type :: work_t
        real(kind = rk) :: a2 = 0.0_rk
        real(kind = rk) :: z_shift_intrinsic = 0.0_rk

        real(kind = rk), allocatable :: f_grid(:)         !! n_points
        real(kind = rk), allocatable :: fp_grid(:)        !! n_points

        real(kind = rk) :: f_min = 0.0_rk
        logical :: beak_ok = .false.

        real(kind = rk), allocatable :: z(:)              !! n_points
        real(kind = rk), allocatable :: rho(:)            !! n_points
        real(kind = rk), allocatable :: drho_dz(:)        !! n_points
        real(kind = rk) :: rho_max = 0.0_rk
        logical :: rho_positive = .false.

        real(kind = rk) :: g0 = 0.0_rk
        real(kind = rk) :: s_opt = 0.0_rk
        real(kind = rk) :: g_opt = 0.0_rk
        real(kind = rk) :: z_shift_total = 0.0_rk
        real(kind = rk) :: r_north = 0.0_rk
        real(kind = rk) :: r_south = 0.0_rk
        logical :: star_ok = .false.
    end type work_t

contains

    !===========================================================================
    ! STATUS
    !===========================================================================

    !> Fixed diagnostic string for a status code, shared or FoS-specific.
    !!
    !! @param[in] code  Any integer; unrecognised values get a fallback message
    !! @return          Blank-padded message, never empty
    pure function status_message(code) result(msg)
        integer(kind = ik), intent(in) :: code
        character(len = STATUS_MESSAGE_LEN) :: msg
        select case (code)
        case (SHAPE_VALID);                        msg = 'valid'
        case (SHAPE_ERROR_TOO_MANY_PARAMS);        msg = 'too many parameters'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);  msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);           msg = 'invalid grid: n_points, theta count, or theta domain'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);      msg = 'params length outside 1..max_params'
        case (SHAPE_ERROR_INVALID_INIT);           msg = 'invalid init arguments'
        case (FOS_ERROR_RHO_NEGATIVE);             msg = 'rho <= 0 away from the poles'
        case (FOS_ERROR_NOT_STAR_CONVEX);          msg = 'shape not star-convex from any origin'
        case (FOS_ERROR_INVALID_C);                msg = 'elongation c below the minimum'
        case (FOS_ERROR_BEAK_SINGULARITY);         msg = 'f_min below the beak threshold'
        case (FOS_ERROR_CONVERGENCE);              msg = 'iteration did not converge'
        case (FOS_ERROR_BUFFER_MISMATCH);          msg = 'output buffer size mismatch'
        case default;                              msg = 'unknown status code'
        end select
    end function status_message

    !===========================================================================
    ! CACHE LIFECYCLE
    !===========================================================================

    !> Builds a read-only cache.
    !!
    !! `max_params` is the consumer's choice: the longest parameter vector the
    !! cache accepts. It sets the table order, `(max_params + 2)/2 + 1`. Calling
    !! this on a live cache releases the old contents through `intent(out)`.
    !!
    !! @param[out] cache       Ready on success; uninitialized otherwise
    !! @param[in]  max_params  1 <= max_params <= FOS_MAX_PARAMS
    !! @param[in]  n_points    u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[in]  thetas      Primary polar angles in [0, pi], at least one
    !! @param[out] status      SHAPE_VALID, SHAPE_ERROR_INVALID_INIT,
    !!                         SHAPE_ERROR_TOO_MANY_PARAMS or SHAPE_ERROR_INVALID_GRID
    pure subroutine cache_init_s(cache, max_params, n_points, thetas, status)

        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(in) :: max_params
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        ! c is mandatory, so a cache for zero parameters has no shape to serve.
        if (max_params < 1_ik) then
            status = SHAPE_ERROR_INVALID_INIT
            return
        end if

        if (max_params > FOS_PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        if (n_points < FOS_N_POINTS_FLOOR) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        ! A cache with no theta nodes can serve no R(theta) output. The builder
        ! below is the permissive one the theta-less one-shot forms share, so
        ! the floor is stated here.
        if (size(thetas, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call cache_build_s(cache, max_params, n_points, thetas, status)

    end subroutine cache_init_s

    !> Releases the cache. Infallible, and safe on an already-free instance.
    !!
    !! @param[in,out] cache  Reset to the uninitialized state
    pure subroutine cache_free_s(cache)

        type(cache_t), intent(inout) :: cache

        call tables_free_s(cache%tables)
        cache%max_params = 0_ik

    end subroutine cache_free_s

    !> Longest parameter vector the cache accepts; 0 when uninitialized.
    pure function cache_max_params_f(cache) result(max_params)

        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: max_params

        max_params = cache%max_params

    end function cache_max_params_f

    !> u-grid resolution of the cache; 0 when uninitialized.
    pure function cache_n_points_f(cache) result(n_points)

        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n_points

        n_points = cache%tables%n_points

    end function cache_n_points_f

    !> Number of primary thetas of the cache; 0 when uninitialized.
    pure function cache_n_thetas_f(cache) result(n_thetas)

        type(cache_t), intent(in) :: cache
        integer(kind = ik) :: n_thetas

        n_thetas = cache%tables%n_theta

    end function cache_n_thetas_f

    !> .true. after a successful `cache_init_s`, .false. after `cache_free_s`.
    pure function cache_is_initialized_f(cache) result(is_initialized)

        type(cache_t), intent(in) :: cache
        logical :: is_initialized

        is_initialized = cache%tables%initialized

    end function cache_is_initialized_f

    !===========================================================================
    ! CACHED TIER
    !===========================================================================

    !> R(theta) at the cache's primary thetas.
    !!
    !! Shares its solve with `cache_radius_and_derivative_s`, so the radii of the
    !! two outputs agree bitwise.
    !!
    !! @param[in]  cache   Initialized cache
    !! @param[in]  params  1 .. max_params parameters
    !! @param[out] radii   R(theta), cache_n_thetas_f(cache) long
    !! @param[out] status  SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_radius_grid_s(cache, params, radii, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        real(kind = rk), allocatable :: dr_sink(:)
        integer :: alloc_stat

        radii = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= cache%tables%n_theta) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        ! The shared solver always produces the derivative; this output does
        ! not report it, so it lands in a write-only sink.
        allocate(dr_sink(cache%tables%n_theta), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        call resolved_radii_s(cache, params, cache%tables%thetas, radii, dr_sink, &
                status)

    end subroutine cache_radius_grid_s

    !> R(theta) and dR/dtheta at the cache's primary thetas.
    !!
    !! @param[in]  cache      Initialized cache
    !! @param[in]  params     1 .. max_params parameters
    !! @param[out] radii      R(theta), cache_n_thetas_f(cache) long
    !! @param[out] dr_dtheta  dR/dtheta, same length
    !! @param[out] status     SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, &
            status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= cache%tables%n_theta &
                .or. size(dr_dtheta, kind = ik) /= cache%tables%n_theta) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call resolved_radii_s(cache, params, cache%tables%thetas, radii, &
                dr_dtheta, status)

    end subroutine cache_radius_and_derivative_s

    !> R(theta) and dR/dtheta at caller-supplied thetas.
    !!
    !! The contract's "radius at caller-chosen thetas": a consumer folds several
    !! grids into one call by passing their concatenation. Given the cache's own
    !! thetas it returns the same bits as `cache_radius_and_derivative_s` — same
    !! pipeline, same bundle, same solver.
    !!
    !! @param[in]  cache      Initialized cache
    !! @param[in]  params     1 .. max_params parameters
    !! @param[in]  thetas     Polar angles in [0, pi], at least one
    !! @param[out] radii      R(theta), size(thetas) long
    !! @param[out] dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out] status     SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_radius_and_derivative_at_thetas_s(cache, params, thetas, &
            radii, dr_dtheta, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: n

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n = size(thetas, kind = ik)
        if (size(radii, kind = ik) /= n .or. size(dr_dtheta, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call check_thetas_s(thetas, status)
        if (status /= SHAPE_VALID) return

        call resolved_radii_s(cache, params, thetas, radii, dr_dtheta, status)

    end subroutine cache_radius_and_derivative_at_thetas_s

    !> Resolved shape: total z-shift and the analytic pole radii.
    !!
    !! NOTE on 103 vs 100: an interior rho <= 0 means f <= 0 there, which the
    !! denser beak scan also sees, so such a shape is rejected here with 103.
    !! The cylindrical output, which never runs the beak scan, reports 100.
    !!
    !! @param[in]  cache    Initialized cache
    !! @param[in]  params   1 .. max_params parameters
    !! @param[out] z_shift  Total shift: intrinsic COM shift + chosen origin
    !! @param[out] r_north  R(0) = c + z_shift
    !! @param[out] r_south  R(pi) = |-c + z_shift|
    !! @param[out] status   SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_shape_s(cache, params, z_shift, r_north, r_south, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south
        integer(kind = ik), intent(out) :: status

        type(work_t) :: work
        integer(kind = ik) :: n_active

        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .true., work, status)
        if (status /= SHAPE_VALID) return

        call gate_s(work, .true., .true., status)
        if (status /= SHAPE_VALID) return

        z_shift = work%z_shift_total
        r_north = work%r_north
        r_south = work%r_south

    end subroutine cache_shape_s

    !> Cylindrical rho(z) profile in the COM frame.
    !!
    !! Beak and star-convexity do NOT gate here — both exist only to guarantee
    !! the R(theta) conversion, and a beak-marginal shape still has a
    !! well-defined cylindrical profile. Status 100 means a SAMPLED interior
    !! node has rho <= 0: a gap narrower than the node spacing passes, so this
    !! is a resolution-dependent check, not a connectivity verdict.
    !!
    !! @param[in]  cache    Initialized cache
    !! @param[in]  params   1 .. max_params parameters
    !! @param[out] z        Axial coordinate, cache_n_points_f(cache) long
    !! @param[out] rho      Cylindrical radius, same length
    !! @param[out] drho_dz  Slope, same length
    !! @param[out] z_shift  Intrinsic COM shift baked into z
    !! @param[out] status   SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_rho_z_grid_s(cache, params, z, rho, drho_dz, z_shift, &
            status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_core_s(cache, params, .true., z, rho, drho_dz, z_shift, &
                status)

    end subroutine cache_rho_z_grid_s

    !> Cylindrical rho(z) profile without the rho-positivity gate.
    !!
    !! Never returns 100. The profile comes back as the rho-grid stage computes
    !! it: rho = 0 and drho/dz = 0 at the two tips and at every node where
    !! f <= 0, which for a separated shape is the void between the fragments.
    !! For a shape `cache_rho_z_grid_s` accepts, the two forms return the same
    !! bits.
    !!
    !! `z_shift` is the closed-form intrinsic shift. For a separated shape it is
    !! NOT the centre of mass of the fragments: the closed form integrates f
    !! with its sign.
    !!
    !! @param[in]  cache    Initialized cache
    !! @param[in]  params   1 .. max_params parameters
    !! @param[out] z        Axial coordinate, cache_n_points_f(cache) long
    !! @param[out] rho      Cylindrical radius, 0 in the void
    !! @param[out] drho_dz  Slope, 0 in the void
    !! @param[out] z_shift  Intrinsic shift baked into z
    !! @param[out] status   SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_rho_z_grid_unchecked_s(cache, params, z, rho, drho_dz, &
            z_shift, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_core_s(cache, params, .false., z, rho, drho_dz, z_shift, &
                status)

    end subroutine cache_rho_z_grid_unchecked_s

    !> Neck of the cylindrical profile: the interior rho minimum between the two
    !! largest rho maxima, refined to machine precision.
    !!
    !! Beak (103) and star-convexity (101) do NOT gate here — the neck is a
    !! property of the rho(z) profile, and scission-adjacent shapes are exactly
    !! the ones the R(theta) gates reject. `found` is .false. — with status
    !! SHAPE_VALID — for a shape with fewer than two rho maxima: having no neck
    !! is an answer, not an error.
    !!
    !! @param[in]  cache     Initialized cache
    !! @param[in]  params    1 .. max_params parameters
    !! @param[out] z_neck    Neck z-position in the COM frame
    !! @param[out] rho_neck  Neck radius, reduced units (R0 = 1)
    !! @param[out] found     .true. iff the profile has a neck
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_neck_s(cache, params, z_neck, rho_neck, found, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        type(work_t) :: work
        integer(kind = ik) :: n_active

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .false., work, status)
        if (status /= SHAPE_VALID) return

        if (.not. work%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            return
        end if

        call refine_neck_s(cache%tables, params(1:n_active), work%rho, &
                work%z_shift_intrinsic, z_neck, rho_neck, found)

    end subroutine cache_neck_s

    !> Diagnostic: the raw, UNGATED star-convexity optimum of a shape.
    !!
    !! Returns g(s*) = min_s max_i[(z_i + s) drho_dz_i - rho_i] and the total
    !! shift z_shift_intrinsic + s*, WITHOUT applying STAR_CONVEXITY_MARGIN, so a
    !! shape the R(theta) outputs reject with 101 is still reported on. The
    !! returned shift is ALWAYS the optimum, never the branch-selected R(theta)
    !! origin (which keeps the COM when the COM is already well-conditioned).
    !!
    !! @param[in]  cache          Initialized cache
    !! @param[in]  params         1 .. max_params parameters
    !! @param[out] z_shift_total  Intrinsic COM shift plus the optimum s*
    !! @param[out] g_opt          g(s*), the star-convexity test value there
    !! @param[out] status         SHAPE_VALID on success, else the rejecting code
    pure subroutine cache_star_convexity_optimum_s(cache, params, z_shift_total, &
            g_opt, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: z_shift_total
        real(kind = rk), intent(out) :: g_opt
        integer(kind = ik), intent(out) :: status

        type(work_t) :: work
        integer(kind = ik) :: n_active

        z_shift_total = 0.0_rk
        g_opt = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .true., work, status)
        if (status /= SHAPE_VALID) return

        call gate_s(work, .true., .false., status)
        if (status /= SHAPE_VALID) return

        z_shift_total = shifted_origin_f(work%z_shift_intrinsic, work%s_opt)
        g_opt = work%g_opt

    end subroutine cache_star_convexity_optimum_s

    !===========================================================================
    ! ONE-SHOT TIER
    !===========================================================================
    ! Each form does its own usage checks, builds a local cache with
    ! max_params = size(params), and calls the cached routine. There is no
    ! second pipeline, which is what makes one-shot == cached a structural
    ! property rather than a tested coincidence.
    !===========================================================================

    !> R(theta) at caller-supplied theta nodes, one-shot.
    !!
    !! @param[in]  params    1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  thetas    Polar angles in [0, pi], at least one
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] radii     R(theta), size(thetas) long
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_radius_grid_standalone_s(params, thetas, n_points, &
            radii, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: radii(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        radii = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(thetas, kind = ik)) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call cache_init_s(cache, size(params, kind = ik), n_points, thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_grid_s(cache, params, radii, status)

    end subroutine compute_radius_grid_standalone_s

    !> R(theta) and dR/dtheta at caller-supplied theta nodes, one-shot.
    !!
    !! @param[in]  params     1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  thetas     Polar angles in [0, pi], at least one
    !! @param[in]  n_points   u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] radii      R(theta), size(thetas) long
    !! @param[out] dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out] status     SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_radius_and_derivative_standalone_s(params, thetas, &
            n_points, radii, dr_dtheta, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call zero_fill_pair_s(radii, dr_dtheta)

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(thetas, kind = ik) &
                .or. size(dr_dtheta, kind = ik) /= size(thetas, kind = ik)) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call cache_init_s(cache, size(params, kind = ik), n_points, thetas, status)
        if (status /= SHAPE_VALID) return

        call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)

    end subroutine compute_radius_and_derivative_standalone_s

    !> Resolved shape — total z-shift and the analytic pole radii — one-shot.
    !!
    !! @param[in]  params    1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_shift   Total shift: intrinsic COM shift + chosen origin
    !! @param[out] r_north   R(0) = c + z_shift
    !! @param[out] r_south   R(pi) = |-c + z_shift|
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_shape_standalone_s(params, n_points, z_shift, r_north, &
            r_south, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_shift
        real(kind = rk), intent(out) :: r_north
        real(kind = rk), intent(out) :: r_south
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        z_shift = 0.0_rk
        r_north = 0.0_rk
        r_south = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        call cache_build_s(cache, size(params, kind = ik), n_points, NO_THETAS, &
                status)
        if (status /= SHAPE_VALID) return

        call cache_shape_s(cache, params, z_shift, r_north, r_south, status)

    end subroutine compute_shape_standalone_s

    !> Cylindrical rho(z) profile in the COM frame, one-shot.
    !!
    !! @param[in]  params    1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z         Axial coordinate, n_points long
    !! @param[out] rho       Cylindrical radius, n_points long
    !! @param[out] drho_dz   Slope, n_points long
    !! @param[out] z_shift   Intrinsic COM shift baked into z
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_rho_z_grid_standalone_s(params, n_points, z, rho, &
            drho_dz, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_one_shot_s(params, n_points, .true., z, rho, drho_dz, &
                z_shift, status)

    end subroutine compute_rho_z_grid_standalone_s

    !> Cylindrical rho(z) profile without the rho-positivity gate, one-shot.
    !! See `cache_rho_z_grid_unchecked_s` for what the void looks like.
    !!
    !! @param[in]  params    1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z         Axial coordinate, n_points long
    !! @param[out] rho       Cylindrical radius, 0 in the void
    !! @param[out] drho_dz   Slope, 0 in the void
    !! @param[out] z_shift   Intrinsic shift baked into z
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_rho_z_grid_unchecked_standalone_s(params, n_points, z, &
            rho, drho_dz, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        call rho_z_grid_one_shot_s(params, n_points, .false., z, rho, drho_dz, &
                z_shift, status)

    end subroutine compute_rho_z_grid_unchecked_standalone_s

    !> Neck of the cylindrical profile, one-shot.
    !!
    !! @param[in]  params    1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  n_points  u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_neck    Neck z-position in the COM frame
    !! @param[out] rho_neck  Neck radius, reduced units (R0 = 1)
    !! @param[out] found     .true. iff the profile has a neck
    !! @param[out] status    SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_neck_standalone_s(params, n_points, z_neck, rho_neck, &
            found, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_neck
        real(kind = rk), intent(out) :: rho_neck
        logical, intent(out) :: found
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        z_neck = 0.0_rk
        rho_neck = 0.0_rk
        found = .false.

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        call cache_build_s(cache, size(params, kind = ik), n_points, NO_THETAS, &
                status)
        if (status /= SHAPE_VALID) return

        call cache_neck_s(cache, params, z_neck, rho_neck, found, status)

    end subroutine compute_neck_standalone_s

    !> Ungated star-convexity optimum, one-shot.
    !!
    !! @param[in]  params         1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  n_points       u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] z_shift_total  Intrinsic COM shift plus the optimum s*
    !! @param[out] g_opt          g(s*), the star-convexity test value there
    !! @param[out] status         SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_star_convexity_optimum_standalone_s(params, n_points, &
            z_shift_total, g_opt, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: z_shift_total
        real(kind = rk), intent(out) :: g_opt
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        z_shift_total = 0.0_rk
        g_opt = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        call cache_build_s(cache, size(params, kind = ik), n_points, NO_THETAS, &
                status)
        if (status /= SHAPE_VALID) return

        call cache_star_convexity_optimum_s(cache, params, z_shift_total, g_opt, &
                status)

    end subroutine compute_star_convexity_optimum_standalone_s

    !===========================================================================
    ! DIAGNOSTICS — outside the tier rules
    !===========================================================================

    !> Raw beak-quantity diagnostic: the beak scan's f_min with its location,
    !! ungated by beak/rho/star by construction. Probe use; the production
    !! verdict on the same quantity is `compute_shape_standalone_s` (103 iff
    !! f_min <= F_MIN_THRESHOLD).
    !!
    !! @param[in]  params        1 .. FOS_MAX_PARAMS parameters
    !! @param[out] f_min         Smallest f over the 1001-pt clamped scan
    !! @param[out] u_at_min      Scan-grid u of the minimum
    !! @param[out] interior_min  .true. iff strict interior scan minimum
    !! @param[out] status        SHAPE_VALID, SHAPE_ERROR_WRONG_PARAM_COUNT,
    !!                           SHAPE_ERROR_TOO_MANY_PARAMS,
    !!                           SHAPE_ERROR_INVALID_GRID or FOS_ERROR_INVALID_C
    pure subroutine compute_f_min_standalone_s(params, f_min, u_at_min, &
            interior_min, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(out) :: f_min
        real(kind = rk), intent(out) :: u_at_min
        logical, intent(out) :: interior_min
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache
        logical :: beak_ok
        integer(kind = ik) :: n_active

        f_min = 0.0_rk
        u_at_min = 0.0_rk
        interior_min = .false.

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        ! The beak grid has its own fixed resolution; the u grid is not read,
        ! so the smallest legal one is built.
        call cache_build_s(cache, size(params, kind = ik), FOS_N_POINTS_FLOOR, &
                NO_THETAS, status)
        if (status /= SHAPE_VALID) return

        if (params(1) <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        n_active = max(1_ik, active_length_f(params))
        call beak_scan_f_min_s(cache%tables, params(1:n_active), f_min, beak_ok, &
                u_at_min, interior_min)

        status = SHAPE_VALID

    end subroutine compute_f_min_standalone_s

    !> Beak- and star-ungated R(theta) conversion with the resolve quantities
    !! surfaced. Probe use ONLY; production consumers use
    !! `compute_radius_and_derivative_standalone_s`.
    !!
    !! Identical to the production conversion except the beak (103) and star
    !! (101) gates are skipped — of the shape gates only interior rho <= 0
    !! (100) can reject. Newton may legitimately fail to converge on a
    !! non-star-convex shape; FOS_ERROR_CONVERGENCE is probe data, not a bug.
    !!
    !! Failure behavior differs from the tiers by design: when the solve returns
    !! 104 the radii and derivatives are zero-filled, but `z_shift` and `g_opt`
    !! KEEP their resolved values — the probe bins a failed conversion by them.
    !!
    !! @param[in]  params     1 .. FOS_MAX_PARAMS parameters
    !! @param[in]  thetas     Polar angles in [0, pi], at least one
    !! @param[in]  n_points   u-grid resolution, >= FOS_N_POINTS_FLOOR
    !! @param[out] radii      R(theta), size(thetas) long
    !! @param[out] dr_dtheta  dR/dtheta, size(thetas) long
    !! @param[out] z_shift    Total shift of the resolved origin
    !! @param[out] g_opt      g(s*), the star-convexity optimum value
    !! @param[out] status     SHAPE_VALID on success, else the rejecting code
    pure subroutine compute_conversion_diagnostic_standalone_s(params, thetas, &
            n_points, radii, dr_dtheta, z_shift, g_opt, status)

        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        real(kind = rk), intent(out) :: z_shift
        real(kind = rk), intent(out) :: g_opt
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache
        type(work_t) :: work
        type(fos_bundle_t) :: bundle
        integer(kind = ik) :: n_active

        call zero_fill_pair_s(radii, dr_dtheta)
        z_shift = 0.0_rk
        g_opt = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(radii, kind = ik) /= size(thetas, kind = ik) &
                .or. size(dr_dtheta, kind = ik) /= size(thetas, kind = ik)) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call cache_init_s(cache, size(params, kind = ik), n_points, thetas, status)
        if (status /= SHAPE_VALID) return

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .true., work, status)
        if (status /= SHAPE_VALID) return

        call gate_s(work, .false., .false., status)
        if (status /= SHAPE_VALID) return

        z_shift = work%z_shift_total
        g_opt = work%g_opt

        bundle = fos_bundle_f(params(1:n_active), work%z_shift_total, work%rho_max)
        call solve_thetas_s(bundle, cache%tables%thetas, radii, dr_dtheta, status)

        if (status /= SHAPE_VALID) call zero_fill_pair_s(radii, dr_dtheta)

    end subroutine compute_conversion_diagnostic_standalone_s

    !===========================================================================
    ! PIPELINE (private)
    !===========================================================================

    !> Builds a cache without the theta floor.
    !!
    !! The theta-less one-shot forms (shape, rho_z_grid, neck, optimum, f_min)
    !! need no theta nodes; `cache_init_s` adds the one-node floor for everyone
    !! else. Preconditions (caller's, unchecked): 1 <= max_params <= L.
    !!
    !! @param[out] cache       Ready on success; uninitialized otherwise
    !! @param[in]  max_params  Longest accepted vector
    !! @param[in]  n_points    u-grid resolution
    !! @param[in]  thetas      Polar angles in [0, pi]; may be empty
    !! @param[out] status      SHAPE_VALID, or SHAPE_ERROR_INVALID_GRID
    pure subroutine cache_build_s(cache, max_params, n_points, thetas, status)

        type(cache_t), intent(out) :: cache
        integer(kind = ik), intent(in) :: max_params
        integer(kind = ik), intent(in) :: n_points
        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        call build_tables_s(cache%tables, n_points, thetas, &
                table_order_f(max_params), status)
        if (status /= SHAPE_VALID) return

        cache%max_params = max_params

    end subroutine cache_build_s

    !> The two usage checks every cached compute starts with: an uninitialized
    !! cache (2), then a parameter count outside 1..max_params (4).
    pure subroutine check_call_s(cache, params, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        if (.not. cache%tables%initialized) then
            status = SHAPE_ERROR_CACHE_NOT_INITIALIZED
            return
        end if

        if (size(params, kind = ik) < 1_ik &
                .or. size(params, kind = ik) > cache%max_params) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if

        status = SHAPE_VALID

    end subroutine check_call_s

    !> The two parameter-count checks every one-shot form starts with: an empty
    !! vector (4), then one longer than the library limit (1) — never silently
    !! truncated.
    pure subroutine check_one_shot_s(params, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(out) :: status

        if (size(params, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_WRONG_PARAM_COUNT
            return
        end if

        if (size(params, kind = ik) > FOS_PARAM_LIMIT) then
            status = SHAPE_ERROR_TOO_MANY_PARAMS
            return
        end if

        status = SHAPE_VALID

    end subroutine check_one_shot_s

    !> Caller-supplied theta set: at least one node, every node in [0, pi].
    pure subroutine check_thetas_s(thetas, status)

        real(kind = rk), intent(in) :: thetas(:)
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: i

        status = SHAPE_VALID

        if (size(thetas, kind = ik) < 1_ik) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        do i = 1_ik, size(thetas, kind = ik)
            if (thetas(i) < 0.0_rk .or. thetas(i) > PI_C) then
                status = SHAPE_ERROR_INVALID_GRID
                return
            end if
        end do

    end subroutine check_thetas_s

    !> The single pipeline: scratch, the degenerate-c gate, then the stage
    !! kernels in producer order — a2, intrinsic z-shift, f grid, [beak scan],
    !! rho grid, [origin resolve].
    !!
    !! Preconditions (caller's, unchecked): `tables` initialized, `params`
    !! already trimmed and 1 .. table-order-compatible entries long.
    !!
    !! @param[in]  tables   The cache's trig tables
    !! @param[in]  params   Trimmed parameter vector
    !! @param[in]  resolve  .true. to also run the beak scan and the origin
    !!                      resolve; .false. stops at the rho(z) grid
    !! @param[out] work     Per-call state, valid only when status is SHAPE_VALID
    !! @param[out] status   SHAPE_VALID, SHAPE_ERROR_INVALID_GRID (scratch) or
    !!                      FOS_ERROR_INVALID_C
    pure subroutine run_pipeline_s(tables, params, resolve, work, status)

        type(tables_t), intent(in) :: tables
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: resolve
        type(work_t), intent(out) :: work
        integer(kind = ik), intent(out) :: status

        integer(kind = ik) :: n
        integer :: alloc_stat

        ! Heap, never automatic: Release puts automatics on the stack, and
        ! n_points is the caller's.
        n = tables%n_points
        allocate(work%f_grid(n), work%fp_grid(n), work%z(n), work%rho(n), &
                work%drho_dz(n), stat = alloc_stat)
        if (alloc_stat /= 0) then
            status = SHAPE_ERROR_INVALID_GRID
            return
        end if

        if (params(1) <= C_MIN) then
            status = FOS_ERROR_INVALID_C
            return
        end if

        call compute_a2_s(params, work%a2, status)
        if (status /= SHAPE_VALID) return

        call compute_z_shift_s(params, work%z_shift_intrinsic, status)
        if (status /= SHAPE_VALID) return

        call compute_f_grid_s(tables, params, work%f_grid, work%fp_grid)

        if (resolve) then
            call beak_scan_f_min_s(tables, params, work%f_min, work%beak_ok)
        end if

        call scale_rho_grid_s(tables, params(1), work%z_shift_intrinsic, &
                work%f_grid, work%fp_grid, work%z, work%rho, work%drho_dz, &
                work%rho_max, work%rho_positive)

        if (resolve) then
            call resolve_origin_s(work%z, work%rho, work%drho_dz, params(1), &
                    work%z_shift_intrinsic, work%g0, work%s_opt, work%g_opt, &
                    work%star_ok, work%z_shift_total, work%r_north, work%r_south)
        end if

        status = SHAPE_VALID

    end subroutine run_pipeline_s

    !> The shape gates in their normative order: beak singularity (103),
    !! interior rho <= 0 (100), star-convexity margin (101).
    !!
    !! Precondition (caller's, unchecked): `work` comes from `run_pipeline_s`
    !! with `resolve = .true.` whenever `gate_beak` or `gate_star` is set.
    !!
    !! @param[in]  work       Resolved per-call state
    !! @param[in]  gate_beak  .true. to reject on the beak verdict
    !! @param[in]  gate_star  .true. to reject on the star-convexity margin
    !! @param[out] status     SHAPE_VALID, or the rejecting code
    pure subroutine gate_s(work, gate_beak, gate_star, status)

        type(work_t), intent(in) :: work
        logical, intent(in) :: gate_beak
        logical, intent(in) :: gate_star
        integer(kind = ik), intent(out) :: status

        status = SHAPE_VALID

        if (gate_beak .and. .not. work%beak_ok) then
            status = FOS_ERROR_BEAK_SINGULARITY
            return
        end if

        if (.not. work%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            return
        end if

        if (gate_star .and. .not. work%star_ok) status = FOS_ERROR_NOT_STAR_CONVEX

    end subroutine gate_s

    !> Shared back half of the three R(theta) outputs: pipeline, gates, bundle,
    !! solve.
    !!
    !! Preconditions (caller's, checked there): the cache and the parameter
    !! count passed `check_call_s`, `thetas` is valid, and both buffers are
    !! size(thetas) long.
    pure subroutine resolved_radii_s(cache, params, thetas, radii, dr_dtheta, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        real(kind = rk), intent(in) :: thetas(:)
        real(kind = rk), intent(out) :: radii(:)
        real(kind = rk), intent(out) :: dr_dtheta(:)
        integer(kind = ik), intent(out) :: status

        type(work_t) :: work
        type(fos_bundle_t) :: bundle
        integer(kind = ik) :: n_active

        call zero_fill_pair_s(radii, dr_dtheta)

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .true., work, status)
        if (status /= SHAPE_VALID) return

        call gate_s(work, .true., .true., status)
        if (status /= SHAPE_VALID) return

        bundle = fos_bundle_f(params(1:n_active), work%z_shift_total, work%rho_max)
        call solve_thetas_s(bundle, thetas, radii, dr_dtheta, status)

        if (status /= SHAPE_VALID) call zero_fill_pair_s(radii, dr_dtheta)

    end subroutine resolved_radii_s

    !> Shared body of the checked and unchecked cylindrical outputs.
    !!
    !! @param[in] checked  .true. to reject a sampled interior rho <= 0 (100)
    pure subroutine rho_z_grid_core_s(cache, params, checked, z, rho, drho_dz, &
            z_shift, status)

        type(cache_t), intent(in) :: cache
        real(kind = rk), intent(in) :: params(:)
        logical, intent(in) :: checked
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        type(work_t) :: work
        integer(kind = ik) :: i, n, n_active

        call zero_fill_grid_s(z, rho, drho_dz)
        z_shift = 0.0_rk

        call check_call_s(cache, params, status)
        if (status /= SHAPE_VALID) return

        n = cache%tables%n_points
        if (size(z, kind = ik) /= n .or. size(rho, kind = ik) /= n &
                .or. size(drho_dz, kind = ik) /= n) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        n_active = max(1_ik, active_length_f(params))
        call run_pipeline_s(cache%tables, params(1:n_active), .false., work, status)
        if (status /= SHAPE_VALID) return

        if (checked .and. .not. work%rho_positive) then
            status = FOS_ERROR_RHO_NEGATIVE
            return
        end if

        do i = 1_ik, n
            z(i) = work%z(i)
            rho(i) = work%rho(i)
            drho_dz(i) = work%drho_dz(i)
        end do
        z_shift = work%z_shift_intrinsic

    end subroutine rho_z_grid_core_s

    !> Shared body of the two one-shot cylindrical forms.
    pure subroutine rho_z_grid_one_shot_s(params, n_points, checked, z, rho, &
            drho_dz, z_shift, status)

        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: n_points
        logical, intent(in) :: checked
        real(kind = rk), intent(out) :: z(:)
        real(kind = rk), intent(out) :: rho(:)
        real(kind = rk), intent(out) :: drho_dz(:)
        real(kind = rk), intent(out) :: z_shift
        integer(kind = ik), intent(out) :: status

        type(cache_t) :: cache

        call zero_fill_grid_s(z, rho, drho_dz)
        z_shift = 0.0_rk

        call check_one_shot_s(params, status)
        if (status /= SHAPE_VALID) return

        if (size(z, kind = ik) /= n_points .or. size(rho, kind = ik) /= n_points &
                .or. size(drho_dz, kind = ik) /= n_points) then
            status = FOS_ERROR_BUFFER_MISMATCH
            return
        end if

        call cache_build_s(cache, size(params, kind = ik), n_points, NO_THETAS, &
                status)
        if (status /= SHAPE_VALID) return

        call rho_z_grid_core_s(cache, params, checked, z, rho, drho_dz, z_shift, &
                status)

    end subroutine rho_z_grid_one_shot_s

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

    !===========================================================================
    ! RAW EVALUATORS
    !===========================================================================
    ! Unvalidated, status-free, outside the tier rules. See the module header.

    !> Computes rho and optionally drho/dz at a given z-coordinate.
    !!
    !! Raw evaluator: no validation, no status. The z argument is in the SHIFTED
    !! frame (the frame the R(theta) origin sits in), so the FoS parameter is
    !! u = (z - z_shift) / c. Degenerate input — an empty vector, c <= C_MIN, or
    !! |u| at a tip within roundoff — returns the documented fallback rho = 0,
    !! drho_dz = 0 rather than a rejection. The same fallback holds wherever
    !! f(u) <= 0, the void between the fragments of a separated shape included.
    !!
    !! @param[in]  params   FoS parameters: params(1) = c, params(k-1) = a_k, k >= 3
    !! @param[in]  z        Axial coordinate in the shifted frame
    !! @param[in]  z_shift  Shift baked into that frame
    !! @param[out] rho      Cylindrical radius, 0 outside the body
    !! @param[out] drho_dz  Optional slope, 0 where rho is 0
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
        if (abs(u) >= 1.0_rk - FOS_U_TIP_TOL) return

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

    !> a2 from the volume constraint. Private helper of get_fos_coefficient_f;
    !! the status-reporting public form is `compute_a2_s`.
    pure function fos_a2_f(params) result(a2)
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

    end function fos_a2_f

    !> Coefficient a_k of the FoS expansion, read out of a parameter vector.
    !!
    !! Raw evaluator: a_1 is 0 by definition, a_2 comes from the volume
    !! constraint, a_k for k >= 3 lives at params(k - 1), and an index past the
    !! end of the vector reads as 0 (zero-extension, never an error).
    !!
    !! @param[in] params  FoS parameters
    !! @param[in] k       Coefficient index
    !! @return            a_k
    pure function get_fos_coefficient_f(params, k) result(a_k)
        real(kind = rk), intent(in) :: params(:)
        integer(kind = ik), intent(in) :: k
        real(kind = rk) :: a_k
        integer(kind = ik) :: idx

        if (k < 2_ik) then
            a_k = 0.0_rk
        else if (k == 2_ik) then
            a_k = fos_a2_f(params)
        else
            idx = k - 1_ik
            if (idx <= size(params, kind = ik)) then
                a_k = params(idx)
            else
                a_k = 0.0_rk
            end if
        end if

    end function get_fos_coefficient_f

    !> Computes f(u) and its derivatives at an arbitrary u.
    !!
    !! Raw evaluator: live trig, no tables, no validation. The tabled kernel
    !! `compute_f_grid_s` reproduces it node-for-node on the u-grid; this form
    !! is what any diagnostic uses off-grid.
    !!
    !! @param[in]  params  FoS parameters (zero-extended past its end)
    !! @param[in]  u       Reduced axial coordinate in [-1, 1]
    !! @param[out] f       f(u) = 1 - u^2 - sum_k [a_2k cos(w_k u) + a_2k+1 sin(psi_k u)]
    !! @param[out] fp      Optional df/du
    !! @param[out] fpp     Optional d2f/du2
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
        k_max = min((n_params + 2_ik) / 2_ik + 1_ik, FOS_MAX_K)

        sum_f = 0.0_rk
        sum_fp = 0.0_rk
        sum_fpp = 0.0_rk

        do k = 1_ik, k_max
            a_even = get_fos_coefficient_f(params, 2_ik * k)
            a_odd = get_fos_coefficient_f(params, 2_ik * k + 1_ik)
            if (abs(a_even) < FOS_COEFF_NEGLIGIBLE .and. &
                    abs(a_odd) < FOS_COEFF_NEGLIGIBLE) cycle

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

end module fos_parameterization_mod
