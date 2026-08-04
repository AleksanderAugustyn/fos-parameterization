!> Public surface of the Fourier-over-Spheroid (FoS) nuclear shape library.
!!
!! This module is the library's front door. It owns no shape algorithm: every
!! kernel, table, cache and standalone form lives in
!! `fos_parameterization_workers_mod` and is re-exported here, so a consumer
!! never has to name the workers module or `shape_core_mod`. The only code
!! written here is the status-message lookup and three raw evaluators (see
!! "Raw evaluators" below).
!!
!! ## Three tiers
!!
!! 1. **Standalone (tier 1)** — `compute_*_standalone_s`: one call in, one
!!    answer out. Tables are sized per call, which is what lifts the parameter
!!    cap to N_max = `FOS_MAX_PARAMS` (50). Nothing is retained; there is
!!    nothing to free.
!! 2. **Cached (tier 2)** — `cache_*_s` on a `cache_t`: a per-shape working
!!    level that recomputes only what the changed parameters invalidated.
!!    Capped at the engine's `SHAPE_CACHE_MAX_PARAMS`.
!! 3. **Shared tables** — a `tables_t` built once and bound to any number of
!!    caches through `cache_init_shared_s`.
!!
!! Physics policy is NOT this library's job. Rejections are mathematics
!! (c > 0, rho > 0) plus the numerical constraints the R(theta) conversion
!! itself needs (beak `F_MIN_THRESHOLD`, `STAR_CONVEXITY_MARGIN`). Filtering
!! mathematically valid but physically meaningless shapes belongs to the
!! consumer.
!!
!! ## Lifetime and thread rules (normative)
!!
!! - **Never copy-assign a `cache_t`.** Intrinsic assignment copies the owned
!!   `tables` pointer shallowly, so the copy double-frees it. One cache has one
!!   owner and one thread.
!! - **One cache per thread.** A `cache_t` is mutated by every compute call.
!!   `tables_t` is read-only after init and may be shared by any number of
!!   threads' caches.
!! - **`cache_free_s` is mandatory — there is no finalizer.** A cache that owns
!!   its tables leaks them if it goes out of scope or is re-initialized unfreed.
!! - **A shared `tables_t` must be declared `target` and must outlive every
!!   cache bound to it** (`cache_init_shared_s` stores a pointer, not a copy).
!!   Free the caches first, then the tables.
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
!! the engine once, then recomputes only the intermediates its own output
!! needs, in that order. `cache_recompute_count_f` reports the counters.
!!
!! ## Coordinate systems and shifts
!!
!! FoS defines the shape on u = z/c in [-1, 1]: z = c*u in reduced units
!! (R0 = 1), and rho^2 = f(u)/c. The intrinsic z-shift places the centre of
!! mass at the origin. The R(theta) conversion additionally needs a star-convex
!! origin, so the resolve step keeps the COM origin when it is already
!! well-conditioned and otherwise moves to the star-convexity optimum s*. The
!! z_shift every output reports is the TOTAL shift applied to the shape's
!! z-coordinates; positive means shifted toward +z.
!!
!! ## Raw evaluators (outside the tier rules)
!!
!! `compute_fos_f_and_derivatives_s`, `get_fos_coefficient_f` and
!! `compute_rho_at_z_s` are raw, unvalidated evaluators kept public for
!! diagnostics. They are NOT part of the tier surface: no length checking
!! beyond zero-extension of short vectors, no status argument, and degenerate
!! input gives the documented fallback value rather than a rejection. The
!! exact-`n_params` rule of the cached tier does not apply to them.
module fos_parameterization_mod

    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use shape_core_mod, only: SHAPE_VALID, SHAPE_ERROR_TOO_MANY_PARAMS, &
            SHAPE_ERROR_CACHE_NOT_INITIALIZED, SHAPE_ERROR_INVALID_GRID, &
            SHAPE_ERROR_WRONG_PARAM_COUNT, SHAPE_ERROR_INVALID_INIT, &
            SHAPE_ERROR_TABLES_NOT_INITIALIZED
    use fos_parameterization_workers_mod, only: tables_t, tables_init_s, tables_free_s, &
            compute_a2_s, compute_z_shift_s, &
            cache_t, cache_init_s, cache_init_shared_s, cache_free_s, &
            cache_shape_s, cache_rho_z_grid_s, cache_radius_grid_s, &
            cache_radius_and_derivative_s, cache_radius_and_derivative_at_thetas_s, &
            cache_neck_s, cache_star_convexity_optimum_s, cache_recompute_count_f, &
            FOS_MAX_PARAMS, compute_radius_grid_standalone_s, &
            compute_radius_and_derivative_standalone_s, compute_shape_standalone_s, &
            compute_rho_z_grid_standalone_s, compute_neck_standalone_s, &
            compute_star_convexity_optimum_standalone_s, &
            FOS_N_POINTS_FLOOR, FOS_MAX_K, FOS_COEFF_NEGLIGIBLE, FOS_U_TIP_TOL, &
            C_MIN, F_MIN_THRESHOLD, STAR_CONVEXITY_MARGIN, &
            FOS_ERROR_RHO_NEGATIVE, FOS_ERROR_NOT_STAR_CONVEX, FOS_ERROR_INVALID_C, &
            FOS_ERROR_BEAK_SINGULARITY, FOS_ERROR_CONVERGENCE, FOS_ERROR_BUFFER_MISMATCH

    implicit none

    private

    !---------------------------------------------------------------------------
    ! Shared tables — build once, bind to many caches
    !---------------------------------------------------------------------------
    !! A `tables_t` is immutable after `tables_init_s` and therefore safe to
    !! share across threads. A tables object handed to `cache_init_shared_s`
    !! MUST be declared with the `target` attribute and MUST outlive every cache
    !! bound to it: the cache stores a pointer, not a copy. Free the caches
    !! first, then the tables.
    public :: tables_t
    public :: tables_init_s
    public :: tables_free_s
    public :: FOS_N_POINTS_FLOOR

    !---------------------------------------------------------------------------
    ! Cached tier — one cache, one owner, one thread
    !---------------------------------------------------------------------------
    !! Never copy-assign a cache (double-free of the owned tables pointer); one
    !! cache per thread; `cache_free_s` is mandatory, there is no finalizer.
    !! `cache_init_shared_s` additionally requires that the caller's tables be
    !! declared `target` and outlive every cache bound to it.
    public :: cache_t
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

    !---------------------------------------------------------------------------
    ! Standalone (tier-1) forms — no engine, no cache, N_max = FOS_MAX_PARAMS
    !---------------------------------------------------------------------------
    public :: FOS_MAX_PARAMS
    public :: compute_radius_grid_standalone_s
    public :: compute_radius_and_derivative_standalone_s
    public :: compute_shape_standalone_s
    public :: compute_rho_z_grid_standalone_s
    public :: compute_neck_standalone_s
    public :: compute_star_convexity_optimum_standalone_s

    !---------------------------------------------------------------------------
    ! Status-reporting scalar helpers
    !---------------------------------------------------------------------------
    public :: compute_a2_s
    public :: compute_z_shift_s

    !---------------------------------------------------------------------------
    ! Raw evaluators — outside the tier rules (see the module header)
    !---------------------------------------------------------------------------
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
    !! module. FOS_VALID is an alias of SHAPE_VALID, so success compares equal
    !! across both libraries.
    public :: SHAPE_VALID
    public :: SHAPE_ERROR_TOO_MANY_PARAMS
    public :: SHAPE_ERROR_CACHE_NOT_INITIALIZED
    public :: SHAPE_ERROR_INVALID_GRID
    public :: SHAPE_ERROR_WRONG_PARAM_COUNT
    public :: SHAPE_ERROR_INVALID_INIT
    public :: SHAPE_ERROR_TABLES_NOT_INITIALIZED

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
        case (SHAPE_ERROR_TOO_MANY_PARAMS);        msg = 'too many parameters for this tier'
        case (SHAPE_ERROR_CACHE_NOT_INITIALIZED);  msg = 'cache not initialized'
        case (SHAPE_ERROR_INVALID_GRID);           msg = 'invalid grid: n_points, theta count, or theta domain'
        case (SHAPE_ERROR_WRONG_PARAM_COUNT);      msg = 'params length differs from n_params'
        case (SHAPE_ERROR_INVALID_INIT);           msg = 'invalid init arguments'
        case (SHAPE_ERROR_TABLES_NOT_INITIALIZED); msg = 'tables not initialized'
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
    ! RAW EVALUATORS
    !===========================================================================
    ! Unvalidated, status-free, outside the tier rules. See the module header.

    !> Computes rho and optionally drho/dz at a given z-coordinate.
    !!
    !! Raw evaluator: no validation, no status. The z argument is in the SHIFTED
    !! frame (the frame the R(theta) origin sits in), so the FoS parameter is
    !! u = (z - z_shift) / c. Degenerate input — an empty vector, c <= C_MIN, or
    !! |u| at a tip within roundoff — returns the documented fallback rho = 0,
    !! drho_dz = 0 rather than a rejection.
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
    !! is what the Newton solver and any diagnostic use off-grid.
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
