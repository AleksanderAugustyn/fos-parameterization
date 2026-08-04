/**
 * @file fos_parameterization.h
 * @brief C API for the Fortran Fourier-over-Spheroid (FoS) shape library (v2.0.0).
 *
 * All lengths are in reduced units (R0 = 1); params = [c, a3, a4, a5, ...]
 * (a2 is fixed by the volume constraint).
 *
 * Two handle types, from shared to per-shape:
 *   - tables handle (`fos_param_tables_create`) — everything determined by the
 *     u-grid resolution `n_points` and the primary theta set. Built once,
 *     shared read-only.
 *   - cache handle (`fos_param_cache_create*`) — per-shape working state: the
 *     recompute engine and every per-shape buffer. Mutated by every compute
 *     call.
 *
 * Two tiers of computation:
 *   - Cached (`fos_param_cache_*`): create a cache once, call it for many
 *     shapes; only the intermediates invalidated by the changed parameters are
 *     recomputed. `n_params` must equal the cache's own n_params EXACTLY
 *     (otherwise FOS_ERROR_WRONG_PARAM_COUNT), and n_params is capped at
 *     FOS_PARAM_CACHE_MAX_PARAMS.
 *   - Flat tier-1 (`fos_param_<compute>`): one-off calls that build, use and
 *     discard their own tables. n_params is capped at FOS_PARAM_MAX_PARAMS.
 *
 * Thread safety:
 *   A cache handle is THREAD-CONFINED. It is mutated by every compute call;
 *   concurrent use of one cache from more than one thread is undefined. Give
 *   every thread its own cache. A tables handle is immutable after creation
 *   and may be shared across threads for concurrent reads, including as the
 *   backing tables of per-thread caches created with
 *   fos_param_cache_create_shared(). The flat tier-1 calls hold no state and
 *   are safe from any thread.
 *
 * Lifetime:
 *   Destruction is MANDATORY — there is no garbage collection and no
 *   finalizer. A tables handle passed to fos_param_cache_create_shared() MUST
 *   outlive every cache built from it: the cache holds a reference, not a
 *   copy. Destroy order: caches first, their tables last. Both destroy
 *   functions accept NULL.
 *
 * Diagnostics:
 *   There are no message buffers. `_create` functions return NULL on failure.
 *   The cached computes return the status code directly; the flat tier-1 calls
 *   report through a trailing `int32_t* status` out-parameter, which may be
 *   NULL to ignore it. FOS_VALID (0) means success.
 *   fos_param_status_message() maps a code to a fixed, static,
 *   null-terminated string; the returned pointer is owned by the library,
 *   never freed by the caller, and safe to read from any thread.
 *
 * Failure behavior:
 *   On any nonzero status, every numeric output is zero-filled and a cache is
 *   returned to a cold state — the next call recomputes from scratch. No
 *   partially updated results are ever visible. A NULL cache handle passed to
 *   a compute function returns FOS_ERROR_CACHE_NOT_INITIALIZED (2).
 *
 * Buffer sizes:
 *   Every cached compute takes the size of its output buffer explicitly. The
 *   size must equal the handle's own extent — n_theta for the radius grids,
 *   n_points for the rho(z) grid, size(thetas) for the at-thetas form.
 *   A mismatch returns FOS_ERROR_BUFFER_MISMATCH (105) with the outputs
 *   zero-filled; nothing is written past the caller's stated size.
 *
 *   A stated size MUST also be the ACTUAL extent of the buffer you pass. It is
 *   a contract, not a bound the library can verify: on the 105 path the
 *   library zero-fills exactly the stated number of elements, so a size larger
 *   than your real buffer is undefined behavior no check can catch. What 105
 *   guarantees is that a size which is merely wrong — but honest about your
 *   own memory — is reported rather than computed on. Negative sizes are read
 *   as zero and therefore also report 105. Internal marshalling buffers are
 *   heap-allocated, so an outsized size argument can never overflow the stack.
 *
 * Frames:
 *   The R(theta) outputs and fos_param_*shape / *star_convexity_optimum
 *   report the TOTAL z-shift (intrinsic COM shift + any extra star-convexity
 *   shift). The rho(z) grid and the neck report in the COM frame (intrinsic
 *   shift only), which is also what fos_param_z_shift returns.
 *
 * Precondition — finite input:
 *   `params` and `thetas` must be finite. Non-finite input is undefined
 *   behavior: the library cannot detect NaN under fast-math, so validation
 *   comparisons silently pass and the call returns FOS_VALID (0) with NaN
 *   outputs. Screen inputs before calling.
 */

#ifndef FOS_PARAMETERIZATION_H
#define FOS_PARAMETERIZATION_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* --- Limits --- */
/** Highest n_params the flat tier-1 calls accept. */
#define FOS_PARAM_MAX_PARAMS 50
/** Highest n_params a cache (cached tier) accepts. */
#define FOS_PARAM_CACHE_MAX_PARAMS 8
/** Lowest u-grid resolution a tables/cache handle may be built for. */
#define FOS_PARAM_N_POINTS_FLOOR 100

/* --- Shared contract status codes (0-99, identical numbers in every
 *     shape-parameterization library) --- */
#define FOS_VALID                          0
#define FOS_ERROR_TOO_MANY_PARAMS          1
#define FOS_ERROR_CACHE_NOT_INITIALIZED    2
#define FOS_ERROR_INVALID_GRID             3
#define FOS_ERROR_WRONG_PARAM_COUNT        4
#define FOS_ERROR_INVALID_INIT             5
#define FOS_ERROR_TABLES_NOT_INITIALIZED   6

/* --- FoS-specific status codes (>= 100) --- */
#define FOS_ERROR_RHO_NEGATIVE             100
#define FOS_ERROR_NOT_STAR_CONVEX          101
#define FOS_ERROR_INVALID_C                102
#define FOS_ERROR_BEAK_SINGULARITY         103
#define FOS_ERROR_CONVERGENCE              104
#define FOS_ERROR_BUFFER_MISMATCH          105

/* ===================================================================== */
/* Diagnostics                                                            */
/* ===================================================================== */

/**
 * Static, null-terminated description of a status code. Never NULL; unknown
 * codes get a fallback string. The pointer is owned by the library.
 */
const char *fos_param_status_message(int32_t status);

/* ===================================================================== */
/* Handle lifecycle                                                       */
/* ===================================================================== */

/**
 * Build the shared, immutable tables level.
 *
 * @param n_points  u-grid resolution, >= FOS_PARAM_N_POINTS_FLOOR
 * @param thetas    n_theta polar angles in [0, pi], at least one
 * @param n_theta   Number of primary thetas
 * @return          Handle, or NULL on failure
 */
void *fos_param_tables_create(int32_t n_points, const double *thetas, int32_t n_theta);

/** Release a tables handle. NULL-safe. Every cache bound to it must be gone. */
void fos_param_tables_destroy(void *tables);

/**
 * Create a cache owning private tables built to the given resolution and
 * thetas.
 *
 * @param n_params  1 .. FOS_PARAM_CACHE_MAX_PARAMS
 * @param n_theta   Number of thetas, at least 1 (same floor as
 *                  fos_param_tables_create)
 * @return          Handle, or NULL on failure
 */
void *fos_param_cache_create(int32_t n_params, int32_t n_points,
                             const double *thetas, int32_t n_theta);

/**
 * Create a cache bound to caller-owned shared tables, which must outlive it.
 * The cache stores a reference, not a copy.
 *
 * @return  Handle, or NULL on failure (including a NULL tables handle)
 */
void *fos_param_cache_create_shared(void *tables, int32_t n_params);

/** Release a cache. NULL-safe. Shared tables are left to their owner. */
void fos_param_cache_destroy(void *cache);

/* ===================================================================== */
/* Cached computes (tier 2)                                               */
/* ===================================================================== */

/** R(theta) at the handle's own thetas. n_radii must equal n_theta. */
int32_t fos_param_cache_radius_grid(void *cache, const double *params, int32_t n_params,
                                    double *radii, int32_t n_radii);

/** R and dR/dtheta at the handle's own thetas. n_radii must equal n_theta. */
int32_t fos_param_cache_radius_and_derivative(void *cache, const double *params,
                                              int32_t n_params, double *radii,
                                              double *dr_dtheta, int32_t n_radii);

/**
 * R and dR/dtheta at caller-supplied thetas in [0, pi] (uncached evaluation
 * against the cache's resolved shape). `radii` and `dr_dtheta` are n_thetas
 * long.
 */
int32_t fos_param_cache_radius_and_derivative_at_thetas(
        void *cache, const double *params, int32_t n_params,
        const double *thetas, int32_t n_thetas, double *radii, double *dr_dtheta);

/**
 * Per-shape resolve step: TOTAL z-shift (intrinsic COM + star-convexity) and
 * the analytic pole radii r_north = R(0), r_south = R(pi) in that frame.
 */
int32_t fos_param_cache_shape(void *cache, const double *params, int32_t n_params,
                              double *z_shift, double *r_north, double *r_south);

/**
 * Cylindrical profile rho(z) with exact drho/dz in the COM frame; rho = 0 at
 * the poles. n_z must equal the handle's n_points. Star-convexity and the beak
 * gate do NOT apply here, so plotters can render shapes the R(theta)
 * conversion rejects.
 */
int32_t fos_param_cache_rho_z_grid(void *cache, const double *params, int32_t n_params,
                                   double *z, double *rho, double *drho_dz,
                                   int32_t n_z, double *z_shift);

/**
 * Neck (interior minimum of rho between two maxima), Newton-refined, in the
 * COM frame. A valid shape without a neck returns FOS_VALID with *found = 0
 * and z_neck = rho_neck = 0.
 */
int32_t fos_param_cache_neck(void *cache, const double *params, int32_t n_params,
                             double *z_neck, double *rho_neck, int32_t *found);

/**
 * Star-convexity optimum: the total shift (intrinsic COM + the optimum s*) and
 * g(s*), the test value there.
 */
int32_t fos_param_cache_star_convexity_optimum(void *cache, const double *params,
                                               int32_t n_params,
                                               double *z_shift_total, double *g_opt);

/* ===================================================================== */
/* Flat tier-1 computes — no handle, trailing nullable status              */
/* ===================================================================== */

/** R(theta) at caller thetas; `radii` is n_thetas long. */
void fos_param_radius_grid(const double *params, int32_t n_params,
                           const double *thetas, int32_t n_thetas,
                           int32_t n_points, double *radii, int32_t *status);

/** R and dR/dtheta at caller thetas; both buffers are n_thetas long. */
void fos_param_radius_and_derivative(const double *params, int32_t n_params,
                                     const double *thetas, int32_t n_thetas,
                                     int32_t n_points, double *radii,
                                     double *dr_dtheta, int32_t *status);

/** Resolve step: TOTAL z-shift and the analytic pole radii. */
void fos_param_shape(const double *params, int32_t n_params, int32_t n_points,
                     double *z_shift, double *r_north, double *r_south,
                     int32_t *status);

/** Cylindrical rho(z) in the COM frame; all three buffers are n_points long. */
void fos_param_rho_z_grid(const double *params, int32_t n_params, int32_t n_points,
                          double *z, double *rho, double *drho_dz,
                          double *z_shift, int32_t *status);

/** Neck in the COM frame; *found is 0 or 1. */
void fos_param_neck(const double *params, int32_t n_params, int32_t n_points,
                    double *z_neck, double *rho_neck, int32_t *found, int32_t *status);

/** Star-convexity optimum: total shift and g(s*). */
void fos_param_star_convexity_optimum(const double *params, int32_t n_params,
                                      int32_t n_points, double *z_shift_total,
                                      double *g_opt, int32_t *status);

/** Intrinsic COM z-shift (closed form). */
void fos_param_z_shift(const double *params, int32_t n_params,
                       double *z_shift, int32_t *status);

/** a2 from the volume-conservation constraint. */
void fos_param_a2(const double *params, int32_t n_params, double *a2, int32_t *status);

/* ===================================================================== */
/* Raw evaluator — outside the tier rules                                 */
/* ===================================================================== */

/**
 * rho and drho/dz at one z, unvalidated and status-free. `z` is in the SHIFTED
 * frame, so the FoS coordinate is u = (z - z_shift) / c. Degenerate input —
 * an empty vector, c at or below the minimum, or |u| at a tip within roundoff
 * — yields the documented fallback rho = 0, drho_dz = 0 rather than an error.
 * The parameter vector is zero-extended past its end; there is no length cap.
 */
void fos_param_rho_at_z(const double *params, int32_t n_params, double z,
                        double z_shift, double *rho, double *drho_dz);

#ifdef __cplusplus
}
#endif

#endif /* FOS_PARAMETERIZATION_H */
