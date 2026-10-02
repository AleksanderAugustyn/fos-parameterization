/**
 * @file fos_parameterization.h
 * @brief C API for the Fortran Fourier-over-Spheroid (FoS) shape library (v3.0.0).
 *
 * All lengths are in reduced units (R0 = 1); params = [c, a3, a4, a5, ...]
 * (a2 is fixed by the volume constraint).
 *
 * Two tiers of computation:
 *   - One-shot (`fos_param_<compute>`): one call computes a shape. Tables are
 *     built, used and discarded inside the call. n_params is 1 ..
 *     FOS_PARAM_MAX_PARAMS.
 *   - Read-only cache (`fos_param_cache_*`): build a cache once with
 *     fos_param_cache_create(), then request any number of shapes against it.
 *     The cache holds only what depends on the resolution — the u grid, the
 *     Fourier basis and the primary thetas. Nothing derived from a shape
 *     parameter is stored, so no call influences a later one.
 *
 *   The two tiers return the same status and bitwise-identical outputs for any
 *   cache with max_params >= n_params at the same resolution.
 *
 * Parameters:
 *   `max_params`, chosen at cache creation, is the longest vector the cache
 *   accepts (1 .. FOS_PARAM_MAX_PARAMS). Per call, n_params may be anything in
 *   1 .. max_params; missing trailing parameters are zero, and a short vector
 *   gives the same bits as its zero-padded form.
 *
 * Thread safety:
 *   A cache is IMMUTABLE after creation. Every compute takes it `const`, and
 *   any number of threads may compute on one handle concurrently.
 *   fos_param_cache_create() and fos_param_cache_destroy() on a handle must
 *   not race with any other call on it. The one-shot calls hold no state and
 *   are safe from any thread.
 *
 * Lifetime:
 *   Destruction is the caller's job — there is no garbage collection and no
 *   finalizer. fos_param_cache_destroy() accepts NULL.
 *
 * Diagnostics:
 *   There are no message buffers. fos_param_cache_create() returns NULL on
 *   failure and reports the rejecting code through its `status` argument,
 *   which may be NULL to ignore it. The cached computes return the status code
 *   directly; the one-shot calls report through a trailing `int32_t* status`
 *   out-parameter, which may also be NULL. FOS_VALID (0) means success.
 *   fos_param_status_message() maps a code to a fixed, static,
 *   null-terminated string; the returned pointer is owned by the library,
 *   never freed by the caller, and safe to read from any thread.
 *
 * Order of the checks (the first failing one sets the status):
 *   create        max_params < 1 -> 5; max_params > FOS_PARAM_MAX_PARAMS -> 1;
 *                 then the grid -> 3 (n_points below the floor, no thetas, a
 *                 theta outside [0, pi], or a grid the heap cannot satisfy)
 *   cached call   NULL handle -> 2; n_params outside 1 .. max_params -> 4;
 *                 output size -> 105; at-thetas only: no thetas or a theta
 *                 outside [0, pi] -> 3; then the shape codes
 *   one-shot      n_params = 0 -> 4; n_params > FOS_PARAM_MAX_PARAMS -> 1;
 *                 output size -> 105; the grid -> 3; then the shape codes
 *   Shape codes, in order: 102 (c too small), then 103 / 100 / 101 as each
 *   output gates them, then 104.
 *
 * Failure behavior:
 *   On any nonzero status, every numeric output is zero-filled. No state
 *   exists to invalidate; the next call is unaffected.
 *
 * Buffer sizes:
 *   Every cached compute takes the size of its output buffer explicitly. The
 *   size must equal the handle's own extent — n_theta for the radius grids,
 *   n_points for the rho(z) grids, n_thetas for the at-thetas form.
 *   A mismatch returns FOS_ERROR_BUFFER_MISMATCH (105) with the outputs
 *   zero-filled; nothing is written past the caller's stated size. A wrong
 *   size is checked before the shape is judged: a beak-invalid shape passed
 *   with the wrong n_radii reports 105, not 103.
 *
 *   A stated size MUST also be the ACTUAL extent of the buffer you pass. It is
 *   a contract, not a bound the library can verify: on the 105 path the
 *   library zero-fills exactly the stated number of elements, so a size larger
 *   than your real buffer is undefined behavior no check can catch. A negative
 *   output size is not an extent: every cached compute reports it as 105 and
 *   writes nothing. A negative n_params in a cached call is an empty vector
 *   (4). The one-shot calls have no handle extent to mismatch, so a negative
 *   n_params or n_thetas there is a bad C argument and reports
 *   FOS_ERROR_INVALID_INIT (5). Internal marshalling buffers are
 *   heap-allocated, so an outsized size argument can never overflow the stack.
 *
 * Frames:
 *   The R(theta) outputs and fos_param_*shape / *star_convexity_optimum
 *   report the TOTAL z-shift (intrinsic COM shift + any extra star-convexity
 *   shift). The rho(z) grids and the neck report in the COM frame (intrinsic
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
/** Highest n_params either tier accepts, and the highest max_params of a cache. */
#define FOS_PARAM_MAX_PARAMS 50
/** Lowest u-grid resolution a cache or a one-shot call accepts. */
#define FOS_PARAM_N_POINTS_FLOOR 100

/* --- Shared contract status codes (0-99, identical numbers in every
 *     shape-parameterization library). Code 6 is retired and never reused. --- */
#define FOS_VALID                          0
#define FOS_ERROR_TOO_MANY_PARAMS          1
#define FOS_ERROR_CACHE_NOT_INITIALIZED    2
#define FOS_ERROR_INVALID_GRID             3
#define FOS_ERROR_WRONG_PARAM_COUNT        4
#define FOS_ERROR_INVALID_INIT             5

/* --- FoS-specific status codes (>= 100) --- */
#define FOS_ERROR_RHO_NEGATIVE             100
#define FOS_ERROR_NOT_STAR_CONVEX          101
#define FOS_ERROR_INVALID_C                102
#define FOS_ERROR_BEAK_SINGULARITY         103
#define FOS_ERROR_CONVERGENCE              104
#define FOS_ERROR_BUFFER_MISMATCH          105

/** Opaque read-only cache. Created once, shared across threads. */
typedef struct fos_param_cache fos_param_cache_t;

/* ===================================================================== */
/* Diagnostics                                                            */
/* ===================================================================== */

/**
 * Static, null-terminated description of a status code. Never NULL; unknown
 * codes get a fallback string. The pointer is owned by the library.
 */
const char *fos_param_status_message(int32_t status);

/* ===================================================================== */
/* Cache lifecycle                                                        */
/* ===================================================================== */

/**
 * Build a read-only cache.
 *
 * @param max_params  Longest parameter vector the cache will accept,
 *                    1 .. FOS_PARAM_MAX_PARAMS. It sets the table order.
 * @param n_points    u-grid resolution, >= FOS_PARAM_N_POINTS_FLOOR
 * @param thetas      n_theta primary polar angles in [0, pi]
 * @param n_theta     Number of primary thetas, at least 1
 * @param status      Receives 0 or the rejecting code (5, 1 or 3). May be NULL.
 * @return            Handle, or NULL on failure
 */
fos_param_cache_t *fos_param_cache_create(int32_t max_params, int32_t n_points,
                                          const double *thetas, int32_t n_theta,
                                          int32_t *status);

/** Release a cache. NULL-safe. Must not race with any other call on it. */
void fos_param_cache_destroy(fos_param_cache_t *cache);

/* ===================================================================== */
/* Cached computes — const handle, safe to call concurrently              */
/* ===================================================================== */

/** R(theta) at the handle's own thetas. n_radii must equal n_theta. */
int32_t fos_param_cache_radius_grid(const fos_param_cache_t *cache,
                                    const double *params, int32_t n_params,
                                    double *radii, int32_t n_radii);

/** R and dR/dtheta at the handle's own thetas. n_radii must equal n_theta. */
int32_t fos_param_cache_radius_and_derivative(const fos_param_cache_t *cache,
                                              const double *params, int32_t n_params,
                                              double *radii, double *dr_dtheta,
                                              int32_t n_radii);

/**
 * R and dR/dtheta at caller-supplied thetas in [0, pi], at least one. `radii`
 * and `dr_dtheta` are n_thetas long. Pass the concatenation of several grids
 * to evaluate them in one call. Given the handle's own thetas this returns the
 * same bits as fos_param_cache_radius_and_derivative().
 */
int32_t fos_param_cache_radius_and_derivative_at_thetas(
        const fos_param_cache_t *cache, const double *params, int32_t n_params,
        const double *thetas, int32_t n_thetas, double *radii, double *dr_dtheta);

/**
 * Resolve step: TOTAL z-shift (intrinsic COM + star-convexity) and the
 * analytic pole radii r_north = R(0), r_south = R(pi) in that frame.
 */
int32_t fos_param_cache_shape(const fos_param_cache_t *cache, const double *params,
                              int32_t n_params, double *z_shift, double *r_north,
                              double *r_south);

/**
 * Cylindrical profile rho(z) with exact drho/dz in the COM frame; rho = 0 at
 * the poles. n_z must equal the handle's n_points. Star-convexity and the beak
 * gate do NOT apply here, so plotters can render shapes the R(theta)
 * conversion rejects. Status 100 means a SAMPLED interior node has rho <= 0: a
 * gap narrower than the node spacing passes, so this is a resolution-dependent
 * check, not a connectivity verdict.
 */
int32_t fos_param_cache_rho_z_grid(const fos_param_cache_t *cache, const double *params,
                                   int32_t n_params, double *z, double *rho,
                                   double *drho_dz, int32_t n_z, double *z_shift);

/**
 * Cylindrical profile without the rho-positivity gate: never returns 100.
 * rho = 0 and drho_dz = 0 at the two tips and at every node where f <= 0,
 * which for a separated shape is the void between the fragments. For a shape
 * fos_param_cache_rho_z_grid() accepts, the two return the same bits.
 * `z_shift` is the closed-form intrinsic shift; for a separated shape it is
 * NOT the centre of mass of the fragments.
 */
int32_t fos_param_cache_rho_z_grid_unchecked(const fos_param_cache_t *cache,
                                             const double *params, int32_t n_params,
                                             double *z, double *rho, double *drho_dz,
                                             int32_t n_z, double *z_shift);

/**
 * Neck (interior minimum of rho between two maxima), Newton-refined, in the
 * COM frame. A valid shape without a neck returns FOS_VALID with *found = 0
 * and z_neck = rho_neck = 0.
 */
int32_t fos_param_cache_neck(const fos_param_cache_t *cache, const double *params,
                             int32_t n_params, double *z_neck, double *rho_neck,
                             int32_t *found);

/**
 * Star-convexity optimum: the total shift (intrinsic COM + the optimum s*) and
 * g(s*), the test value there. Not gated on the star-convexity margin.
 */
int32_t fos_param_cache_star_convexity_optimum(const fos_param_cache_t *cache,
                                               const double *params, int32_t n_params,
                                               double *z_shift_total, double *g_opt);

/* ===================================================================== */
/* One-shot computes — no handle, trailing nullable status                */
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

/**
 * Cylindrical rho(z) without the rho-positivity gate; all three buffers are
 * n_points long. See fos_param_cache_rho_z_grid_unchecked().
 */
void fos_param_rho_z_grid_unchecked(const double *params, int32_t n_params,
                                    int32_t n_points, double *z, double *rho,
                                    double *drho_dz, double *z_shift, int32_t *status);

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

/** a2 from the volume-conservation constraint. An empty vector is the sphere. */
void fos_param_a2(const double *params, int32_t n_params, double *a2, int32_t *status);

/* ===================================================================== */
/* Raw evaluator — outside the tier rules                                 */
/* ===================================================================== */

/**
 * rho and drho/dz at one z, unvalidated and status-free. `z` is in the SHIFTED
 * frame, so the FoS coordinate is u = (z - z_shift) / c. Degenerate input —
 * an empty vector, c at or below the minimum, or |u| at a tip within roundoff
 * — yields the documented fallback rho = 0, drho_dz = 0 rather than an error.
 * The same fallback holds wherever f(u) <= 0, the void between the fragments
 * of a separated shape included. The parameter vector is zero-extended past
 * its end; there is no length cap.
 */
void fos_param_rho_at_z(const double *params, int32_t n_params, double z,
                        double z_shift, double *rho, double *drho_dz);

#ifdef __cplusplus
}
#endif

#endif /* FOS_PARAMETERIZATION_H */
