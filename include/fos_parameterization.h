/**
 * @file fos_parameterization.h
 * @brief C API for the Fortran Fourier-over-Spheroid (FoS) shape library (v1.1.0).
 *
 * Standalone functions only: FoS has no shape-independent precompute, so
 * there is no cache/handle tier. All lengths are in reduced units (R0 = 1);
 * params = [c, a3, a4, a5, ...] (a2 is fixed by the volume constraint).
 *
 * Frames: fos_compute_radius_grid returns the TOTAL z-shift baked into the
 * R(theta) grid (intrinsic COM shift + any extra star-convexity shift).
 * fos_compute_rho_profile and fos_compute_neck report in the COM frame
 * (intrinsic shift applied); fos_z_shift returns the intrinsic shift only.
 *
 * Status-returning functions: 0 (FOS_VALID) means success. The message
 * buffer is always null-terminated (truncated to fit). Recommended size: 256.
 * On any failure, numeric output buffers are zero-filled.
 *
 * FOS_ERROR_INVALID_ARGUMENTS (5) is C-API-level only (bad n_grid / n_z);
 * codes 100+ mirror the Fortran FOS_* parameters exactly.
 *
 * Internal rho(z) resolution: the library floors it at 100 nodes, so every
 * call here uses max(100, requested). 1.x used the request verbatim; below
 * the floor the resolved z-shift therefore moves slightly.
 *
 * Thread safety: the library is stateless; all functions are safe to call
 * concurrently.
 */

#ifndef FOS_PARAMETERIZATION_H
#define FOS_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Status codes (100+ mirror the Fortran FOS_* parameters) --- */
#define FOS_VALID                    0
#define FOS_ERROR_INVALID_ARGUMENTS  5
#define FOS_ERROR_RHO_NEGATIVE       100
#define FOS_ERROR_NOT_STAR_CONVEX    101
#define FOS_ERROR_INVALID_C          102
#define FOS_ERROR_BEAK_SINGULARITY   103
#define FOS_ERROR_CONVERGENCE        104
#define FOS_ERROR_BUFFER_MISMATCH    105

/**
 * Validate the shape and compute R(theta) on n_grid uniform points over [0, pi].
 *
 * @param params           FoS parameters [c, a3, a4, ...]
 * @param n_params         Number of parameters (>= 1)
 * @param n_grid           Number of theta grid points (>= 2)
 * @param radii            Output buffer, n_grid doubles
 * @param z_shift          Output: total z-shift baked into the grid
 * @param message_buf_len  Size of message_buf including null terminator
 * @param message_buf      Buffer for the validation message (empty on success)
 * @return                 FOS_VALID (0) on success, error code otherwise
 */
int fos_compute_radius_grid(
        const double* params, int n_params, int n_grid,
        double* radii, double* z_shift,
        int message_buf_len, char* message_buf);

/**
 * Compute the cylindrical profile rho(z) with exact drho/dz in the COM frame.
 * z spans [-c + z_shift_intrinsic, c + z_shift_intrinsic] on n_z uniform
 * points; rho = 0 at the poles. Succeeds for any c > 0 — no star-convexity
 * requirement, so plotters can render shapes the R(theta) conversion rejects.
 *
 * @param n_z  Number of profile points (>= 2)
 */
int fos_compute_rho_profile(
        const double* params, int n_params, int n_z,
        double* z, double* rho, double* drho_dz,
        int message_buf_len, char* message_buf);

/**
 * Find the neck (interior minimum of rho between two maxima), Newton-refined.
 * Positions in the COM frame. A valid shape without a neck returns FOS_VALID
 * with *found = 0 and z_neck = rho_neck = 0.
 */
int fos_compute_neck(
        const double* params, int n_params,
        double* z_neck, double* rho_neck, int* found);

/**
 * Per-shape resolve step: validity, TOTAL z-shift (intrinsic COM +
 * star-convexity), and analytic pole radii r_north = R(0), r_south = R(pi)
 * in the shifted frame. Feed the resulting z_shift to
 * fos_compute_radius_and_derivative_at_thetas for any number of theta sets.
 * On any failure all numeric outputs are zero-filled.
 *
 * @param n_rho_grid  Internal rho(z) grid size (>= 2, raised to the 100 floor)
 */
int fos_compute_shape(
        const double* params, int n_params, int n_rho_grid,
        double* z_shift, double* r_north, double* r_south,
        int message_buf_len, char* message_buf);

/**
 * Batch R(theta) and analytic dR/dtheta at caller-supplied thetas in [0, pi].
 * z_shift must come from fos_compute_shape. No message buffer.
 *
 * This call now VALIDATES the shape, which 1.x did not: it builds an internal
 * rho(z) grid for the Newton bracket bound, and that grid's gates propagate.
 * So on top of the withdrawn 1.x unit-sphere fallback for degenerate params
 * (an empty vector or c <= C_MIN is now FOS_ERROR_INVALID_C), a pinched
 * interior returns FOS_ERROR_RHO_NEGATIVE (100) and more than 50 parameters
 * returns SHAPE_ERROR_TOO_MANY_PARAMS (1) — each with the outputs zero-filled,
 * where 1.x evaluated unconditionally and returned numbers.
 *
 * @param thetas     n_thetas angles in [0, pi]
 * @param radii      Output buffer, n_thetas doubles
 * @param dr_dtheta  Output buffer, n_thetas doubles
 * @return           FOS_VALID, FOS_ERROR_INVALID_ARGUMENTS if n_thetas < 1,
 *                   or the rejecting status code
 */
int fos_compute_radius_and_derivative_at_thetas(
        const double* params, int n_params,
        const double* thetas, int n_thetas, double z_shift,
        double* radii, double* dr_dtheta);

/**
 * Intrinsic COM z-shift (closed form). Returns 0 for invalid params, which now
 * includes more than 50 parameters (the library's N_max): 1.x summed the first
 * 50 and returned that truncated value, 2.0 rejects and this sentinel-returning
 * wrapper reports the rejection as 0.
 */
double fos_z_shift(const double* params, int n_params);

/**
 * a2 from the volume-conservation constraint. Returns 0 for empty params and,
 * as for fos_z_shift, for more than 50 parameters — 1.x returned a truncated
 * sum there.
 */
double fos_a2(const double* params, int n_params);

#ifdef __cplusplus
}
#endif

#endif /* FOS_PARAMETERIZATION_H */
