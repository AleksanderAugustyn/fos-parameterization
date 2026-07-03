/**
 * @file fos_parameterization.h
 * @brief C API for the Fortran Fourier-over-Spheroid (FoS) shape library (v1.0.1).
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
 * codes 0-4 mirror the Fortran FOS_* parameters exactly.
 *
 * Thread safety: the library is stateless; all functions are safe to call
 * concurrently.
 */

#ifndef FOS_PARAMETERIZATION_H
#define FOS_PARAMETERIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

/* --- Status codes (0-4 mirror Fortran FOS_* parameters) --- */
#define FOS_VALID                    0
#define FOS_ERROR_RHO_NEGATIVE       1
#define FOS_ERROR_NOT_STAR_CONVEX    2
#define FOS_ERROR_INVALID_C          3
#define FOS_ERROR_BEAK_SINGULARITY   4
#define FOS_ERROR_INVALID_ARGUMENTS  5

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

/** Intrinsic COM z-shift (closed form). Returns 0 for invalid params. */
double fos_z_shift(const double* params, int n_params);

/** a2 from the volume-conservation constraint. Returns 0 for empty params. */
double fos_a2(const double* params, int n_params);

#ifdef __cplusplus
}
#endif

#endif /* FOS_PARAMETERIZATION_H */
