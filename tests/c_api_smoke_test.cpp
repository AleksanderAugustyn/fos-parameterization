// Smoke test for the 2.0.0 C API. Exercises the SHARED library — the same
// binary the Python bindings load.
//
// Covers the whole surface shape: handle lifecycle (tables + caches), the
// cached computes with their explicit buffer sizes, the buffer-mismatch code
// (105), the flat tier-1 calls with their nullable trailing status, the raw
// rho(z) evaluator, and the static status-message table.
#include "fos_parameterization.h"

#include <cmath>
#include <cstdio>
#include <numbers>
#include <vector>

namespace {
int failures = 0;

void check(const bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}

constexpr int n_theta = 64;
constexpr int n_points = 501;

// Uniform theta nodes over [0, pi], endpoints pinned: under fast-math the
// product (n-1)*pi/(n-1) can land one ulp above pi, which the domain check
// would then reject.
std::vector<double> theta_grid(const int n) {
    std::vector<double> thetas(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        thetas[static_cast<std::size_t>(i)] =
                static_cast<double>(i) * std::numbers::pi / static_cast<double>(n - 1);
    }
    thetas[0] = 0.0;
    thetas[static_cast<std::size_t>(n - 1)] = std::numbers::pi;
    return thetas;
}

bool all_finite_positive(const std::vector<double>& v) {
    for (const double x : v) {
        if (!std::isfinite(x) || x <= 0.0) return false;
    }
    return true;
}

bool all_zero(const std::vector<double>& v) {
    for (const double x : v) {
        if (x != 0.0) return false;
    }
    return true;
}
}  // namespace

int main() {
    const std::vector<double> thetas = theta_grid(n_theta);
    const std::vector<double> params{1.5, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002};
    const int n_params = static_cast<int>(params.size());
    const std::vector<double> sphere{1.0};

    // --- Tables lifecycle ---------------------------------------------------
    void* tables = fos_param_tables_create(n_points, thetas.data(), n_theta);
    check(tables != nullptr, "tables_create(501, 64 thetas) is non-NULL");

    // Below the 100-node floor: NULL handle, invalid-grid code available via
    // the shared macro set.
    void* bad_tables = fos_param_tables_create(4, thetas.data(), n_theta);
    check(bad_tables == nullptr, "tables_create below the node floor is NULL");

    // --- Cache lifecycle ----------------------------------------------------
    void* cache = fos_param_cache_create_shared(tables, n_params);
    check(cache != nullptr, "cache_create_shared is non-NULL");

    void* own_cache = fos_param_cache_create(n_params, n_points, thetas.data(), n_theta);
    check(own_cache != nullptr, "cache_create (private tables) is non-NULL");

    // The cached tier caps n_params at 8.
    void* too_many = fos_param_cache_create_shared(tables, 9);
    check(too_many == nullptr, "cache_create_shared with 9 params is NULL");

    void* null_tables_cache = fos_param_cache_create_shared(nullptr, n_params);
    check(null_tables_cache == nullptr, "cache_create_shared(NULL tables) is NULL");

    // n_theta floor is 1, the same floor tables_create applies.
    void* no_theta_cache = fos_param_cache_create(n_params, n_points, thetas.data(), 0);
    check(no_theta_cache == nullptr, "cache_create with n_theta = 0 is NULL");
    void* neg_theta_cache = fos_param_cache_create(n_params, n_points, thetas.data(), -4);
    check(neg_theta_cache == nullptr, "cache_create with negative n_theta is NULL");
    void* neg_theta_tables = fos_param_tables_create(n_points, thetas.data(), -4);
    check(neg_theta_tables == nullptr, "tables_create with negative n_theta is NULL");

    // --- Cached computes ----------------------------------------------------
    std::vector<double> radii(n_theta, -1.0);
    int s = fos_param_cache_radius_grid(cache, params.data(), n_params,
                                        radii.data(), n_theta);
    check(s == FOS_VALID, "cache_radius_grid returns FOS_VALID");
    check(all_finite_positive(radii), "cache_radius_grid radii finite and positive");
    // Kept because the rejection probes below zero-fill `radii`.
    const std::vector<double> radii_ref = radii;

    // Wrong buffer size -> 105, outputs zero-filled.
    std::vector<double> short_radii(n_theta - 1, -1.0);
    s = fos_param_cache_radius_grid(cache, params.data(), n_params,
                                    short_radii.data(), n_theta - 1);
    check(s == FOS_ERROR_BUFFER_MISMATCH, "wrong n_radii -> FOS_ERROR_BUFFER_MISMATCH (105)");
    check(all_zero(short_radii), "buffer mismatch zero-fills the output");

    // A negative size is read as zero, so it mismatches the handle's extent
    // like any other wrong size — 105, not a crash and not an allocation.
    // (No large-stated-size case here: the contract is that the stated size IS
    // the buffer's real extent, so overstating it is undefined behavior by
    // construction. What is covered is that an outsized size cannot overflow
    // the stack, which the heap marshalling buffers guarantee structurally.)
    s = fos_param_cache_radius_grid(cache, params.data(), n_params,
                                    short_radii.data(), -5);
    check(s == FOS_ERROR_BUFFER_MISMATCH, "negative n_radii -> FOS_ERROR_BUFFER_MISMATCH (105)");

    double neg_shift = -1.0;
    std::vector<double> gz_probe(1, -1.0), rho_probe(1, -1.0), drho_probe(1, -1.0);
    s = fos_param_cache_rho_z_grid(cache, params.data(), n_params, gz_probe.data(),
                                   rho_probe.data(), drho_probe.data(), -1, &neg_shift);
    check(s == FOS_ERROR_BUFFER_MISMATCH, "negative n_z -> FOS_ERROR_BUFFER_MISMATCH (105)");

    // Wrong parameter count -> 4 (the cached tier is exact-length).
    s = fos_param_cache_radius_grid(cache, params.data(), n_params - 1,
                                    radii.data(), n_theta);
    check(s == FOS_ERROR_WRONG_PARAM_COUNT, "short params -> FOS_ERROR_WRONG_PARAM_COUNT (4)");

    // Precedence 4 > 105: both wrong at once still reports the parameter count.
    s = fos_param_cache_radius_grid(cache, params.data(), n_params - 1,
                                    short_radii.data(), n_theta - 1);
    check(s == FOS_ERROR_WRONG_PARAM_COUNT, "wrong params + wrong n_radii -> 4");

    // NULL handle -> 2, never a dereference.
    s = fos_param_cache_radius_grid(nullptr, params.data(), n_params,
                                    radii.data(), n_theta);
    check(s == FOS_ERROR_CACHE_NOT_INITIALIZED, "NULL cache -> FOS_ERROR_CACHE_NOT_INITIALIZED (2)");

    std::vector<double> radii_d(n_theta, -1.0), dr(n_theta, -1.0);
    s = fos_param_cache_radius_and_derivative(cache, params.data(), n_params,
                                              radii_d.data(), dr.data(), n_theta);
    check(s == FOS_VALID, "cache_radius_and_derivative returns FOS_VALID");
    check(all_finite_positive(radii_d), "cache_radius_and_derivative radii positive");

    const std::vector<double> extra_thetas{0.4, 1.5707963267948966, 2.7};
    std::vector<double> extra_r(3, -1.0), extra_dr(3, -1.0);
    s = fos_param_cache_radius_and_derivative_at_thetas(
            cache, params.data(), n_params, extra_thetas.data(), 3,
            extra_r.data(), extra_dr.data());
    check(s == FOS_VALID, "cache_radius_and_derivative_at_thetas returns FOS_VALID");
    check(all_finite_positive(extra_r), "at_thetas radii positive");

    // A negative size is not an extent on this form either: 105, nothing
    // written. (Clamping it to zero would make the stated and the expected
    // extent agree and return FOS_VALID having computed nothing.)
    s = fos_param_cache_radius_and_derivative_at_thetas(
            cache, params.data(), n_params, extra_thetas.data(), -3,
            extra_r.data(), extra_dr.data());
    check(s == FOS_ERROR_BUFFER_MISMATCH, "negative n_thetas -> FOS_ERROR_BUFFER_MISMATCH (105)");

    // A wrong buffer size is the CALLER's error and outranks the SHAPE's.
    // Probe the beak vector at the right size first, so the 105 below is known
    // to be reported instead of a live 103 and not merely instead of success.
    std::vector<double> beak_params(static_cast<std::size_t>(n_params), 0.0);
    beak_params[0] = 2.0;
    beak_params[2] = 0.7495;  // the frozen 1.x beak probe vector
    s = fos_param_cache_radius_grid(cache, beak_params.data(), n_params,
                                    radii.data(), n_theta);
    check(s == FOS_ERROR_BEAK_SINGULARITY, "beak vector at the right size -> 103");
    s = fos_param_cache_radius_grid(cache, beak_params.data(), n_params,
                                    short_radii.data(), n_theta - 1);
    check(s == FOS_ERROR_BUFFER_MISMATCH,
          "beak vector with wrong n_radii -> 105 (size outranks shape)");

    double z_shift = -1.0, r_north = -1.0, r_south = -1.0;
    s = fos_param_cache_shape(cache, params.data(), n_params, &z_shift, &r_north, &r_south);
    check(s == FOS_VALID, "cache_shape returns FOS_VALID");
    check(r_north > 0.0 && r_south > 0.0, "cache_shape pole radii positive");

    std::vector<double> gz(n_points, -1.0), grho(n_points, -1.0), gdrho(n_points, -1.0);
    double grid_shift = -1.0;
    s = fos_param_cache_rho_z_grid(cache, params.data(), n_params, gz.data(),
                                   grho.data(), gdrho.data(), n_points, &grid_shift);
    check(s == FOS_VALID, "cache_rho_z_grid returns FOS_VALID");
    check(grho[0] == 0.0 && grho[n_points - 1] == 0.0, "rho = 0 at the poles");

    s = fos_param_cache_rho_z_grid(cache, params.data(), n_params, gz.data(),
                                   grho.data(), gdrho.data(), n_points - 1, &grid_shift);
    check(s == FOS_ERROR_BUFFER_MISMATCH, "wrong n_z -> FOS_ERROR_BUFFER_MISMATCH (105)");

    double z_neck = -1.0, rho_neck = -1.0;
    int found = -1;
    s = fos_param_cache_neck(cache, params.data(), n_params, &z_neck, &rho_neck, &found);
    check(s == FOS_VALID, "cache_neck returns FOS_VALID");
    check(found == 0 || found == 1, "cache_neck found is a 0/1 flag");

    double z_shift_total = -1.0, g_opt = -1.0;
    s = fos_param_cache_star_convexity_optimum(cache, params.data(), n_params,
                                               &z_shift_total, &g_opt);
    check(s == FOS_VALID, "cache_star_convexity_optimum returns FOS_VALID");

    // A private-tables cache must agree with the shared-tables one.
    std::vector<double> own_radii(n_theta, -1.0);
    s = fos_param_cache_radius_grid(own_cache, params.data(), n_params,
                                    own_radii.data(), n_theta);
    check(s == FOS_VALID, "private-tables cache returns FOS_VALID");
    bool same = true;
    for (int i = 0; i < n_theta; ++i) {
        same = same && own_radii[static_cast<std::size_t>(i)] == radii_ref[static_cast<std::size_t>(i)];
    }
    check(same, "private and shared tables give identical radii");

    // --- Flat tier-1 calls --------------------------------------------------
    int st = -1;
    std::vector<double> flat_radii(n_theta, -1.0);
    fos_param_radius_grid(params.data(), n_params, thetas.data(), n_theta,
                          n_points, flat_radii.data(), &st);
    check(st == FOS_VALID, "flat radius_grid returns FOS_VALID");
    check(all_finite_positive(flat_radii), "flat radius_grid radii positive");

    std::vector<double> flat_dr(n_theta, -1.0);
    fos_param_radius_and_derivative(params.data(), n_params, thetas.data(), n_theta,
                                    n_points, flat_radii.data(), flat_dr.data(), &st);
    check(st == FOS_VALID, "flat radius_and_derivative returns FOS_VALID");

    z_shift = r_north = r_south = -1.0;
    fos_param_shape(params.data(), n_params, n_points, &z_shift, &r_north, &r_south, &st);
    check(st == FOS_VALID, "flat shape returns FOS_VALID");
    check(r_north > 0.0 && r_south > 0.0, "flat shape pole radii positive");

    fos_param_rho_z_grid(params.data(), n_params, n_points, gz.data(), grho.data(),
                         gdrho.data(), &grid_shift, &st);
    check(st == FOS_VALID, "flat rho_z_grid returns FOS_VALID");

    fos_param_neck(params.data(), n_params, n_points, &z_neck, &rho_neck, &found, &st);
    check(st == FOS_VALID, "flat neck returns FOS_VALID");

    fos_param_star_convexity_optimum(params.data(), n_params, n_points,
                                     &z_shift_total, &g_opt, &st);
    check(st == FOS_VALID, "flat star_convexity_optimum returns FOS_VALID");

    double zs = -1.0;
    fos_param_z_shift(params.data(), n_params, &zs, &st);
    check(st == FOS_VALID, "flat z_shift returns FOS_VALID");

    double a2 = -1.0;
    // a4-only shape: a2 = a4/3.
    const std::vector<double> a4_only{2.0, 0.0, 0.5};
    fos_param_a2(a4_only.data(), 3, &a2, &st);
    check(st == FOS_VALID && std::fabs(a2 - 0.5 / 3.0) < 1e-15, "flat a2 = a4/3");

    // A negative size on a flat call is a bad C argument, not a handle
    // mismatch: FOS_ERROR_INVALID_INIT (5).
    fos_param_radius_grid(params.data(), n_params, thetas.data(), -3, n_points,
                          flat_radii.data(), &st);
    check(st == FOS_ERROR_INVALID_INIT, "flat radius_grid, negative n_thetas -> 5");
    fos_param_z_shift(params.data(), -2, &zs, &st);
    check(st == FOS_ERROR_INVALID_INIT, "flat z_shift, negative n_params -> 5");

    // Rejections still travel through the flat tier.
    const std::vector<double> bad_c{0.0, 0.0, 0.0};
    fos_param_z_shift(bad_c.data(), 3, &zs, &st);
    check(st == FOS_ERROR_INVALID_C, "flat z_shift with c = 0 -> FOS_ERROR_INVALID_C (102)");
    check(zs == 0.0, "rejected flat z_shift zero-fills its output");

    // A NULL status pointer is an absent optional dummy on the Fortran side.
    fos_param_z_shift(params.data(), n_params, &zs, nullptr);
    fos_param_shape(params.data(), n_params, n_points, &z_shift, &r_north, &r_south, nullptr);
    check(true, "NULL status pointer does not crash");

    // --- Raw evaluator ------------------------------------------------------
    double rho = -1.0, drho_dz = -1.0;
    fos_param_rho_at_z(sphere.data(), 1, 0.0, 0.0, &rho, &drho_dz);
    check(std::fabs(rho - 1.0) < 1e-14, "rho_at_z(sphere, z = 0) = 1");

    // --- Status messages ----------------------------------------------------
    const char* msg = fos_param_status_message(FOS_ERROR_INVALID_C);
    check(msg != nullptr && msg[0] != '\0', "status_message(102) is non-empty");
    const char* ok_msg = fos_param_status_message(FOS_VALID);
    check(ok_msg != nullptr && ok_msg[0] != '\0', "status_message(0) is non-empty");
    const char* unknown = fos_param_status_message(-12345);
    check(unknown != nullptr && unknown[0] != '\0', "status_message(unknown) falls back");

    // --- Teardown: caches first, tables last --------------------------------
    fos_param_cache_destroy(cache);
    fos_param_cache_destroy(own_cache);
    fos_param_cache_destroy(nullptr);  // NULL-safe
    fos_param_tables_destroy(tables);
    fos_param_tables_destroy(nullptr);  // NULL-safe

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
