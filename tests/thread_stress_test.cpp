// Contract family 3: concurrency. ONE read-only cache, shared by every thread.
//
// A cache is immutable after creation, so any number of threads may compute on
// one handle at once. The main thread first walks every shape of every thread
// serially through the shared handle and records all eight cached outputs.
// Then 8 threads repeat exactly their own walk concurrently on the SAME handle.
// Concurrency is the only difference, and it may not move a bit: every status
// must match and every output must be identical byte for byte.
//
// A data race on the shared cache, or any state leaking from one call to the
// next, shows up as a status mismatch or a differing byte.
#include "fos_parameterization.h"

#include <cstdio>
#include <cstring>
#include <numbers>
#include <thread>
#include <vector>

namespace {
constexpr int n_theta = 64;
constexpr int n_other = 9;
constexpr int n_points = 501;
constexpr int n_params = 7;
constexpr int n_threads = 8;
constexpr int n_steps = 100;

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

// Deterministic per-(thread, step) shape. Amplitudes stay small enough that
// every shape resolves; the walk differs per thread, and the vector length
// cycles through 3, 5 and 7 so short vectors are exercised concurrently too.
std::vector<double> shape_for(const int thread_id, const int step) {
    const double t = static_cast<double>(thread_id + 1) * 0.01;
    const double u = static_cast<double>(step % 41) * 0.005;
    std::vector<double> p{1.30 + u + t, 0.05 * t, 0.10 - 0.5 * u, 0.02, 0.01, 0.005, 0.002};
    p.resize(static_cast<std::size_t>(3 + 2 * (step % 3)));
    return p;
}

// One step's outputs: all eight cached computes.
struct StepResult {
    int status[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    std::vector<double> radii, radii_d, dr, at_radii, at_dr;
    std::vector<double> z, rho, drho, raw_z, raw_rho, raw_drho;
    double scalars[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    int found = -1;
};

using Walk = std::vector<StepResult>;

// The one code path both the serial reference and the threads use, so any
// difference is concurrency — never a different computation.
void run_walk(const fos_param_cache_t* cache, const std::vector<double>& other_thetas,
              const int thread_id, Walk& out) {
    out.assign(static_cast<std::size_t>(n_steps), StepResult{});
    for (int step = 0; step < n_steps; ++step) {
        StepResult& r = out[static_cast<std::size_t>(step)];
        const std::vector<double> params = shape_for(thread_id, step);
        const int n = static_cast<int>(params.size());

        r.radii.assign(n_theta, 0.0);
        r.radii_d.assign(n_theta, 0.0);
        r.dr.assign(n_theta, 0.0);
        r.at_radii.assign(n_other, 0.0);
        r.at_dr.assign(n_other, 0.0);
        r.z.assign(n_points, 0.0);
        r.rho.assign(n_points, 0.0);
        r.drho.assign(n_points, 0.0);
        r.raw_z.assign(n_points, 0.0);
        r.raw_rho.assign(n_points, 0.0);
        r.raw_drho.assign(n_points, 0.0);

        r.status[0] = fos_param_cache_radius_grid(cache, params.data(), n,
                                                  r.radii.data(), n_theta);
        r.status[1] = fos_param_cache_radius_and_derivative(
                cache, params.data(), n, r.radii_d.data(), r.dr.data(), n_theta);
        r.status[2] = fos_param_cache_radius_and_derivative_at_thetas(
                cache, params.data(), n, other_thetas.data(), n_other,
                r.at_radii.data(), r.at_dr.data());
        r.status[3] = fos_param_cache_shape(cache, params.data(), n, &r.scalars[0],
                                            &r.scalars[1], &r.scalars[2]);
        r.status[4] = fos_param_cache_rho_z_grid(cache, params.data(), n, r.z.data(),
                                                 r.rho.data(), r.drho.data(), n_points,
                                                 &r.scalars[3]);
        r.status[5] = fos_param_cache_rho_z_grid_unchecked(
                cache, params.data(), n, r.raw_z.data(), r.raw_rho.data(),
                r.raw_drho.data(), n_points, &r.scalars[4]);
        r.status[6] = fos_param_cache_neck(cache, params.data(), n, &r.scalars[5],
                                           &r.scalars[6], &r.found);
        r.status[7] = fos_param_cache_star_convexity_optimum(
                cache, params.data(), n, &r.scalars[7], &r.scalars[8]);
    }
}

bool same_bytes(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size()
           && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

bool same_step(const StepResult& a, const StepResult& b) {
    return std::memcmp(a.status, b.status, sizeof(a.status)) == 0
           && std::memcmp(a.scalars, b.scalars, sizeof(a.scalars)) == 0
           && a.found == b.found
           && same_bytes(a.radii, b.radii) && same_bytes(a.radii_d, b.radii_d)
           && same_bytes(a.dr, b.dr) && same_bytes(a.at_radii, b.at_radii)
           && same_bytes(a.at_dr, b.at_dr) && same_bytes(a.z, b.z)
           && same_bytes(a.rho, b.rho) && same_bytes(a.drho, b.drho)
           && same_bytes(a.raw_z, b.raw_z) && same_bytes(a.raw_rho, b.raw_rho)
           && same_bytes(a.raw_drho, b.raw_drho);
}

bool all_valid(const StepResult& r) {
    for (const int s : r.status) {
        if (s != FOS_VALID) return false;
    }
    return true;
}
}  // namespace

int main() {
    const std::vector<double> thetas = theta_grid(n_theta);
    const std::vector<double> other_thetas = theta_grid(n_other);

    int status = -1;
    fos_param_cache_t* cache =
            fos_param_cache_create(n_params, n_points, thetas.data(), n_theta, &status);
    if (cache == nullptr) {
        std::printf("FAIL: cache_create returned NULL (status %d)\n", status);
        return 1;
    }

    // Serial reference first, on the same handle the threads will share.
    std::vector<Walk> serial(static_cast<std::size_t>(n_threads));
    for (int t = 0; t < n_threads; ++t) {
        run_walk(cache, other_thetas, t, serial[static_cast<std::size_t>(t)]);
    }

    std::vector<Walk> threaded(static_cast<std::size_t>(n_threads));
    {
        std::vector<std::jthread> pool;
        pool.reserve(static_cast<std::size_t>(n_threads));
        for (int t = 0; t < n_threads; ++t) {
            pool.emplace_back([cache, &other_thetas, &threaded, t] {
                run_walk(cache, other_thetas, t, threaded[static_cast<std::size_t>(t)]);
            });
        }
    }  // all joined here

    int failures = 0;
    int valid_steps = 0;
    for (int t = 0; t < n_threads; ++t) {
        const Walk& mt = threaded[static_cast<std::size_t>(t)];
        const Walk& ref = serial[static_cast<std::size_t>(t)];
        if (mt.size() != ref.size()) {
            std::printf("FAIL: thread %d recorded %zu steps, reference %zu\n", t,
                        mt.size(), ref.size());
            ++failures;
            continue;
        }
        for (int step = 0; step < n_steps; ++step) {
            const StepResult& a = mt[static_cast<std::size_t>(step)];
            const StepResult& b = ref[static_cast<std::size_t>(step)];
            if (!same_step(a, b)) {
                std::printf("FAIL: thread %d step %d differs from the serial reference\n",
                            t, step);
                ++failures;
            }
            if (all_valid(b)) ++valid_steps;
        }
    }

    // A rejected shape zero-fills every buffer, and two zero buffers always
    // match — the comparison must not be able to pass vacuously.
    if (valid_steps != n_threads * n_steps) {
        std::printf("FAIL: only %d of %d steps were fully VALID\n", valid_steps,
                    n_threads * n_steps);
        ++failures;
    }

    fos_param_cache_destroy(cache);

    std::printf("thread_stress_test: %d thread(s) x %d steps on one shared cache, "
                "%d failure(s)\n", n_threads, n_steps, failures);
    return failures == 0 ? 0 : 1;
}
