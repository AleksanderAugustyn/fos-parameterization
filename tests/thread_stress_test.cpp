// Thread model of the 2.0.0 C API, positive test: SHARE TABLES, CONFINE CACHES.
//
// One tables handle is built on the main thread and read concurrently by 8
// worker threads. Every thread owns a cache created with
// fos_param_cache_create_shared over that one tables handle; no cache handle
// ever crosses a thread boundary. Each thread walks 200 private shapes through
// its incremental cache, recording radius_grid, shape and rho_z_grid for every
// step.
//
// After the join, the SAME walks are replayed serially on the main thread
// through a SINGLE cache — same shared tables, same call order, same
// arithmetic. Concurrency and cache history are the only differences, and
// neither may move a bit: the required maximum absolute difference is exactly
// 0.0, and every status must match.
//
// A data race on the shared tables, or cache state leaking between threads,
// shows up as a status mismatch or a nonzero difference. Destruction order is
// the documented one: every cache first, the tables last.
#include "fos_parameterization.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <numbers>
#include <thread>
#include <vector>

namespace {

constexpr int n_theta = 64;
constexpr int n_points = 501;
constexpr int n_params = 7;
constexpr int n_threads = 8;
constexpr int n_steps = 200;

// Uniform theta nodes over [0, pi], endpoints pinned (see the smoke test).
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
// every shape resolves; the walk differs per thread so the caches never see
// the same parameter sequence.
std::vector<double> shape_for(const int thread_id, const int step) {
    const double t = static_cast<double>(thread_id + 1) * 0.01;
    const double u = static_cast<double>(step % 41) * 0.005;
    return {1.30 + u + t, 0.05 * t, 0.10 - 0.5 * u, 0.02, 0.01, 0.005, 0.002};
}

// One step's outputs, recorded for the bit-comparison.
struct StepResult {
    int radius_status = -1;
    int shape_status = -1;
    int rho_status = -1;
    std::vector<double> radii;
    double z_shift = 0.0;
    double r_north = 0.0;
    double r_south = 0.0;
    std::vector<double> z;
    std::vector<double> rho;
    std::vector<double> drho_dz;
    double grid_shift = 0.0;
};

using Walk = std::vector<StepResult>;

// The one code path both the threaded run and the serial replay use, so any
// difference is concurrency or cache history — never a different computation.
void run_walk(void* cache, const int thread_id, Walk& out) {
    out.assign(static_cast<std::size_t>(n_steps), StepResult{});
    for (int step = 0; step < n_steps; ++step) {
        StepResult& r = out[static_cast<std::size_t>(step)];
        const std::vector<double> params = shape_for(thread_id, step);

        r.radii.assign(static_cast<std::size_t>(n_theta), 0.0);
        r.radius_status = fos_param_cache_radius_grid(cache, params.data(), n_params,
                                                      r.radii.data(), n_theta);

        r.shape_status = fos_param_cache_shape(cache, params.data(), n_params,
                                               &r.z_shift, &r.r_north, &r.r_south);

        r.z.assign(static_cast<std::size_t>(n_points), 0.0);
        r.rho.assign(static_cast<std::size_t>(n_points), 0.0);
        r.drho_dz.assign(static_cast<std::size_t>(n_points), 0.0);
        r.rho_status = fos_param_cache_rho_z_grid(cache, params.data(), n_params,
                                                  r.z.data(), r.rho.data(),
                                                  r.drho_dz.data(), n_points,
                                                  &r.grid_shift);
    }
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double worst = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        worst = std::max(worst, std::fabs(a[i] - b[i]));
    }
    return worst;
}

double scalar_diff(const double a, const double b) { return std::fabs(a - b); }

}  // namespace

int main() {
    const std::vector<double> thetas = theta_grid(n_theta);

    void* tables = fos_param_tables_create(n_points, thetas.data(), n_theta);
    if (tables == nullptr) {
        std::printf("FAIL: tables_create returned NULL\n");
        return 1;
    }

    std::vector<Walk> threaded(static_cast<std::size_t>(n_threads));
    // char, not bool: std::vector<bool> packs bits, so concurrent writes to
    // distinct "elements" would be a genuine data race.
    std::vector<char> created(static_cast<std::size_t>(n_threads), 0);

    const auto worker = [&tables, &threaded, &created](const int id) {
        void* cache = fos_param_cache_create_shared(tables, n_params);
        if (cache == nullptr) return;
        created[static_cast<std::size_t>(id)] = 1;
        run_walk(cache, id, threaded[static_cast<std::size_t>(id)]);
        fos_param_cache_destroy(cache);
    };

    {
        std::vector<std::jthread> pool;
        pool.reserve(static_cast<std::size_t>(n_threads));
        for (int t = 0; t < n_threads; ++t) pool.emplace_back(worker, t);
    }  // all joined here

    int failures = 0;
    for (int t = 0; t < n_threads; ++t) {
        if (!created[static_cast<std::size_t>(t)]) {
            std::printf("FAIL: thread %d could not create its cache\n", t);
            ++failures;
        }
    }

    // Serial replay: ONE cache walks every thread's parameter sequence in turn.
    void* serial_cache = fos_param_cache_create_shared(tables, n_params);
    if (serial_cache == nullptr) {
        std::printf("FAIL: serial cache_create_shared returned NULL\n");
        fos_param_tables_destroy(tables);
        return 1;
    }

    double worst = 0.0;
    int valid_computes = 0;
    for (int t = 0; t < n_threads; ++t) {
        Walk serial;
        run_walk(serial_cache, t, serial);

        const Walk& mt = threaded[static_cast<std::size_t>(t)];
        if (mt.size() != serial.size()) {
            std::printf("FAIL: thread %d recorded %zu steps, replay %zu\n", t,
                        mt.size(), serial.size());
            ++failures;
            continue;
        }

        for (int step = 0; step < n_steps; ++step) {
            const StepResult& a = mt[static_cast<std::size_t>(step)];
            const StepResult& b = serial[static_cast<std::size_t>(step)];

            if (a.radius_status != b.radius_status || a.shape_status != b.shape_status
                || a.rho_status != b.rho_status) {
                std::printf("FAIL: thread %d step %d status (%d,%d,%d) vs replay "
                            "(%d,%d,%d)\n", t, step, a.radius_status, a.shape_status,
                            a.rho_status, b.radius_status, b.shape_status, b.rho_status);
                ++failures;
                continue;
            }
            if (a.radius_status == FOS_VALID && a.shape_status == FOS_VALID
                && a.rho_status == FOS_VALID) {
                ++valid_computes;
            }

            worst = std::max(worst, max_abs_diff(a.radii, b.radii));
            worst = std::max(worst, max_abs_diff(a.z, b.z));
            worst = std::max(worst, max_abs_diff(a.rho, b.rho));
            worst = std::max(worst, max_abs_diff(a.drho_dz, b.drho_dz));
            worst = std::max(worst, scalar_diff(a.z_shift, b.z_shift));
            worst = std::max(worst, scalar_diff(a.r_north, b.r_north));
            worst = std::max(worst, scalar_diff(a.r_south, b.r_south));
            worst = std::max(worst, scalar_diff(a.grid_shift, b.grid_shift));
        }
    }

    if (worst != 0.0) {
        std::printf("FAIL: max abs diff %.17g, required exactly 0\n", worst);
        ++failures;
    }
    // A rejected shape zero-fills every buffer, and two zero buffers always
    // match — the comparison must not be able to pass vacuously.
    if (valid_computes != n_threads * n_steps) {
        std::printf("FAIL: only %d of %d steps were fully VALID\n", valid_computes,
                    n_threads * n_steps);
        ++failures;
    }

    // Documented teardown order: every cache first, the tables last.
    fos_param_cache_destroy(serial_cache);
    fos_param_tables_destroy(tables);

    std::printf("thread_stress_test: %d thread(s) x %d steps, max abs diff %.17g, "
                "%d failure(s)\n", n_threads, n_steps, worst, failures);
    return failures == 0 ? 0 : 1;
}
