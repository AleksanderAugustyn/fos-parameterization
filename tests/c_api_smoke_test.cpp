// Smoke test for the C API. Exercises the SHARED library — the same binary
// the Python bindings load.
#include "fos_parameterization.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {
int failures = 0;

void check(bool ok, const char* label) {
    if (!ok) {
        std::printf("FAIL: %s\n", label);
        ++failures;
    }
}
}  // namespace

int main() {
    std::array<char, 256> buf{};
    const int nbuf = static_cast<int>(buf.size());

    // --- fos_compute_radius_grid: sphere ---
    std::vector<double> sphere{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> radii(181, -1.0);
    double z_shift = -1.0;
    int s = fos_compute_radius_grid(sphere.data(), 7, 181, radii.data(), &z_shift, nbuf, buf.data());
    check(s == FOS_VALID, "C: sphere returns FOS_VALID");
    bool all_one = true;
    for (double r : radii) all_one = all_one && (std::fabs(r - 1.0) < 5e-12);
    check(all_one, "C: sphere R(theta) = 1");
    check(std::fabs(z_shift) < 1e-14, "C: sphere z_shift = 0");
    check(buf[0] == '\0', "C: success message is empty");

    // --- error paths ---
    std::vector<double> bad_c{0.0, 0.0, 0.0};
    std::vector<double> radii_fail(41, -7.0);
    z_shift = -7.0;
    s = fos_compute_radius_grid(bad_c.data(), 3, 41, radii_fail.data(), &z_shift, nbuf, buf.data());
    check(s == FOS_ERROR_INVALID_C, "C: c = 0 -> FOS_ERROR_INVALID_C");
    check(buf[0] != '\0', "C: failure message is non-empty");
    bool zeros = (z_shift == 0.0);
    for (double r : radii_fail) zeros = zeros && (r == 0.0);
    check(zeros, "C: failed call zero-fills radii and z_shift");

    s = fos_compute_radius_grid(sphere.data(), 7, 1, radii_fail.data(), &z_shift, nbuf, buf.data());
    check(s == FOS_ERROR_INVALID_ARGUMENTS, "C: n_grid = 1 -> FOS_ERROR_INVALID_ARGUMENTS");

    // Message truncation: tiny buffer must still be null-terminated within bounds.
    std::array<char, 8> tiny{};
    tiny.fill('X');
    s = fos_compute_radius_grid(bad_c.data(), 3, 41, radii_fail.data(), &z_shift,
                                static_cast<int>(tiny.size()), tiny.data());
    check(std::memchr(tiny.data(), '\0', tiny.size()) != nullptr,
          "C: truncated message is null-terminated within the buffer");

    // Beak error code surfaces through the C API.
    std::vector<double> beak{2.0, 0.0, 0.7497, 0.0, 0.0, 0.0, 0.0};
    s = fos_compute_radius_grid(beak.data(), 7, 41, radii_fail.data(), &z_shift, nbuf, buf.data());
    check(s == FOS_ERROR_BEAK_SINGULARITY, "C: near-scission beak -> code 4");

    // --- fos_compute_rho_profile: sphere, odd n_z so z = 0 is on-grid ---
    const int nz = 201;
    std::vector<double> z(nz, -1.0), rho(nz, -1.0), drho(nz, -1.0);
    s = fos_compute_rho_profile(sphere.data(), 7, nz, z.data(), rho.data(), drho.data(), nbuf, buf.data());
    check(s == FOS_VALID, "C: sphere profile FOS_VALID");
    check(std::fabs(z[nz / 2]) < 1e-14, "C: sphere profile midpoint z = 0");
    check(std::fabs(rho[nz / 2] - 1.0) < 1e-12, "C: sphere profile rho(0) = 1");
    check(rho[0] == 0.0 && rho[nz - 1] == 0.0, "C: sphere profile rho = 0 at poles");

    s = fos_compute_rho_profile(sphere.data(), 7, 1, z.data(), rho.data(), drho.data(), nbuf, buf.data());
    check(s == FOS_ERROR_INVALID_ARGUMENTS, "C: n_z = 1 -> FOS_ERROR_INVALID_ARGUMENTS");

    // --- fos_compute_neck: symmetric family closed form ---
    // c=2, a4=0.5: rho_neck = sqrt((1 - 4*0.5/3)/2) = sqrt(1/6), z_neck = 0.
    std::vector<double> necked{2.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0};
    double z_neck = -1.0, rho_neck = -1.0;
    int found = -1;
    s = fos_compute_neck(necked.data(), 7, &z_neck, &rho_neck, &found);
    check(s == FOS_VALID, "C: neck call FOS_VALID");
    check(found == 1, "C: neck found");
    check(std::fabs(rho_neck - std::sqrt(1.0 / 6.0)) < 1e-9, "C: rho_neck analytic");
    check(std::fabs(z_neck) < 1e-9, "C: z_neck = 0");

    s = fos_compute_neck(sphere.data(), 7, &z_neck, &rho_neck, &found);
    check(s == FOS_VALID && found == 0, "C: sphere has no neck, still FOS_VALID");

    s = fos_compute_neck(bad_c.data(), 3, &z_neck, &rho_neck, &found);
    check(s == FOS_ERROR_INVALID_C, "C: neck with c = 0 -> FOS_ERROR_INVALID_C");

    // --- fos_z_shift / fos_a2 ---
    check(fos_z_shift(necked.data(), 7) == 0.0, "C: symmetric shape z_shift = 0");
    check(std::fabs(fos_a2(necked.data(), 7) - 0.5 / 3.0) < 1e-15, "C: a2 = a4/3 for a4-only");
    check(fos_a2(sphere.data(), 7) == 0.0, "C: sphere a2 = 0");

    std::printf("c_api_smoke_test: %d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
