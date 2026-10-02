# fos-parameterization

Fourier-over-Spheroid (FoS) nuclear shape parameterization: a Fortran 2018
library with a C API and Python bindings. It turns a parameter vector
`[c, a3, a4, ...]` into the nuclear surface, either as the cylindrical profile
rho(z) or as R(theta) about a star-convex origin.

The library follows the two-tier shape parameterization contract shared with
[beta-parameterization](https://github.com/AleksanderAugustyn/beta-parameterization):
one-shot functions, plus a read-only cache that is built once and shared across
threads.

## The shape

In reduced units (R0 = 1) the surface of an axially symmetric nucleus is

    rho^2(z) = f(u) / c,        u = (z - z_shift) / c,     -1 <= u <= 1

    f(u) = 1 - u^2 - sum_{k>=1} [ a_2k cos((2k - 1) pi u / 2) + a_2k+1 sin(k pi u) ]

- `c` is the elongation: the half-length of the shape in units of R0.
- `a3, a4, a5, ...` are the deformation coefficients. Odd coefficients break
  the left-right symmetry; `a4` controls the neck.
- `a2` is not a parameter. Volume conservation fixes it:
  `a2 = a4/3 - a6/5 + a8/7 - ...`.
- The parameter vector is `params = [c, a3, a4, a5, ...]`, 1 to 50 entries.
  Missing trailing parameters are zero.

`z_shift` places the centre of mass at the origin. The R(theta) outputs need a
star-convex origin as well, so they keep the centre of mass when it is
well-conditioned and otherwise move to the star-convexity optimum; they report
the total shift they used.

The library rejects only what the mathematics requires: `c > 0`, `rho > 0` in
the interior, and the two numerical conditions of the R(theta) conversion (the
beak threshold and the star-convexity margin). Filtering shapes that are valid
but physically meaningless is the consumer's job.

## Two tiers

1. **One-shot.** One call computes a shape. All workspace is internal and
   discarded on return. For scripts and one-off plots.
2. **Read-only cache.** Build a cache once, then request any number of shapes
   against it. The cache holds only what depends on the resolution (the u grid,
   the Fourier basis, the primary thetas). Nothing derived from a shape
   parameter is stored, so no call influences a later one, and any number of
   threads may compute on one cache at once. For plotting loops and hot paths.

The two tiers return the same status and bitwise-identical outputs. A short
parameter vector gives the same bits as its zero-padded form.

`max_params`, chosen when the cache is built, is the longest vector that cache
accepts; each call may pass any length from 1 to `max_params`.

| Output | Cached | One-shot |
|---|---|---|
| R(theta) at the primary thetas | `cache_radius_grid_s` | `compute_radius_grid_standalone_s` |
| R(theta) and dR/dtheta | `cache_radius_and_derivative_s` | `compute_radius_and_derivative_standalone_s` |
| R(theta) and dR/dtheta at caller thetas | `cache_radius_and_derivative_at_thetas_s` | (the one-shot forms take caller thetas) |
| total z-shift and pole radii | `cache_shape_s` | `compute_shape_standalone_s` |
| rho(z) profile | `cache_rho_z_grid_s` | `compute_rho_z_grid_standalone_s` |
| rho(z) profile, separated shapes included | `cache_rho_z_grid_unchecked_s` | `compute_rho_z_grid_unchecked_standalone_s` |
| neck position and radius | `cache_neck_s` | `compute_neck_standalone_s` |
| star-convexity optimum | `cache_star_convexity_optimum_s` | `compute_star_convexity_optimum_standalone_s` |

The C functions are `fos_param_cache_<output>` and `fos_param_<output>`; the
Python names are the methods of `Cache` and the module-level functions.

### Fortran

```fortran
program fos_example
    use precision_utilities_mod, only: ik, rk
    use mathematical_and_physical_constants_mod, only: PI_C
    use fos_parameterization_mod, only: cache_t, cache_init_s, cache_free_s, &
            cache_radius_and_derivative_s, compute_radius_grid_standalone_s, &
            status_message, SHAPE_VALID
    implicit none

    integer(kind = ik), parameter :: N_THETA = 64_ik, N_POINTS = 501_ik
    real(kind = rk) :: thetas(N_THETA), radii(N_THETA), dr_dtheta(N_THETA)
    real(kind = rk) :: params(3)
    type(cache_t) :: cache
    integer(kind = ik) :: i, status

    do i = 1_ik, N_THETA
        thetas(i) = real(i, rk) * PI_C / real(N_THETA + 1_ik, rk)
    end do
    params = [1.5_rk, 0.1_rk, 0.05_rk]      ! c, a3, a4

    ! Tier 1: one-shot
    call compute_radius_grid_standalone_s(params, thetas, N_POINTS, radii, status)
    if (status /= SHAPE_VALID) print *, trim(status_message(status))

    ! Tier 2: build the cache once, then compute; it may be shared across threads
    call cache_init_s(cache, 8_ik, N_POINTS, thetas, status)
    call cache_radius_and_derivative_s(cache, params, radii, dr_dtheta, status)
    if (status /= SHAPE_VALID) print *, trim(status_message(status))
    call cache_free_s(cache)
end program fos_example
```

### C

```c
#include "fos_parameterization.h"
#include <stdio.h>

int main(void) {
    enum { N_THETA = 64, N_POINTS = 501 };
    double thetas[N_THETA], radii[N_THETA], dr_dtheta[N_THETA];
    const double params[3] = {1.5, 0.1, 0.05};      /* c, a3, a4 */
    int32_t status;

    for (int i = 0; i < N_THETA; ++i)
        thetas[i] = (i + 1) * 3.14159265358979323846 / (N_THETA + 1);

    /* Tier 1: one-shot */
    fos_param_radius_grid(params, 3, thetas, N_THETA, N_POINTS, radii, &status);
    if (status != FOS_VALID) puts(fos_param_status_message(status));

    /* Tier 2: build the cache once, then compute */
    fos_param_cache_t *cache =
            fos_param_cache_create(8, N_POINTS, thetas, N_THETA, &status);
    if (cache == NULL) {
        puts(fos_param_status_message(status));
        return 1;
    }
    status = fos_param_cache_radius_and_derivative(cache, params, 3, radii,
                                                   dr_dtheta, N_THETA);
    fos_param_cache_destroy(cache);
    return status == FOS_VALID ? 0 : 1;
}
```

### Python

```python
import fos_parameterization as fp

params = [1.5, 0.1, 0.05]                    # c, a3, a4
thetas = fp.theta_grid(64)

# Tier 1: one-shot
res = fp.radius_and_derivative(params, thetas, n_points=501)
if not res.ok:
    print(res.message)

# Tier 2: build the cache once, then compute
with fp.Cache(max_params=8, n_points=501, thetas=thetas) as cache:
    for c in (1.2, 1.5, 1.8):
        res = cache.radius_and_derivative([c, 0.1, 0.05])
        shape = cache.shape([c, 0.1, 0.05])
        print(c, res.status.name, shape.z_shift)
```

A validation failure is a result, not an exception: the outputs come back
zero-filled with a nonzero `status`. `FosParamError` is raised only for usage
errors such as a rejected cache construction or a closed handle.

## Status codes

Every fallible routine reports an integer status; `0` is success. On any
nonzero status every output is zero-filled.

| Code | Name | Meaning |
|---|---|---|
| 1 | `too_many_params` | more than 50 parameters |
| 2 | `cache_not_initialized` | compute on a cache that was not built |
| 3 | `invalid_grid` | `n_points < 100`, no thetas, or a theta outside `[0, pi]` |
| 4 | `wrong_param_count` | empty vector, or longer than the cache's `max_params` |
| 5 | `invalid_init` | `max_params < 1` |
| 100 | `rho_negative` | rho <= 0 at an interior node |
| 101 | `not_star_convex` | no origin makes the shape single-valued in R(theta) |
| 102 | `invalid_c` | elongation `c` not positive |
| 103 | `beak_singularity` | `f` too close to zero for the R(theta) conversion |
| 104 | `convergence` | the radius solve did not converge |
| 105 | `buffer_mismatch` | an output array of the wrong size |

Codes 1-5 are shared with the other shape libraries. Inputs must be finite:
the release build cannot detect NaN.

## Build and install

**Python.** Wheels for Linux x86-64 (manylinux2014) are on PyPI:

    pip install fos-parameterization

The source distribution builds anywhere with gfortran and CMake.

**CMake.** Fetch the library into a Fortran or C/C++ project:

```cmake
include(FetchContent)
FetchContent_Declare(
        fos-parameterization
        GIT_REPOSITORY https://github.com/AleksanderAugustyn/fos-parameterization.git
        GIT_TAG 3.0.0
)
set(FOS_PARAM_BUILD_TESTS OFF)
FetchContent_MakeAvailable(fos-parameterization)

target_link_libraries(my_target PRIVATE FosParameterization::fos_parameterization)
```

Targets: `FosParameterization::fos_parameterization` (static, Fortran modules),
`FosParameterization::fos_parameterization_shared` (shared library) and
`FosParameterization::fos_parameterization_cxx` (the C header plus the shared
library).

**From a checkout.**

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure -E geometry_sweep

`geometry_sweep` validates the conversion over half a million shapes and takes
minutes; run it on its own with `ctest --test-dir build -R geometry_sweep`.

## Version pins

| Dependency | Version |
|---|---|
| [gcc-compiler-options](https://github.com/AleksanderAugustyn/gcc-compiler-options) | 2.0.0 |
| [fortran-foundations](https://github.com/AleksanderAugustyn/fortran-foundations) | 3.0.0 |

Both are fetched by CMake. The library needs gfortran 10 or newer; it is
developed with gfortran 13. Python 3.9 or newer, NumPy 1.21 or newer.

## Changelog and license

See the
[changelog](https://github.com/AleksanderAugustyn/fos-parameterization/blob/master/CHANGELOG.md).
MIT license.
