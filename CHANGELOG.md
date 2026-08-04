# Changelog

Notable changes to fos-parameterization. Versions follow semantic versioning;
the Fortran, C and Python surfaces are versioned together.

## 2.0.0 — 2026-08-04

Adoption of the shared shape-parameterization contract
(fortran-foundations 2.4.0 `shape_core_mod`, template beta-parameterization
3.0.0). Every public surface changed. **No 1.x consumer compiles unchanged**,
and several 1.x call sites that compiled would now be silently wrong if the
symbols had been kept, so nothing was kept as a deprecated alias.

Migration is a rewrite of the call sites, not a rename pass: read the module
header of `fos_parameterization_mod` (tiers, lifetime rules, dependency map)
before porting.

### Breaking — Fortran surface

The flat `compute_fos_*_s` surface is gone. Every shape output now comes from
one of three tiers: **tier 2** (`cache_t`, incremental, ≤ 8 parameters),
**tier 1** (`*_standalone_s`, stateless, ≤ 50 parameters), and shared
`tables_t` bound to many caches.

| 1.x | 2.0.0 cached (tier 2) | 2.0.0 standalone (tier 1) |
|---|---|---|
| `compute_fos_radius_grid_s` | `cache_radius_grid_s` | `compute_radius_grid_standalone_s` |
| `compute_fos_radius_and_derivative_s` (elemental core) | `cache_radius_and_derivative_s` | `compute_radius_and_derivative_standalone_s` |
| `compute_fos_shape_s` | `cache_shape_s` | `compute_shape_standalone_s` |
| `compute_rho_z_grid_s` (+ `rho_z_grid_t`) | `cache_rho_z_grid_s` | `compute_rho_z_grid_standalone_s` |
| `compute_fos_neck_s` | `cache_neck_s` | `compute_neck_standalone_s` |
| `compute_fos_star_convexity_optimum_s` | `cache_star_convexity_optimum_s` | `compute_star_convexity_optimum_standalone_s` |
| `compute_fos_radius_and_derivative_at_thetas_s` | `cache_radius_and_derivative_at_thetas_s` | — (per-call thetas are a cached extension) |
| `compute_fos_a2_f` (function) | `compute_a2_s(params, a2, status)` | same |
| `compute_fos_z_shift_f` (function, sentinel 0.0) | `compute_z_shift_s(params, z_shift, status)` | same |

New: `tables_t`, `tables_init_s`, `tables_free_s`, `cache_t`, `cache_init_s`,
`cache_init_shared_s`, `cache_free_s`, `cache_recompute_count_f`,
`status_message`, `FOS_N_POINTS_FLOOR` (100), `FOS_MAX_PARAMS` (50).

Privatized (their verdicts now surface only as status codes):
`validate_rho_grid_s`, `check_star_convexity_s`, `make_fos_shape_f`,
`fos_shape_t`, `rho_z_grid_t`, `compute_radius_at_theta_s`,
`compute_radius_fos_with_zshift_s`.

Kept public and unchanged, now documented as **raw evaluators outside the tier
rules** (no validation, no status, documented fallbacks, zero-extension of
short vectors): `compute_fos_f_and_derivatives_s`, `get_fos_coefficient_f`,
`compute_rho_at_z_s`.

`message` arguments are removed from every routine. Use
`status_message(code)`, which returns a fixed 64-character string.

### Breaking — status codes renumbered

The 0–99 range now belongs to the shared contract, so every FoS-specific code
moved into 100+.

| 1.x | 2.0.0 | name |
|---|---|---|
| 0 | 0 | `SHAPE_VALID` (`FOS_VALID` kept as an alias) |
| 1 | **100** | `FOS_ERROR_RHO_NEGATIVE` |
| 2 | **101** | `FOS_ERROR_NOT_STAR_CONVEX` |
| 3 | **102** | `FOS_ERROR_INVALID_C` |
| 4 | **103** | `FOS_ERROR_BEAK_SINGULARITY` |
| — | **104** | `FOS_ERROR_CONVERGENCE` (new) |
| — | **105** | `FOS_ERROR_BUFFER_MISMATCH` (new) |

Shared codes re-exported from `shape_core_mod`, all new to this library:
1 `SHAPE_ERROR_TOO_MANY_PARAMS`, 2 `SHAPE_ERROR_CACHE_NOT_INITIALIZED`,
3 `SHAPE_ERROR_INVALID_GRID`, 4 `SHAPE_ERROR_WRONG_PARAM_COUNT`,
5 `SHAPE_ERROR_INVALID_INIT`, 6 `SHAPE_ERROR_TABLES_NOT_INITIALIZED`.

The renumbering is silent for any consumer that compares against literals:
1.x's "rho ≤ 0" (1) is now "too many parameters", and 1.x's "invalid c" (3) is
now "invalid theta grid". Compare against the named constants.

100+ is append-only from this release on.

### Breaking — behavior at the gates

- **Gate precedence.** A shape that is both ρ-negative and beak-failing now
  returns **103**; 1.x returned 100 (its code 1), because 1.x scanned ρ before
  f_min and 2.0.0 checks the beak verdict first. 100 is still reachable: on the
  cylindrical path (`rho_z_grid`, `neck`), where the beak never gates, and on
  the R(θ) path when `n_points` > 1001, since the ρ scan then samples nodes the
  fixed 1001-node beak scan does not.
- **Beak boundary.** A shape whose f_min is exactly `F_MIN_THRESHOLD` (1e-3) is
  now rejected: the accept condition is `f_min > F_MIN_THRESHOLD`, where 1.x
  rejected on `f_min < F_MIN_THRESHOLD` and accepted equality. Measure-zero in
  practice, listed for completeness.
- **Gating asymmetry (new capability).** Beak and star-convexity exist only to
  guarantee the R(θ) conversion, so they no longer gate cylindrical output. A
  beak-marginal or non-star-convex shape renders through `rho_z_grid` and
  `neck` while the R(θ) and shape outputs reject it.
- **Oversize parameter vectors.** `fos_param_a2` / `fos_param_z_shift` (Fortran
  `compute_a2_s` / `compute_z_shift_s`) return status 1
  (`SHAPE_ERROR_TOO_MANY_PARAMS`) for more than 50 parameters. 1.x returned a
  **truncated sum** with no indication. Cached calls require
  `size(params)` exactly equal to the cache's `n_params` (status 4 otherwise);
  standalone calls cap at 50; nothing is ever truncated.
- **Degenerate parameters.** The 1.x silent unit-sphere fallback is withdrawn:
  an empty vector or c ≤ `C_MIN` is status 102 with zero-filled outputs,
  uniformly across every entry point.
- **Rejections zero-fill.** Any rejection zero-fills the entire actual extent of
  every output array and invalidates the cache, so the next call runs cold. No
  partially-updated state, no cached rejections.
- **`status_message` returns fixed strings.** 1.x formatted diagnostics with the
  offending values embedded (`'... f_min = 3.2E-04 ...'`). 2.0.0 returns one
  static string per code — parseable, allocation-free, and safe to return
  through the C API as a `const char *`.

### Fixed

- **Newton radius solves use an analytic bracket.** 1.x bracketed from twice the
  pole extents with 8 doublings and returned the last iterate on failure. An
  extreme-oblate but perfectly valid shape (c = 2e-10, equatorial radius
  √(1/c) ≈ 7.07e4) escaped the bracket and returned ≈ 1e-7 — a wrong radius
  with status 0. 2.0.0 brackets at 2·√(ρ_max² + max(z_max², z_min²)) from cached
  state, so the sign change is always enclosed, and reports **104** if any node
  fails to converge. A wrong value with a success status is no longer
  reachable through the public surface.
- **C API marshalling buffers are heap-allocated.** The 1.x/beta idiom of stack
  automatics could overflow the stack and SIGSEGV in Release builds when a
  caller declared a large size. Buffer sizes are now validated against the
  cache's init-time resolution first (**105** on mismatch) and the marshalling
  arrays are allocatable.

### Breaking — thread and lifetime model

- One `cache_t` per thread. A cache is mutated by every compute call.
- `tables_t` is immutable after `tables_init_s` and may be shared by any number
  of threads' caches through `cache_init_shared_s`.
- A shared `tables_t` must be declared `target` and must outlive every cache
  bound to it (the cache stores a pointer, not a copy). Free the caches first.
- Never copy-assign a `cache_t`: intrinsic assignment copies the owned tables
  pointer shallowly and the copy double-frees it.
- `cache_free_s` / `tables_free_s` are mandatory — there is no finalizer.

### Breaking — C API

Every symbol is renamed `fos_*` → `fos_param_*` and reshaped into
`<prefix>_<noun>_<verb>`:

| 1.x | 2.0.0 |
|---|---|
| `fos_compute_radius_grid` | `fos_param_cache_radius_grid` / `fos_param_radius_grid` |
| `fos_compute_rho_profile` | `fos_param_cache_rho_z_grid` / `fos_param_rho_z_grid` |
| `fos_compute_neck` | `fos_param_cache_neck` / `fos_param_neck` |
| `fos_compute_shape` | `fos_param_cache_shape` / `fos_param_shape` |
| `fos_compute_radius_and_derivative_at_thetas` | `fos_param_cache_radius_and_derivative_at_thetas` |
| `fos_z_shift`, `fos_a2` (returned a double, sentinel on failure) | `fos_param_z_shift`, `fos_param_a2` (out-params + nullable `int32_t *status`) |

New: `fos_param_tables_create/_destroy`, `fos_param_cache_create`,
`fos_param_cache_create_shared`, `fos_param_cache_destroy`,
`fos_param_cache_radius_and_derivative`,
`fos_param_cache_star_convexity_optimum`,
`fos_param_radius_and_derivative`, `fos_param_star_convexity_optimum`,
`fos_param_rho_at_z` (a genuinely new entry point — 1.x exported no
point-evaluator), `fos_param_status_message`.

Handle-creating calls return `NULL` on failure. Cached compute calls return the
status directly and take explicit buffer sizes; a size that does not match the
cache's init-time resolution returns 105 with zero-filled outputs, so no write
ever exceeds the declared capacity. Cache handles are thread-confined; tables
handles are shareable.

### Breaking — Python

The wheel exposes the full 2.0 surface. `RhoProfileResult` and `ShapeResult`
are gone; `rho_profile` is now `rho_z_grid`.

- `Cache` class (`radius_grid`, `radius_and_derivative`,
  `radius_and_derivative_at_thetas`, `shape`, `rho_z_grid`, `neck`,
  `star_convexity_optimum`), context-manager capable, thread-confined.
- Module-level tier-1 functions: `radius_grid`, `radius_and_derivative`,
  `shape`, `rho_z_grid`, `neck`, `star_convexity_optimum`, `rho_at_z`,
  `z_shift`, `a2`, `theta_grid`, `status_message`.
- Result dataclasses with `.ok` / `.message`: `RadiusGridResult`,
  `RadiusDerivativeResult`, `FosShape`, `RhoZGridResult`, `NeckResult`,
  `StarConvexityResult`.
- Complete `Status` enum — shared codes 1–6 as well as 0 and 100–105.
- Validation failures return zero-filled result objects. Exceptions
  (`FosParamError`) are raised only for handle-create failure, use of a closed
  handle, and non-1-D input.

### Consumer trap — theta grids and `-ffast-math`

θ nodes must lie in [0, π]; the evaluator reconstructs
sin θ = +√(1 − cos²θ), so anything outside silently aliases back into the
domain and the library rejects it with status 3 instead. A closed grid built as
`(i-1)·π/(n-1)` can land one ulp **above** π under `-ffast-math` and be
rejected. Pin the endpoints explicitly, or use an open grid such as the wheel's
`theta_grid(n)` = i·π/(n+1). Recorded as atlas item
`2026-08-04-fos-theta-grid-ulp-overshoot`.

### Internal

- `n_points` floor of 100 enforced at tables init (`SHAPE_ERROR_INVALID_GRID`).
- Cached tier caps at 8 parameters (`SHAPE_ERROR_TOO_MANY_PARAMS`); tier 1 at
  50, where all 50 slots demonstrably influence the surface.
- Eight tracked intermediates with a normative dependency map (module header);
  `cache_recompute_count_f` exposes the recompute counters.
- Tests: contract families (bitwise warm ≡ cold, minimality, interleaving,
  boundaries), gating asymmetry, thread stress, a dual-configuration pytest
  suite run against the installed wheel, and a PES-box sweep gate asserting
  incremental ≡ cold bit-for-bit across ~9,000 shapes with every cached slot
  moving.
- Bitwise goldens re-captured once at adoption: 2.0.0 evaluates f(u) from the
  tabled trig basis rather than live per-call trig, so 1.x bit patterns do not
  carry over. Physical-invariant tolerances are unchanged.
- Pins: GCC-Compiler-Options 2.0.0, Fortran-Foundations 2.4.0 (was 2.3.2).

## 1.3.2 and earlier

Not tracked in this file. See the git history.
