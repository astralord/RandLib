# Roadmap status

Status of each task in [MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md). This file stands in for
the GitHub milestones and issues of plan task 0.5 until they are created. [NEXT_STEPS.md](NEXT_STEPS.md)
describes how to continue.

Legend: ✅ done · 🟡 partly done (see note) · ⏳ not started · 🚫 needs the maintainer

## Phase 0: Preparation

| Task | Status | Notes |
|---|---|---|
| 0.1 Coordinate upstream | 🚫 | Tracking issue on `astralord/RandLib` is for the maintainer. Branch model agreed: `1.x` for 1.2, `2.x` for v2, `master` untouched. |
| 0.2 Tag the baseline | ✅ | `v1.1.0-baseline` at 6536dad. |
| 0.3 Golden outputs | ✅ | `tools/golden/` records 14.6k values for 160 cases from the baseline library. |
| 0.4 Reference data | 🟡 | `tools/refdata/` (mpmath, 40 digits). The skewed stable and geometric stable laws with α = 1 are still being computed (Phase 5). |
| 0.5 Project tracking | 🚫 | Milestones and issues need repository access; this file tracks progress meanwhile. |

## Phase 1: Stabilize (v1.2, branch `1.x`)

| Task | Status | Notes |
|---|---|---|
| 1.1 Repository hygiene | ✅ | |
| 1.2 CMake build | ✅ | Presets, install/export, `find_package`, consumer example built in CI. |
| 1.3 Compiler portability (#18) | 🟡 | GCC 11–14, Clang 15–18 and AppleClang build warning-free in CI; the MSVC 14.44 duplicate-instantiation fix awaits CI. |
| 1.4 Correctness fixes | ✅ | B1–B9 as planned, B10 deferred to v2 by design, and B11–B37 found by the test suite; see CHANGELOG.md. |
| 1.5 Continuous integration | 🟡 | GCC, Clang, sanitizer, coverage and static-analysis jobs pass; the MSVC 14.44 build fix and final AppleClang statistical result await CI. |
| 1.6 Test suite (#12) | ✅ | Golden, reference, coherence, mode/median, reproducibility, engines, parameters, regressions, sampling (KS/χ², moments, empirical CF) and fitting. |
| 1.7 Release v1.2 | 🟡 | CHANGELOG written; `v1.2.0` is tagged once CI is green on all compilers. |

## Phase 2: Modern architecture (v2.0-alpha, branch `2.x`)

| Task | Status | Notes |
|---|---|---|
| 2.1 Header-only layout and namespace | ✅ | `include/randlib/...`, namespace `randlib`, `RANDLIB_PRECOMPILED`, generated ziggurat tables. |
| 2.2 URBG interoperability | ✅ | `operator()(G &)` on every distribution; Pcg32, Pcg64, Xoshiro256pp, SplitMix64, legacy engines, `make_streams`. |
| 2.3 Remove hidden state | ✅ | CI job fails on static or `thread_local` variables that are not constant. |
| 2.4 API cleanup | ✅ | Non-const `Reseed`, `std::span`, `[[nodiscard]]`, `enum class Support`, `fit_error`, `final`, named factories, private type-erased sampling adapters, and internal numerical helpers isolated in `math::detail`. |
| 2.5 Numerical building blocks | ✅ | `math/quadrature.hpp` (adaptive Gauss–Kronrod, infinite limits), `math/roots.hpp` (TOMS 748, Newton, Halley; used by the quantile search), `math/optimize.hpp` (Brent, Nelder–Mead, BFGS), and all special functions are templates on their floating-point type with separate `long double` coverage. |
| 2.6 Tooling | ✅ | Clang 18 formatting and pre-commit, blocking clang-tidy over the explicit instantiations and public headers, and opt-in Google Benchmark sampling comparisons. |

## Phase 3: Fitting and features (v2.0-beta)

| Task | Status | Notes |
|---|---|---|
| 3.1 Fitting API | 🟡 | `FitResult`, constrained fixed-parameter MLE, AIC/BIC, numerical-Hessian covariance, method-of-moments entry point, conjugate Bayes, model selection and goodness-of-fit are implemented. The free API covers all 31 non-stable Tier 1 types plus GPD and GEV, including bounded searches for integer-valued parameters. Seeded recovery grids and SciPy comparison fixtures are committed. Remaining: stable-family fits (Phase 5) and deprecated member wrappers. |
| 3.2 Tail estimators and extreme-value distributions | ✅ | Generalized Pareto and GEV laws have mpmath references and MLEs; Hill/Pickands diagnostics, peaks-over-threshold fitting, and censored/truncated MLE are implemented. |
| 3.3 Close #11 | 🟡 | The unused private `POISSON` selector is removed; exact Bernoulli-sum, waiting-time and rejection sampling remain. Closing #11 with that rationale needs the maintainer. |

## Phase 4: Documentation, packaging, release (v2.0)

| Task | Status | Notes |
|---|---|---|
| 4.1 Documentation | ⏳ | |
| 4.2 Migration guide | ⏳ | |
| 4.3 Packaging | ⏳ | |
| 4.4 Release v2.0.0 | 🚫 | Release candidates, announcements and the upstream vcpkg PR need the maintainer. |

## Phase 5: Hard numerics (v2.1)

| Task | Status | Notes |
|---|---|---|
| 5.1 Stable distributions (#1) | 🟡 | v1.2 fixed the CF sign, the α = 1 location, the B8 floor and the tails. Open: S0, the region near α = 1, the slow far tails of the α = 1 geometric stable density. |
| 5.2 Heavy-tailed characteristic functions (#6, #10) | 🟡 | Heavy-tailed CFs are accurate to about 1e-3, the log-normal one to 6e-5 at σ² = 4 (tracked in `tests/data/reference_tolerances.csv`). |
| 5.3 Marcum Q (#4) | 🟡 | v1.2 fixed the noncentral χ² density; the upper tail is still 1 − Marcum P. |
| 5.4 Noncentral t (#7) | ⏳ | |
| 5.5 Log-stable (#5) | ⏳ | |
