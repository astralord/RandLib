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
| 0.4 Reference data | 🟡 | `tools/refdata/` uses mpmath at 40 digits. Stable α=1 and focused geometric-stable far-tail evidence are committed; exhaustive geometric α=1 generation was time-boxed because the nested integral remains exceptionally expensive. |
| 0.5 Project tracking | 🚫 | Milestones and issues need repository access; this file tracks progress meanwhile. |

## Phase 1: Stabilize (v1.2, branch `1.x`)

| Task | Status | Notes |
|---|---|---|
| 1.1 Repository hygiene | ✅ | |
| 1.2 CMake build | ✅ | Presets, install/export, `find_package`, consumer example built in CI. |
| 1.3 Compiler portability (#18) | ✅ | GCC 11–14, Clang 15–18, AppleClang and MSVC 14.44 build warning-free in the current CI run. |
| 1.4 Correctness fixes | ✅ | B1–B9 as planned and B11–B38 shipped in v1.2; B10 was deferred to v2 by design. The post-release B39–B55 fitting, entropy, parameter-boundary, sampler, covariance, joint-CDF and discrete-moment fixes have landed for the unreleased v1.2.1. See CHANGELOG.md. |
| 1.5 Continuous integration | ✅ | CI run 39018516580 passes B39 on GCC 11–14, Clang 15–18, AppleClang, MSVC, sanitizer, coverage, static analysis, engine variants, portable special functions, and install/consumer checks. |
| 1.6 Test suite (#12) | ✅ | Golden, reference, coherence, mode/median, reproducibility, engines, parameters, regressions, sampling (KS/χ², moments, empirical CF) and fitting. |
| 1.7 Release v1.2 | ✅ | `v1.2.0` tags the all-green tested commit f25fd1a. |

## Phase 2: Modern architecture (v2.0-alpha, branch `2.x`)

| Task | Status | Notes |
|---|---|---|
| 2.1 Header-only layout and namespace | ✅ | `include/randlib/...`, namespace `randlib`, `RANDLIB_PRECOMPILED`, generated ziggurat tables. |
| 2.2 URBG interoperability | ✅ | `operator()(G &)` on every distribution; Pcg32, Pcg64, Xoshiro256pp, SplitMix64, legacy engines, `make_streams`. |
| 2.3 Remove hidden state | ✅ | CI job fails on static or `thread_local` variables that are not constant. |
| 2.4 API cleanup | ✅ | Non-const `Reseed`, `std::span`, `[[nodiscard]]`, `enum class Support`, `fit_error`, `final`, named factories, private type-erased sampling adapters, and internal numerical helpers isolated in `math::detail`. |
| 2.5 Numerical building blocks | ✅ | `math/quadrature.hpp` (adaptive Gauss–Kronrod, infinite limits), `math/roots.hpp` (TOMS 748, Newton, Halley; used by the quantile search), `math/optimize.hpp` (Brent, Nelder–Mead, BFGS), and all special functions are templates on their floating-point type with separate `long double` coverage. |
| 2.6 Tooling | ✅ | Clang 18 formatting and pre-commit, blocking clang-tidy over the explicit instantiations and public headers, and opt-in Google Benchmark sampling comparisons. CI run 39018516566 passes all 17 jobs at 2a2194c; `v2.0.0-alpha.1` tags that exact tested commit. |

## Phase 3: Fitting and features (v2.0-beta)

| Task | Status | Notes |
|---|---|---|
| 3.1 Fitting API | ✅ | `FitResult`, constrained fixed-parameter MLE, AIC/BIC, numerical-Hessian covariance, method-of-moments entry point, conjugate Bayes, model selection and goodness-of-fit are implemented. The free API covers 36 distribution types, including bounded searches for integer-valued parameters and the stable family. Seeded recovery grids and SciPy comparison fixtures are committed, legacy mutating `Fit*` members are deprecated until v3, and `v2.0.0-beta.1` tags the all-green tested commit. |
| 3.2 Tail estimators and extreme-value distributions | ✅ | Generalized Pareto and GEV laws have mpmath references and MLEs; Hill/Pickands diagnostics, peaks-over-threshold fitting, and censored/truncated MLE are implemented. |
| 3.3 Close #11 | 🟡 | The unused private `POISSON` selector is removed; exact Bernoulli-sum, waiting-time and rejection sampling remain. Closing #11 with that rationale needs the maintainer. |

## Phase 4: Documentation, packaging, release (v2.0)

| Task | Status | Notes |
|---|---|---|
| 4.1 Documentation | ✅ | README, strict Doxygen + Pages, fitting/streams/stable/accuracy guides, legacy-source archive, and five CI-compiled examples are done. Run 39021654307 builds and deploys the documentation from `2.x`; the public site serves that deployment. |
| 4.2 Migration guide | ✅ | Every planned breaking change has before/after guidance, including namespace/layout, engines, `Reseed`, spans, fitting, support enums, final classes and reproducibility. |
| 4.3 Packaging | ✅ | The local-source vcpkg overlay and Conan 2 test package are validated; deterministic amalgamation compiles on GCC 13/14 and Clang 18; tag releases attach the header and checksum. The upstream vcpkg PR remains a maintainer action. |
| 4.4 Release v2.0.0 | 🟡 | Release automation is ready, including prerelease classification and RC changelog fallback. The exact-head `v2.0.0-rc.1` tag, two-week soak, announcements and upstream vcpkg PR remain. The already-published alpha/beta GitHub releases need their prerelease flags corrected once by the maintainer. |
| 4.5 Mathematical reference | ✅ | TeX definitions, moments, transforms, support and implementation methods cover all 73 registered public distribution names (71 concrete types plus the Pascal/Polya aliases). The strict docs target runs an executable coverage gate. Plots are intentionally deferred. Stable Count retains its detailed standalone guide on the isolated v2.2 branch and receives a reference entry when that code is merged. |

## Phase 5: Hard numerics (v2.1)

| Task | Status | Notes |
|---|---|---|
| 5.1 Stable distributions (#1) | ✅ | Nolan S0/S1, the skewed near-α=1 probability kernel, TOMS 748 quantiles, McCulloch fitting and constrained MLE are implemented. Eighth-order alpha-one tail expansions remove the geometric-stable bottleneck. The 917-line stable grid and focused former-gap/geometric-tail cases use 40-digit mpmath evidence; exhaustive geometric α=1 generation was time-boxed. |
| 5.2 Heavy-tailed characteristic functions (#6, #10) | ✅ | Type-generic tanh-sinh and QAWF-style oscillatory quadrature are implemented. Positive semi-infinite CFs use half-period integration with Wynn epsilon acceleration; Pareto, Frechet, inverse-gamma, Fisher F, beta-prime and large-variance log-normal references pass the standard `1e-8` absolute target with no overrides. |
| 5.3 Marcum Q (#4) | ✅ | A saddle-centered Poisson mixture evaluates Marcum Q directly without tail cancellation. Noncentral χ² meets the Tier 1 bulk/tail targets with no overrides; the new Tier 2 Rice law has direct probability functions, sampling, moments and 40-digit mpmath references. |
| 5.4 Noncentral t (#7) | ✅ | The v2 port provides direct CDF/survival series, hypergeometric, recurrence and integral densities, large-noncentrality/extreme-tail paths, sampling, quantiles, exact moments, a conditional-normal CF and constrained MLE. Three Tier 1 cases have a 243-value mpmath reference table, supplemented by focused 60-digit integral and exact-series regressions. |
| 5.5 Log-stable (#5) | ✅ | `LogStableRand` provides the exponential stable transform, direct tails, sampling, quantiles, median, Jacobian-correct modes, exact finite moments, S0/S1 handling and log-transformed stable fitting. Exact quarter-turn angles and an analytic finite integration-endpoint limit make the short tail portable on binary64-only platforms. Four Tier 1 cases have a 198-value 40-digit mpmath table. Nonzero CFs deliberately throw `not_implemented` as planned. |

## Phase 6: Stable-count family and stable priors (v2.2)

| Task | Status | Notes |
|---|---|---|
| 6.1 Stable Count distribution | ⏳ | Planned from the published Wright/one-sided-stable definition. The 2021 experimental commit is specification evidence only: its cutoff-based sampler and mutable static state will not be ported. Reference-backed probability functions, exact moments, a standard-URBG sampler and fitting are required. |
| 6.2 Stable-prior semantics | ⏳ | Research gate: identify the exact stable likelihood parameter and prove posterior closure before adding `fit_bayes<StableRand>`. If the cited relation is subordination rather than Bayesian conjugacy, document it as such and keep Stable Count standalone. |
