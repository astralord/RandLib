# Next steps: handoff for continuing the modernization

Written so that work can resume with no other context. Read together with
[MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md) (the plan, revised on 2026-10-17) and
[ROADMAP_STATUS.md](ROADMAP_STATUS.md) (status per task). State as of 2026-10-17.

## 1. Repository, branches and how changes reach GitHub

- Upstream: `github.com/astralord/RandLib`, owned by the maintainer, who drives this work. `master`
  (6536dad, tagged `v1.1.0-baseline`) must never be pushed to.
- Local clone: `/root/RandLib`. Its checked-out branch `1.x` is **stale** (6536dad plus old uncommitted
  copies of the work); do not use it and do not reset it without asking (resetting was refused once).
- Worktrees:

  | Path | Branch | Pushed as | Contents |
  |---|---|---|---|
  | `/root/RandLib-1.x` | `history` | `origin/1.x` | v1.2: C++17, original layout, fixes, CMake, tests, CI |
  | `/root/RandLib-2.x` | `2.x` | `origin/2.x` | v2: header-only C++20 (Phases 2–5 landed; Phase 6 planned) |

- The agent's GitHub account (`nebius-llm-robot`) has pull access only. **The maintainer pushes**, from a
  VS Code terminal on the pod:

  ```sh
  cd /root/RandLib
  git push origin history:refs/heads/1.x 2.x
  ```

  Only add commits on top of what is pushed; never amend or rebase pushed commits.
- 1.x fixes are made on `history` first and then merged into `2.x` (`git merge --no-ff history` in
  `/root/RandLib-2.x`). Conflicts are expected on files that 2.x deleted (`distributions/`, `math/`):
  `git rm` them, and port the change into the matching header under `include/randlib/` by hand (git's
  rename detection often does it already). Data files (`tests/data`, `tools/golden/reviewed_deviations.csv`)
  usually just need both sides kept.
- CI: `.github/workflows/ci.yml` on each branch. Read results with
  `gh run list -R astralord/RandLib` and `gh run view <id> -R astralord/RandLib --log-failed`.

## 2. Building and testing

```sh
# 2.x (header-only); Catch2 is fetched by CMake
cmake -S /root/RandLib-2.x -B /tmp/b2 -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-14 -DRANDLIB_WARNINGS_AS_ERRORS=ON
cmake --build /tmp/b2
/tmp/b2/tests/randlib_unit_tests          # ~10 s
/tmp/b2/tests/randlib_statistical_tests   # ~25 min on the current pod
# 1.x: the same with -S /root/RandLib-1.x
```

Compilers on the pod: g++-13, g++-14, clang++-18. MSVC and AppleClang are only exercised by CI.

Test families (both branches): `tests/unit` (golden v1.1 comparison, mpmath reference values, coherence,
mode/median, parameters, regressions B1–B39, reproducibility, engines, numerics on 2.x) and
`tests/statistical` (KS/χ² sampling, moments, empirical CF, fitting). Tolerances follow plan section 9.

### Golden values and known deviations

`tests/data/golden_v1/*.csv` are v1.1 outputs. Every value that legitimately changed is listed in
`tests/data/golden_v1/known_deviations.csv`, generated from hand-reviewed rules plus values verified
against the mpmath references:

```sh
RANDLIB_GOLDEN_REPORT=/tmp/mismatches.csv <build>/tests/randlib_unit_tests "[golden]"
tools/golden/classify_deviations.py /tmp/mismatches.csv > tests/data/golden_v1/known_deviations.csv
```

The report contains only rows not covered by `known_deviations.csv`, so for a full regeneration first
reduce that file to the reviewed rules (copy `tools/golden/reviewed_deviations.csv` without its `#`
lines). Small additions can be inserted by hand after the reviewed rules. Every new rule needs evidence
in its reason (mpmath value, exact formula, regression test).

### Reference values

`tools/refdata/generate.py` (mpmath, 40 digits, runs for hours) writes `tests/data/reference/*.csv`
and caches per case in `/root/RandLib/build/refdata-cache`. The run started on 2026-10-10 is no
longer running. Exhaustive skewed α=1 stable and geometric-stable generation was time-boxed; focused
40-digit former-gap and geometric-tail references were committed instead. The cache can be reused if
the exhaustive grid is revisited, but no reference job is currently running and this is not a release
gate.

## 3. Immediate next actions

### 3.1 v1.2 release and 1.2.1 follow-up (branch `history` → `origin/1.x`)

CI run 39002870454 passes all 16 jobs at f25fd1a, and the pushed `v1.2.0` tag points to that exact
tested commit. Commit ba0592b is the first post-release fix: B39 keeps negative-binomial fitting
inside the positive shape domain and is recorded under the unreleased 1.2.1 changelog section.
Replacement run 39018516580 passes B39 in all 16 jobs, including AppleClang, MSVC, GCC 11–14,
Clang 15–18, sanitizers, engine variants, portable special functions, static analysis and the
install consumer. The subsequent 4ad71f0 commit only records that result.

### 3.2 Phase 2 remainder (branch `2.x`)

- **CI portability**: pushed commits fix the Clang/AppleClang `-Wshadow` errors in `Pcg64`, discarded
  `[[nodiscard]]` results, MSVC's variadic-macro, secure-CRT, unsigned-negation, uninitialized-local
  and unused-variable diagnostics, the coverage parser error, Apple libc++'s missing `std::jthread`,
  and a GCC 12 nested-vector false positive. The engine-adapter test checks exact engine advancement
  while allowing eight ULPs of platform-dependent floating-point variation; AppleClang measured a
  maximum of six. CI run 38957273851 at pushed head 2d038bd passes every job except AppleClang and
  MSVC unit tests, which agree on one log-stable `1e-10` short-tail quantile. The local head reproduces
  their exact result with GCC's binary64-long-double mode and fixes it with exact quarter-turn
  precomputations plus the analytic finite limit of the totally skewed stable integrand at π/2.
  Native and binary64-long-double probes now agree with the committed 40-digit mpmath reference.
  Warning-clean Release builds and the complete 168-test matrix pass under GCC 13 and Clang 18;
  GCC 14 passes the unit/example and isolated-discrete suites. The expanded v2 continuous sampling
  grid keeps a one-hour bound, while every discrete case remains bounded at 60 seconds. The B38
  sampler fix from 1.x, including checked negative-binomial accumulation, is ported into the
  header-only implementations and uses a standard engine in its v2 regression. CI run 39002870514
  passes every job except MSVC: its deterministic negative-binomial MLE crossed the positive shape
  boundary, spent 36 minutes in Newton iteration, then threw. B39 landed on `history` first and is
  merged into `2.x`: the fitter works in log-shape coordinates, and Newton rejects non-finite states
  after bounded backtracking. Replacement CI run 39018516566 passes all 17 jobs at 2a2194c,
  including the formerly failing MSVC fit, AppleClang, GCC 11–14, Clang 15–18, sanitizer, blocking
  clang-tidy, coverage, precompiled/single-header, portable-special-function, no-static-state and
  install-consumer lanes.
- **2.5 Numerical building blocks**: done. `math/quadrature.hpp`, `math/roots.hpp`,
  `math/optimize.hpp`, and the functions in `math/functions.hpp`, `gamma.hpp`, `beta.hpp` and
  `special_functions.hpp` are templates on their floating-point type. Mixed integer/floating call
  sites preserve the floating type, and `tests/unit/test_numerics.cpp` independently instantiates
  and checks the `long double` path with exact identities and closed forms.
- **2.6 Tooling**:
  - Done: `.clang-format` matching the existing style, a pinned clang-format pre-commit hook,
    and one formatting pass over `include/` and `tests/`. Include sorting is deliberately disabled
    because the converted headers still contain dependency-sensitive implementation include blocks.
  - Done: clang-tidy analyzes `src/precompiled.cpp`, whose explicit instantiations cover the
    public distribution headers, and promotes every diagnostic to an error in blocking CI.
  - Done: opt-in `benchmarks/` with pinned Google Benchmark 1.9.5. Configure with
    `-DRANDLIB_BUILD_BENCHMARKS=ON`; it compares normal, exponential, gamma and Poisson sampling
    against `std::*_distribution` and measures stable sampling. Benchmarks are tracked, not gated.

    ```sh
    cmake -S /root/RandLib-2.x -B /tmp/randlib-bench -DCMAKE_BUILD_TYPE=Release \
      -DRANDLIB_BUILD_TESTS=OFF -DRANDLIB_BUILD_BENCHMARKS=ON
    cmake --build /tmp/randlib-bench --target randlib_benchmarks
    /tmp/randlib-bench/benchmarks/randlib_benchmarks
    ```
- **Clean-ups left by the automated header conversion** (worth doing with or before 2.6):
  - Done: `variate(detail::EngineRef)` and bulk `sample` overrides are private adapters at
    the start of each class; the redundant access-specifier toggles are removed.
  - Done: implementation helpers moved from the public `randlib::math` namespace into
    `randlib::math::detail`, including the special-function iteration constants, `erfinvAux*`,
    Lambert W, Marcum P/Q and portable Bessel helpers.
  - Done: `areClose` and the Newton tolerance are private math details; the public
    `Pair`/`Triplet` aliases moved from `math/constants.hpp` to `core/fwd.hpp`; the global pair
    arithmetic overloads were removed in favor of an explicit exponential-family inner product.
- The annotated `v2.0.0-alpha.1` tag points to the exact all-green commit 2a2194c. Because Phases
  3–5 landed while the cross-platform gate was being repaired, this is a delayed prerelease
  checkpoint rather than the historical end of the Phase 2 commits.

### 3.3 Phase 3: fitting and features (v2.0-beta)

Follow plan 3.1–3.3. Design notes gathered so far:
- Landed: `math::numericalHessian` and the first `include/randlib/fit/` vertical slice.
  `FitResult`, `FitOptions`, fixed-parameter constrained MLE, AIC/BIC and observed-information
  covariance are covered. The free API covers 36 distribution types, including fixed structural
  parameters, generic multi-parameter MLEs and bounded searches for integer-valued parameters.
  The free API now also has a named method-of-moments entry point, conjugate Bayesian updates,
  AIC/BIC model selection, and continuous/discrete goodness-of-fit tests. Conditional Bayesian
  overloads take the model's fixed parameter as their third argument (variance for normal,
  variance of the logarithms for log-normal, trial count for binomial, shape/count for
  negative-binomial, shape for Weibull, scale for Pareto, and minimum for uniform). Stable-family
  fitting landed in Phase 5. Seeded recovery grids cover representative closed-form, generic and
  discrete fits at two sample sizes. `tools/fitting/generate_scipy_reference.py` regenerates the
  fixed datasets and SciPy MLE comparisons; install the pinned developer dependencies with
  `python3 -m pip install -r tools/requirements.txt`. The legacy mutating member functions (`Fit`,
  `FitScale`, `FitShapeBayes`, and related overloads) now emit compiler deprecation diagnostics that
  point to the non-mutating free API. They remain source-compatible until v3; the temporary
  `RANDLIB_DISABLE_DEPRECATED_FIT_WARNINGS` opt-out exists for compatibility test suites. Stable-law
  fitting remains explicitly assigned to Phase 5.
- New free API in `include/randlib/fit/` (`result.hpp`, `mle.hpp`, `moments.hpp`, `bayes.hpp`,
  `gof.hpp`, `select.hpp`): `fit<D>(data, Method)`, `FitResult<D>` with log-likelihood, AIC/BIC,
  standard errors and covariance from a numerical Hessian (add `math::numericalHessian` next to
  `numericalGradient` in `math/optimize.hpp`), `fit_bayes`, `fit_best<D...>(data, Criterion)`,
  `goodness_of_fit` (KS via `KolmogorovSmirnovRand`, Anderson–Darling, Cramér–von Mises, χ² for
  discrete).
- Generic MLE: a per-distribution parameter description (names, bounds, getter/setter) plus
  `math::nelderMead` / `math::bfgs` on an unconstrained transform (log for positive parameters,
  logit for (0, 1)). Closed-form fits override it.
- Tests: parameter recovery over a grid (seeded), comparison to SciPy `fit` results on fixed datasets
  stored as CSV (generated by `tools/fitting/generate_scipy_reference.py`).
- 3.2 is done: `fit/tail.hpp` provides Hill and Pickands estimators, plot-ready Hill
  diagnostics and peaks-over-threshold fitting. Generalized Pareto and generalized extreme-value
  distributions have closed-form probability functions, sampling, moments, generic MLEs and
  committed 40-digit mpmath reference data. `fit/censored.hpp` supports exact, left-, right- and
  interval-censored observations plus one- and two-sided truncated MLE through the shared
  constrained-likelihood pipeline.
- 3.3 implementation is done: the unused private binomial `POISSON` selector is removed on both
  branches. The maintainer still needs to close #11 with the exact-sampling rationale.
- The beta changelog section is separate from the alpha notes. The annotated `v2.0.0-beta.1` tag
  points to the same all-green commit 2a2194c; all Phase 3 engineering precedes it.

### 3.4 Phase 4: documentation, packaging, release (v2.0)

Landed: the README rewrite and comparison table; strict Doxygen Awesome generation and Pages
deployment; fitting, streams, stable-law, accuracy, packaging and migration guides; and the archived
legacy LaTeX source with its generated PDF removed. Five warning-clean examples (standard engine,
parallel streams, model selection, stable fitting and Bayesian update) build and run through CTest.

Packaging is also landed and locally validated. `packaging/vcpkg/randlib` installs the current checkout
as an overlay and its installed `RandLib::RandLib` target built the consumer example with vcpkg
2026-09-26. `conan create . --build=missing` passed with Conan 2.33.0, packaging all 95 headers and
running `test_package`. `tools/amalgamate.py` produces a deterministic `randlib_single.hpp`; the
941,344-byte result compiled and ran with g++-13, g++-14 and clang++-18 under warnings-as-errors.
The release workflow compiles that artifact, extracts the exact tag's changelog section, and attaches
the header and SHA-256 checksum. The upstream vcpkg port must replace the checkout source with the
final release archive and checksum.

Documentation run 39021654307 builds and deploys successfully from `2.x`. The public site at
<https://astralord.github.io/RandLib/> serves that deployment, so task 4.1 is complete.

The first tag runs exposed two release-workflow gaps: alpha and beta were published as ordinary
GitHub releases, and an RC tag had no exact changelog heading. The local release-tooling fix marks
every hyphenated SemVer tag as a prerelease, lets an RC reuse the final version's changelog section,
and adds focused extraction tests to CI. After that fix is pushed and the exact head passes the full
matrix, create annotated tag `v2.0.0-rc.1`. Do not start the two-week soak from the earlier alpha or
beta tags.

The maintainer must correct the two existing GitHub release records once; this does not move tags or
rewrite history:

```sh
gh release edit v2.0.0-alpha.1 -R astralord/RandLib --prerelease
gh release edit v2.0.0-beta.1 -R astralord/RandLib --prerelease
```

The remaining Phase 4 actions are the release-candidate soak, upstream vcpkg PR and release
announcements.

### 3.5 Phase 5: hard numerics (v2.1), deliberately last

The Phase 5 engineering work is complete. The remaining actions are the maintainer-owned issue
closures, CI observation and releases listed below.

Landed for 5.1: `StableRand::Parameterization::{S0,S1}` with S1 as the source-compatible
default. The skewed near-α=1 density/CDF gap is implemented and checked against 40-digit mpmath
values; TOMS 748 quantiles, McCulloch quantile fitting and constrained MLE are covered by seeded
recovery tests. Eighth-order alpha-one tail expansions make the complete geometric-stable coherence
grid practical; two 40-digit mpmath tail points are committed and exhaustive generation was
time-boxed. Phase 5.1 is complete.
Phase 5.2 is complete. `math::integrateTanhSinh` handles finite-interval endpoint singularities, and
`math::integrateOscillatory` combines half-period integration with Wynn epsilon acceleration. The
heavy-tailed and log-normal CF reference grids pass the standard `1e-8` absolute target; their
twelve tolerance overrides were removed.
Phase 5.3 is complete. Marcum Q is evaluated directly as a Poisson mixture centered at the combined
Poisson/incomplete-gamma saddle point, so representable noncentral χ² survival probabilities remain
accurate down to `1e-300`. Its three tolerance overrides were removed. `RiceRand` uses the same
kernel and has direct probability functions, transformed-normal sampling, closed-form mean/variance,
and a committed 237-value mpmath reference table.
Phase 5.4 is complete. `NoncentralTRand` uses centered Benton–Krishnamoorthy series for direct lower
and upper probabilities, recurrence and conditional-normal integral densities, transformed-normal
sampling, generic quantiles, exact moments, a conditional-normal characteristic function and
constrained MLE. Saddle-split integration and a small-denominator expansion cover large
noncentrality and extreme representable tails. A committed 243-value mpmath table covers three Tier
1 cases; focused regressions use direct 60-digit integrals and an exact normal-moment tail series.
The moderate-argument hypergeometric density avoids probability subtraction when `long double` is
binary64.
Phase 5.5 is complete. `LogStableRand` is the exponential transform of `StableRand`, preserving the
selected S0/S1 convention. It provides transformed density/tails, sampling, quantiles, median,
Jacobian-correct modes, stable-family quantile/MLE fitting and exact Laplace-exponent moments where
they exist (`beta=-1`, plus the Gaussian case). Other non-Gaussian skewnesses report divergent
moments, and the deliberately unsupported nonzero characteristic function throws
`not_implemented`. Four Tier 1 cases have a committed 198-value 40-digit mpmath table. Exact
quarter-turn angle reduction and an analytic finite integration-endpoint limit keep the totally
skewed short-tail quantiles portable on binary64-only platforms.

### 3.6 Phase 6: Stable Count and stable-prior research (v2.2)

The maintainer requested Stable Count as the next planned stable-family distribution. The starting
evidence is experimental commit
[`c070098`](https://github.com/astralord/RandLib/commit/c0700988a88767cde9967b41e0cc459b8434cc6a)
and Lihn's revised [Wright-function paper](https://arxiv.org/abs/2405.04693). The commit must not be
ported mechanically: its sampler stops a positive-stable renewal sum at the arbitrary value
`100000`, discards stable draws above `100000`, and uses mutable static state. It also carries the
old `[0.1, 1)` and near-one limitations that Phase 5 removed from the stable kernel.

Plan task 6.1 now specifies a fresh `StableCountRand` implementation from
`N_alpha(x) = L_alpha(1/x) / (Gamma(1/alpha + 1) x)`, with explicit location/scale semantics,
direct tails, quantiles, exact log-gamma moments, the `alpha = 1/2` Gamma reduction, a validated
standard-URBG sampler, 40-digit mpmath data and generic MLE. It starts in Tier 2 and is promoted only
after its independent numerical and sampling evidence satisfies the Tier 1 policy. The paper's
four-parameter generalized Stable Count family is a later proposal, not part of this first slice.

Plan task 6.2 treats the phrase "conjugate prior" as a research claim, not yet an API contract. Before
adding `fit_bayes<StableRand>`, write the normalized likelihood, identify whether the prior is on a
count, count root, scale or reciprocal scale, and prove the posterior stays in the Stable Count
family for arbitrary sample size. The paper also uses inverse-variable, product-mixture and
subordination identities; those are useful, but they are not by themselves Bayesian conjugacy. If
posterior closure is not established within the time box, ship Stable Count standalone and record
the missing derivation instead of exposing a misleading closed-form update.

## 4. Things that need the maintainer

Plan 0.1 and 0.5 (tracking issue, milestones), closing issues (#11 and others as they are fixed),
pushing branches and tags, release announcements, and the upstream vcpkg PR.
