# Next steps: handoff for continuing the modernization

Written so that work can resume with no other context. Read together with
[MODERNIZATION_PLAN.md](MODERNIZATION_PLAN.md) (the plan, revised on 2026-10-11) and
[ROADMAP_STATUS.md](ROADMAP_STATUS.md) (status per task). State as of 2026-10-16.

## 1. Repository, branches and how changes reach GitHub

- Upstream: `github.com/astralord/RandLib`, owned by the maintainer, who drives this work. `master`
  (6536dad, tagged `v1.1.0-baseline`) must never be pushed to.
- Local clone: `/root/RandLib`. Its checked-out branch `1.x` is **stale** (6536dad plus old uncommitted
  copies of the work); do not use it and do not reset it without asking (resetting was refused once).
- Worktrees:

  | Path | Branch | Pushed as | Contents |
  |---|---|---|---|
  | `/root/RandLib-1.x` | `history` | `origin/1.x` | v1.2: C++17, original layout, fixes, CMake, tests, CI |
  | `/root/RandLib-2.x` | `2.x` | `origin/2.x` | v2: header-only C++20 (Phases 2–5 of the plan) |

- The agent's GitHub account (`nebius-llm-robot`) has pull access only. **The maintainer pushes**, from a
  VS Code terminal on the pod:

  ```sh
  cd /root/RandLib
  git push origin history:refs/heads/1.x 2.x
  git push origin v1.2.0            # only once 1.x CI is green, see 3.1
  ```

  Only add commits on top of what is pushed; never amend or rebase pushed commits. Commit messages end
  with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
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
/tmp/b2/tests/randlib_statistical_tests   # ~6 min
# 1.x: the same with -S /root/RandLib-1.x
```

Compilers on the pod: g++-13, g++-14, clang++-18. MSVC and AppleClang are only exercised by CI.

Test families (both branches): `tests/unit` (golden v1.1 comparison, mpmath reference values, coherence,
mode/median, parameters, regressions B1–B37, reproducibility, engines, numerics on 2.x) and
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
longer running; it did not finish the skewed stable laws with α = 1
(`StableRand(1, 0.5, 1, 0)`, `StableRand(1, 0.5, 2, 1)`, `StableRand(0.995, 0.5, 1, 0)`,
`StableRand(1, 0.005, 1, 0)`, `GeometricStableRand(1, 0.5, 1, 0)`). This is Phase 5 work: add a cheaper
mode (no tail quantiles below 1e-6) for α = 1. Once cached, write the CSVs with
`python3 tools/refdata/generate.py --cached-only --cache /root/RandLib/build/refdata-cache` from the
1.x worktree, remove the matching reviewed rules for those cases and regenerate the known deviations.

## 3. Immediate next actions

### 3.1 Finish v1.2 (branch `history` → `origin/1.x`)

1. Push the current `history` head. Commit c92f19a fixes B37 (`doubleFactorial`) and records a
   `2e-12` tolerance for the large-binomial CDF on MSVC; CI measured `1.20259e-12` relative error
   against the existing 40-digit reference. The earlier eafc3de warning fix built successfully on
   MSVC before that reference assertion failed.
2. When every job is green, tag: `git -C /root/RandLib-1.x tag -a v1.2.0 -m "RandLib 1.2.0"`, ask the
   maintainer to push the tag, update ROADMAP_STATUS (1.3, 1.5, 1.7) on both branches.

### 3.2 Phase 2 remainder (branch `2.x`)

- **CI portability**: local commits fix the Clang/AppleClang `-Wshadow` errors in `Pcg64`, discarded
  `[[nodiscard]]` results, MSVC's variadic-macro, secure-CRT, unsigned-negation, uninitialized-local
  and unused-variable diagnostics, and the coverage parser error. The engine-adapter test checks
  exact engine advancement while allowing eight ULPs of platform-dependent floating-point
  variation; AppleClang measured a maximum of six. Check these jobs after the next maintainer push.
- **2.5 Numerical building blocks**: done locally. `math/quadrature.hpp`, `math/roots.hpp`,
  `math/optimize.hpp`, and the functions in `math/functions.hpp`, `gamma.hpp`, `beta.hpp` and
  `special_functions.hpp` are templates on their floating-point type. Mixed integer/floating call
  sites preserve the floating type, and `tests/unit/test_numerics.cpp` independently instantiates
  and checks the `long double` path with exact identities and closed forms.
- **2.6 Tooling**:
  - Done locally: `.clang-format` matching the existing style, a pinned clang-format pre-commit hook,
    and one formatting pass over `include/` and `tests/`. Include sorting is deliberately disabled
    because the converted headers still contain dependency-sensitive implementation include blocks.
  - Done locally: clang-tidy analyzes `src/precompiled.cpp`, whose explicit instantiations cover the
    public distribution headers, and promotes every diagnostic to an error in blocking CI.
  - Done locally: opt-in `benchmarks/` with pinned Google Benchmark 1.9.5. Configure with
    `-DRANDLIB_BUILD_BENCHMARKS=ON`; it compares normal, exponential, gamma and Poisson sampling
    against `std::*_distribution` and measures stable sampling. Benchmarks are tracked, not gated.

    ```sh
    cmake -S /root/RandLib-2.x -B /tmp/randlib-bench -DCMAKE_BUILD_TYPE=Release \
      -DRANDLIB_BUILD_TESTS=OFF -DRANDLIB_BUILD_BENCHMARKS=ON
    cmake --build /tmp/randlib-bench --target randlib_benchmarks
    /tmp/randlib-bench/benchmarks/randlib_benchmarks
    ```
- **Clean-ups left by the automated header conversion** (worth doing with or before 2.6):
  - Done locally: `variate(detail::EngineRef)` and bulk `sample` overrides are private adapters at
    the start of each class; the redundant access-specifier toggles are removed.
  - Done locally: implementation helpers moved from the public `randlib::math` namespace into
    `randlib::math::detail`, including the special-function iteration constants, `erfinvAux*`,
    Lambert W, Marcum P/Q and portable Bessel helpers.
  - Done locally: `areClose` and the Newton tolerance are private math details; the public
    `Pair`/`Triplet` aliases moved from `math/constants.hpp` to `core/fwd.hpp`; the global pair
    arithmetic overloads were removed in favor of an explicit exponential-family inner product.
- After the current heads are pushed and every 2.x CI job is green, tag `v2.0.0-alpha.1` on `2.x`.
  Phase 2 is otherwise complete.

### 3.3 Phase 3: fitting and features (v2.0-beta)

Follow plan 3.1–3.3. Design notes gathered so far:
- Landed locally: `math::numericalHessian` and the first `include/randlib/fit/` vertical slice.
  `FitResult`, `FitOptions`, fixed-parameter constrained MLE, AIC/BIC and observed-information
  covariance are covered. The free API covers 33 distribution types, including fixed structural
  parameters, generic multi-parameter MLEs and bounded searches for integer-valued parameters.
  The free API now also has a named method-of-moments entry point, conjugate Bayesian updates,
  AIC/BIC model selection, and continuous/discrete goodness-of-fit tests. Conditional Bayesian
  overloads take the model's fixed parameter as their third argument (variance for normal,
  variance of the logarithms for log-normal, trial count for binomial, shape/count for
  negative-binomial, shape for Weibull, scale for Pareto, and minimum for uniform). Stable-family
  fitting remains with Phase 5. Seeded recovery grids cover representative closed-form, generic and
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
  stored as CSV (generate them with a script under `tools/`).
- 3.2 is done locally: `fit/tail.hpp` provides Hill and Pickands estimators, plot-ready Hill
  diagnostics and peaks-over-threshold fitting. Generalized Pareto and generalized extreme-value
  distributions have closed-form probability functions, sampling, moments, generic MLEs and
  committed 40-digit mpmath reference data. `fit/censored.hpp` supports exact, left-, right- and
  interval-censored observations plus one- and two-sided truncated MLE through the shared
  constrained-likelihood pipeline.
- 3.3 implementation is done: the unused private binomial `POISSON` selector is removed on both
  branches. The maintainer still needs to close #11 with the exact-sampling rationale.
- After these commits are pushed and every 2.x CI job is green, tag the resulting Phase 3 head
  `v2.0.0-beta.1`. Do not tag the current unpushed local head.

### 3.4 Phase 4: documentation, packaging, release (v2.0)

Landed locally: the README rewrite and comparison table; strict Doxygen Awesome generation and Pages
deployment; fitting, streams, stable-law, accuracy, packaging and migration guides; and the archived
legacy LaTeX source with its generated PDF removed. Four warning-clean examples (standard engine,
parallel streams, model selection and Bayesian update) build and run through CTest. The stable-fitting
example is intentionally coupled to the Phase 5 estimator rather than documenting an API that does
not exist yet.

Packaging is also landed and locally validated. `packaging/vcpkg/randlib` installs the current checkout
as an overlay and its installed `RandLib::RandLib` target built the consumer example with vcpkg
2026-09-26. `conan create . --build=missing` passed with Conan 2.33.0, packaging all 95 headers and
running `test_package`. `tools/amalgamate.py` produces a deterministic `randlib_single.hpp`; the
941,344-byte result compiled and ran with g++-13, g++-14 and clang++-18 under warnings-as-errors.
The release workflow compiles that artifact, extracts the exact tag's changelog section, and attaches
the header and SHA-256 checksum. The upstream vcpkg port must replace the checkout source with the
final release archive and checksum.

Remaining Phase 4 actions need the maintainer: enable GitHub Pages for Actions if it is not already
enabled, push and observe CI, run the release-candidate soak, submit the upstream vcpkg PR, and make
the release announcements. The final stable-fitting example will land as part of Phase 5.1.

### 3.5 Phase 5: hard numerics (v2.1), deliberately last

See plan Phase 5. Known concrete problems: geometric stable with α = 1 and β ≠ 0 takes ~1 min per
density value in the far tails (the coherence test checks only 0.01 < F < 0.99 for it); stable laws
near α = 1 throw `not_implemented`; heavy-tailed CFs accurate to ~1e-3 and the log-normal CF to 6e-5
(overrides in `tests/data/reference_tolerances.csv`); noncentral χ² upper tail via 1 − Marcum P.

Landed locally for 5.1: `StableRand::Parameterization::{S0,S1}` with S1 as the source-compatible
default. The constructor, getters and setters preserve the selected location convention; tests prove
S0/S1 equivalence for the pdf, CDF, survival function, characteristic function and seeded sampler,
including the α = 1 logarithmic conversion. The near-α=1 numerical kernel and stable fitting remain.
The first 5.2 primitive is also landed: `math::integrateTanhSinh` handles finite-interval endpoint
singularities, including an endpoint-distance callable form that avoids cancellation. Oscillatory
quadrature and the heavy-tail characteristic-function migrations remain.

## 4. Things that need the maintainer

Plan 0.1 and 0.5 (tracking issue, milestones), closing issues (#11 and others as they are fixed),
pushing branches and tags, release announcements, and the upstream vcpkg PR.
