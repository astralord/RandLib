# RandLib Modernization Plan

**Goal:** turn RandLib into a modern, header-only C++20 library that people choose for two things no mainstream C++ library offers:
**fitting distributions to data** and **heavy-tailed / stable distributions**. It should work alongside the standard `<random>` header and Boost.Math, not compete with them.

**Scope:** option 2 from the review (focused revival) plus everything from option 1 (minimal upkeep: bug fixes, build, CI, tests).

**Starting point:** `astralord/RandLib` at `6536dad` (June 2023): 146 source files, ~23.5k lines, qmake build, no tests, Travis CI (dead).

Progress against this plan is tracked in [ROADMAP_STATUS.md](ROADMAP_STATUS.md).

---

## Contents

1. [Guiding decisions](#1-guiding-decisions)
2. [Phases and milestones](#2-phases-and-milestones)
3. [Phase 0: Preparation](#phase-0-preparation)
4. [Phase 1: Stabilize (release v1.2)](#phase-1-stabilize-release-v12)
5. [Phase 2: Modern architecture (v2.0-alpha)](#phase-2-modern-architecture-v20-alpha)
6. [Phase 3: Fitting and features (v2.0-beta)](#phase-3-fitting-and-features-v20-beta)
7. [Phase 4: Documentation, packaging, release (v2.0)](#phase-4-documentation-packaging-release-v20)
8. [Phase 5: Hard numerics (v2.1)](#phase-5-hard-numerics-v21)
9. [Phase 6: Stable-count family and stable priors (v2.2)](#phase-6-stable-count-family-and-stable-priors-v22)
10. [Section 8: Open issues mapping](#8-open-issues-mapping)
11. [Section 9: Accuracy and testing policy](#9-accuracy-and-testing-policy)
12. [Section 10: Risks and open decisions](#10-risks-and-open-decisions)
13. [Appendix A: Baseline measurements](#appendix-a-baseline-measurements)
14. [Appendix B: Distribution tiers](#appendix-b-distribution-tiers)
15. [Revision history](#revision-history)

---

## 1. Guiding decisions

| Topic | Decision | Why |
|---|---|---|
| Language standard | **C++20** for v2 (v1.2 stays C++17) | Concepts (`std::uniform_random_bit_generator`), `std::span`, `std::numbers`. All major compilers have supported these since 2021. |
| Packaging | **Header-only** for v2, with an optional precompiled translation unit for `double` | Removes 282 explicit instantiations and the export macros, and lets users use any floating-point type |
| Build system | **CMake ≥ 3.21** only; remove qmake | Industry standard; makes vcpkg and Conan packaging possible |
| Random engines | Accept **any standard URBG**; ship our own engines as a convenience | Interoperates with `<random>`, `std::mt19937`, PCG and other users' code |
| Dependencies | **None** at runtime. Test-only: Catch2 (fetched by CMake). Tooling-only: Python + mpmath for reference data. | Keeps the "no dependencies" selling point |
| Naming | **Keep** existing class and method names (`NormalRand`, `Mean()`, `Fit()`), but inside `namespace randlib` | The value of v2 is semantics, not renaming; this keeps migration cheap |
| Scope control | **No distribution is removed** in v2.0, but distributions are tiered (Appendix B). Tier 1 gets full accuracy guarantees and fitting; Tier 2 is "tested, best effort". | Keeps the plan finishable |
| Versioning | SemVer. **v1.2** = non-breaking fixes. **v2.0** = breaking redesign. | Existing users get fixes without being forced to migrate |
| Branching | `1.x` branch for the stabilization release; v2 is developed on `2.x`, and `master` stays untouched until v2.0 is released | Isolates the breaking work |
| Ordering | **Engineering and refactoring first; hard numerics last** (Phase 5) | The architecture, API, fitting and packaging work is predictable; the numerical research is not, and must not block it |

---

## 2. Phases and milestones

| Phase | Outcome | Release | Rough effort* |
|---|---|---|---|
| 0. Preparation | Baseline frozen, golden data recorded, work tracked | — | 1–2 days |
| 1. Stabilize | Known bugs fixed, CMake, CI on 3 compilers, test suite | **v1.2** | 8–10 days |
| 2. Modern architecture | Header-only C++20, URBG interop, no hidden state, API cleanup, numerical building blocks | **v2.0-alpha** | 10–15 days |
| 3. Fitting and features | Unified fitting API, goodness of fit, tail estimators, extreme-value laws | **v2.0-beta** | 10–14 days |
| 4. Docs and packaging | Docs site, examples, migration guide, vcpkg/Conan, benchmarks | **v2.0** | 4–6 days |
| 5. Hard numerics | Stable laws near α = 1, heavy-tailed characteristic functions, Marcum Q, noncentral t, log-stable | **v2.1** | 15–25 days |
| 6. Stable-count family | Stable Count distribution, validated sampling, and a rigorous stable-prior investigation | **v2.2** | 5–8 days |

\* One person's focused working days. Phases 5 and 6 contain research-grade work and their uncertain
items are time-boxed.

Dependencies: 0 → 1 → 2 → 3 → 4 → 5 → 6. Phase 5 starts only when the engineering work of
Phases 2–4 is done. Phase 6 follows the stable-law kernel completed in Phase 5. Until then, numerical
problems that need research are recorded (as an issue and, where a test must be relaxed, in
`tests/data/reference_tolerances.csv`) and left for the corresponding research phase.

---

## Phase 0: Preparation

### 0.1 Coordinate upstream
- If this is not done by the maintainer: open a tracking issue on `astralord/RandLib` that links this plan. The maintainer said in #18 that he would review PRs.
- Agree on the `1.x` and `main` branch model and on the v2 breaking-change policy before writing code.

### 0.2 Tag the baseline
- Tag `v1.1.0-baseline` at `6536dad` so everything can be compared against it.

### 0.3 Record golden outputs (regression net for the refactor)
- Write `tools/golden/dump_v1.cpp`, linked against the *current* library. For every distribution, over a parameter grid × x grid, it dumps:
  `f`, `logf`, `F`, `S`, `Quantile(p)`, `Mean`, `Variance`, `Skewness`, `ExcessKurtosis`, `Median`, `Mode`, and `CF(t)` where available.
- Store the output as CSV under `tests/data/golden_v1/`.
- After every refactoring step, v2 must reproduce these values within tolerance. Values known to be wrong (Appendix A) are flagged in a `known_deviations.csv`, so tests expect the *fixed* value there instead.

### 0.4 Generate high-precision reference data
- Write `tools/refdata/` Python scripts using **mpmath** at ≥ 30 digits. These produce reference pdf, cdf and quantile values for Tier 1 distributions, stored as CSV in `tests/data/reference/`.
- **Lesson from the review:** do *not* use SciPy as the reference for stable distributions near α=1, because it snaps to the Cauchy result there. Use direct mpmath integration of the characteristic function for the pdf, and integration of the pdf (not Gil–Pelaez with `quadosc`) for the cdf when α<1.
- These scripts are run by developers, not in CI. Their CSV output is committed.

### 0.5 Project tracking
- One GitHub milestone per phase and one issue per numbered task in this plan, labeled `phase-N`, `tier-1`, `breaking`, `numerics`.

**Done when:** the baseline is tagged, golden and reference CSVs are committed, and the milestones exist.

---

## Phase 1: Stabilize (release v1.2)

Non-breaking. Stays C++17 and keeps the current file layout.

### 1.1 Repository hygiene (0.5 day)
- [ ] Remove `RandLib.pro.user`. It contains the author's local paths, and `.gitignore` misses it because it spells the name `Randlib.pro.user`.
- [ ] Remove LaTeX build files: `documentation/document.{aux,toc,synctex.gz}`. Add `*.aux *.toc *.log *.synctex.gz build/` to `.gitignore`.
- [ ] Remove `unit_tests/*` from `.gitignore`, since tests will live in the repo.
- [ ] Move `NoncentralTRand.{h,cpp}` to `attic/`. It uses the pre-template API and doesn't compile; it is ported properly in 3.4.
- [ ] Remove the dead `.travis.yml`.

### 1.2 CMake build (1 day)
- [ ] Add a `CMakeLists.txt` with target `RandLib::RandLib` (a static or shared library in v1.2), `install()` and `export()` rules, and a `RandLibConfig.cmake` so `find_package(RandLib)` works.
- [ ] Add a `CMakePresets.json` with `dev`, `release`, `asan` and `ubsan` presets.
- [ ] Make warnings opt-in: `RANDLIB_WARNINGS_AS_ERRORS=OFF` by default and `ON` in CI. Today `-Werror -Weffc++` in `RandLib.pro` breaks the build on GCC 13.
- [ ] Keep `RandLib.pro` for v1.2 only, marked deprecated, and remove it in v2.

### 1.3 Compiler portability (1 day). Closes #18
- [ ] Replace `constexpr` initializers that call non-constexpr `std::` math (verified fix: 7 lines, after which Clang 18 builds all files):
  - `NormalRand.h:56` and `ExponentialRand.h:53`: make the ziggurat tables `static inline const`, and drop `constexpr` from `createZiggurat()`.
  - `GumbelRand.cpp:103` (`M_LN_LN2`) and `NakagamiRand.cpp:446,459` (`medianCoef`, `skewness`): use `static const`.
- [ ] Remove the non-portable includes `<unistd.h>`, `<pthread.h>` and `<sys/types.h>` from `BasicRandGenerator.cpp`, which only needs `<thread>` and `<chrono>`.
- [ ] Fix the `-Wrange-loop-construct` errors in `BetaPrimeRand.cpp:193,210` and `LogisticRand.cpp`. Use `const RealType&` instead of `const double&`.
- [ ] Fix the `-Winconsistent-missing-override` warnings flagged by Clang (e.g. `PoissonRand::Sample`).
- [ ] Rename macros that collide with POSIX: `M_PI`, `M_E` and others in `math/Constants.h` are redefined as `constexpr long double`. Rename them to `randlib::consts::pi` and so on, keeping the old names behind `#ifndef` guards for v1.2.

### 1.4 Correctness fixes (3 days)

Each fix gets a regression test (see 1.6).

| ID | Bug | Fix |
|---|---|---|
| B1 | **Same seed doesn't reproduce normal samples.** The ziggurat's leftover exponential is a `static thread_local RealType z` (`NormalRand.cpp:125`), shared by every `NormalRand` in the thread. In the review, 995,504 of 1e6 values differed. | Remove the cache, which is a micro-optimization that only matters on the rare base strip. Alternatively, store it in the instance and clear it in `Reseed()`. |
| B2 | **JKiss seeding is weak.** `Reseed()` folds the seed into only `X`, truncated to 32 bits; `Y`, `Z` and `C` are constant. Seeds 7 and 7+2³² give identical streams. | Expand the seed with SplitMix64 into all state words, and reject the all-zero xorshift state. |
| B3 | **PCG seeding doesn't follow the reference.** It sets `state = inc = seed`. | Implement `pcg32_srandom_r`: `state=0; inc=(seq<<1)|1; step; state+=seed; step`. |
| B4 | **Default seed is time-based** (`time(0)` ⊕ thread id ⊕ counter). Two processes collide when address randomization is off. | Use `std::random_device` mixed with a high-resolution clock and a per-process counter, all through SplitMix64. |
| B5 | **The 64-bit generator switch is broken.** The header checks `JLKISS64RAND` while `UniformRand.cpp` checks `RANDLIB_JLKISS64`, and both alternative branches are invalid C++. With the 64-bit engine selected, "uniform" values reach about 4×10⁹. | Write one `uniform01<Real>(gen)` that derives the bit count from `gen.MaxValue()`, and delete the broken `#ifdef` branches. |
| B6 | **Uniforms have only 32 bits of resolution.** `(x+0.5)/2³²` truncates inverse-transform tails; for example `-log(U)` can't exceed ≈22.9. | Produce 53-bit doubles from two 32-bit draws, or one 64-bit draw. For `long double`, use 64 bits. Keep open-interval semantics (0,1). |
| B7 | **The stable constructor rejects valid parameters.** 0.99<α<1.01 with β≠0, and α=1 with 0<\|β\|<0.01, throw an exception, which blocks *sampling* too, even though sampling works there. | Allow construction and sampling. Throw only from `f`, `F` and `Quantile` in the unsupported region, until 3.3 lifts that. |
| B8 | **The stable CDF has a floor bug.** For α<1 and β=1, `F(x)` sticks at 3.8×10⁻⁵ where the density is 0. | Investigate the tail-integration branch for totally skewed α<1 and add tests at the support boundary. |
| B9 | `(signed)B > 0` in the ziggurat sign selection biases toward negative values when the low 32 bits are 0. This effect is negligible. | Use an explicit sign bit: `(B >> 8) & 1`. |
| B10 | `Reseed()` is a `const` method that mutates a `mutable` generator. | Leave it in v1.2 for compatibility and fix it in v2 (2.4). |

### 1.5 Continuous integration (1 day)
- [ ] GitHub Actions matrix:
  - Linux: GCC 11, 12, 13, 14 and Clang 15, 16, 17, 18
  - macOS: AppleClang
  - Windows: MSVC 2022
- [ ] Extra jobs: ASan+UBSan (Clang), coverage (gcov→Codecov, informational), and `clang-tidy` (informational at first).
- [ ] Badges in the README replace the Travis badge.

### 1.6 Test suite (3–4 days). Closes #12

Framework: **Catch2 v3** via `FetchContent`. Layout: `tests/unit/`, `tests/statistical/`, `tests/data/`.

| Test family | What it asserts | Applies to |
|---|---|---|
| Reproducibility | Same seed gives the same sequence, independent of other instances and of call history (B1) | All |
| Seeding | Distinct seeds give distinct streams; seeds differing above bit 32 differ (B2) | Engines |
| Goodness of fit | Fixed seed, n=2×10⁵ samples; KS (continuous) or χ² (discrete) p-value > 10⁻⁴ | All Tier 1, sampled Tier 2 |
| Moments | Sample mean and variance within 5 standard errors of `Mean()` and `Variance()` | Where finite |
| Reference values | `f`, `F` and `S` against `tests/data/reference` within the Tier tolerance (Section 9) | Tier 1 |
| Golden regression | Matches `golden_v1` except `known_deviations` | All |
| Internal coherence | `f = exp(logf)`; `F + S = 1`; `F` monotone; `Quantile(F(x)) ≈ x`; `∫f ≈ F(b)−F(a)` | All |
| Mode and median | `f(mode±ε) ≤ f(mode)`; `F(median) = 0.5` | All |
| Fitting | On data sampled from known parameters, MLE recovers them within CI; log-likelihood is maximal: `L(θ̂±ε) ≤ L(θ̂)` | All with `Fit*` |
| Parameter validation | Invalid parameters throw `std::invalid_argument` | All |

The three checklist items in #12 (MLE optimality, mode, pdf/logpdf coherence) are all covered.

**Determinism:** statistical tests use fixed seeds, so they are deterministic and not flaky. If a test fails after a sampler change, the fix is to investigate, never to reseed until it passes.

### 1.7 Release v1.2 (0.5 day)
- [ ] Write a `CHANGELOG.md` entry. Call out B1–B6 explicitly, since **sampled values change for the same seed**.
- [ ] Tag `v1.2.0` on the `1.x` branch.

**Phase 1 done when:** CI is green on all compilers, the test suite covers every distribution, B1–B9 are fixed, and v1.2.0 is tagged.

---

## Phase 2: Modern architecture (v2.0-alpha)

Breaking changes from here on. C++20.

### 2.1 Header-only layout and namespace (3 days)

```
include/randlib/
  randlib.hpp                 umbrella header (everything)
  engine/   splitmix64.hpp  pcg.hpp  xoshiro.hpp  seeding.hpp  uniform.hpp
  core/     concepts.hpp  distribution.hpp  univariate.hpp  continuous.hpp  discrete.hpp  errors.hpp
  math/     constants.hpp  gamma.hpp  beta.hpp  bessel.hpp  marcum.hpp  lambertw.hpp
            quadrature.hpp  roots.hpp  optimize.hpp
  continuous/ normal.hpp gamma.hpp beta.hpp stable.hpp ...   (one per distribution)
  discrete/   poisson.hpp binomial.hpp ...
  circular/   von_mises.hpp wrapped_exponential.hpp
  singular/   cantor.hpp
  bivariate/  bivariate_normal.hpp normal_inverse_gamma.hpp trinomial.hpp
  fit/        result.hpp mle.hpp moments.hpp bayes.hpp gof.hpp select.hpp
src/precompiled_double.cpp    optional: extern template instantiations for double
tests/  examples/  benchmarks/  tools/  docs/
```

- Everything lives in `namespace randlib`, with internals in `randlib::detail`.
- Remove `RandLib_global.h` and `RANDLIBSHARED_EXPORT`.
- **Compile-time mitigation:** each distribution header includes only what it needs. If the option `RANDLIB_PRECOMPILED=ON` is set, the library declares `extern template` for `double` and builds `src/precompiled_double.cpp` once.
- **Constants:** use `std::numbers` where available. Tables that need math functions (ziggurat, factorials) are generated by `tools/gen_tables.py` as literal arrays, so they are truly `constexpr` on every compiler. This makes the 1.3 workaround permanent.

### 2.2 Random engine interoperability (3 days)

The core of the redesign:

```cpp
namespace randlib {
template <class D>
concept SamplableDistribution = requires(const D& d, std::mt19937_64& g) {
    { d(g) } -> std::convertible_to<typename D::value_type>;
};

template <std::floating_point Real = double>
class NormalRand : public ContinuousDistribution<Real> {
public:
    using value_type = Real;
    template <std::uniform_random_bit_generator G>
    Real operator()(G& g) const;                          // standard-style sampling
    template <std::uniform_random_bit_generator G>
    void Sample(G& g, std::span<Real> out) const;         // bulk sampling
    Real Variate() const;                                 // convenience: per-instance engine
    void Sample(std::span<Real> out) const;
    // ...
};
}
```

- Every sampler is written once against a generic URBG and works with `std::mt19937`, `std::minstd_rand`, PCG or user engines.
- `detail::uniform01<Real>(G&)` handles any URBG range (32-bit, 64-bit and non-power-of-two `min()`/`max()`).
- **Bundled engines** (all satisfy URBG):
  - `Pcg32`, `Pcg64` (reference seeding)
  - `Xoshiro256pp`, with `jump()` and `long_jump()` for **independent parallel streams** (one jump per thread)
  - `SplitMix64` (seeding only)
  - `DefaultEngine = Xoshiro256pp`. JKiss and JLKiss64 are dropped, or kept under `randlib::legacy`.
- **Seeding:** `randlib::seed_from(std::seed_seq&)`, `randlib::random_seed()` (from `std::random_device`), and `randlib::make_streams<Engine>(seed, n)`, which returns n non-overlapping engines.
- Convenience `Variate()` keeps working with a per-instance `DefaultEngine` seeded from `random_seed()`. **Documented as not thread-safe per instance**; use `operator()(G&)` with one engine per thread.

### 2.3 Remove hidden state (1 day)
- [ ] No `static` or `thread_local` mutable state anywhere. Enforce this with a CI grep check.
- [ ] Every sampler's output is a pure function of its parameters and the engine state.

### 2.4 API cleanup (2–3 days)
- [ ] **Unambiguous constructors.** Keep `NormalRand(mean, variance)` for continuity, but add named factories: `NormalRand::FromMeanStd(μ, σ)` and `NormalRand::FromMeanVariance(μ, σ²)`. Apply the same pattern wherever a parameter is a variance, scale or rate (Gamma shape/rate vs shape/scale, LogNormal, etc.).
- [ ] `Reseed()` and the setters become non-`const`.
- [ ] Vector inputs take `std::span<const T>` instead of `const std::vector<T>&`.
- [ ] Getters: `[[nodiscard]]` everywhere; `noexcept` where guaranteed.
- [ ] **Errors:** keep exceptions for invalid parameters (`std::invalid_argument`). Fitting failures use `randlib::fit_error` with an enum code, replacing the string-built `fitErrorDescription`.
- [ ] Remove the `String` typedef; use `std::string`.
- [ ] Replace the `FIT_ERROR_TYPE` and `SUPPORT_TYPE` C-style enums with `enum class`.
- [ ] Keep virtual dispatch for the runtime-polymorphic interface (`UnivariateDistribution<T>*`), but mark leaf classes `final`.

### 2.5 Numerical building blocks (2–3 days)
- [ ] `math/quadrature.hpp`: the adaptive Gauss–Kronrod integrator behind one policy-style entry point
  (`integrate(f, a, b, tol)`), used by every distribution; compare against the golden data and note any
  accuracy changes.
- [ ] `math/roots.hpp`: TOMS 748 (bracketing) and Halley/Newton (with derivatives). Generic quantile
  inversion uses these.
- [ ] `math/optimize.hpp`: Brent (1-D) and Nelder–Mead / BFGS with numeric gradients (n-D) for generic MLE
  in 3.1.
- [ ] All math functions are templates on `Real`. The double path keeps today's precision; the
  `long double` path is tested separately.
- Specialized quadrature for hard integrands (tanh-sinh, oscillatory) moves to Phase 5 (5.2).

### 2.6 Tooling (1 day)
- [ ] `clang-format` config and pre-commit hook; `clang-tidy` set to fail CI on new warnings.
- [ ] Benchmarks with Google Benchmark, in `benchmarks/`, comparing against `std::normal_distribution`, `std::gamma_distribution`, etc., and against Boost.Random if available. These are tracked, not gated.

**Phase 2 done when:** the library is header-only and C++20; every test from Phase 1 is ported and passing; golden data is reproduced; samplers accept any URBG; there is no static mutable state; and v2.0.0-alpha.1 is tagged.

---

## Phase 3: Fitting and features (v2.0-beta)

### 3.1 Fitting as a first-class feature (6–8 days)

Today 17 distributions have about 20 different ad-hoc `Fit*` methods that mutate the object (`FitScale`, `FitShapeBayes`, `FitProbabilityMinimax`, …). v2 unifies them into one API.

```cpp
std::vector<double> data = /* ... */;

// One distribution, chosen method
auto r = randlib::fit<randlib::GammaRand<double>>(data, randlib::Method::MLE);
r.distribution;        // fitted GammaRand<double>
r.log_likelihood;      // double
r.std_errors;          // from observed Fisher information
r.covariance;          // parameter covariance matrix
r.converged, r.iterations;

// Fix some parameters, fit the others
auto r2 = randlib::fit<randlib::GammaRand<double>>(data, {.fixed = {{"rate", 2.0}}});

// Bayesian (conjugate) returns the posterior
auto post = randlib::fit_bayes(data, randlib::NormalInverseGammaRand<double>(0, 1, 1, 1));

// Model selection
auto best = randlib::fit_best<randlib::GammaRand<double>,
                              randlib::WeibullRand<double>,
                              randlib::LogNormalRand<double>>(data, randlib::Criterion::AIC);

// Goodness of fit
auto gof = randlib::goodness_of_fit(r.distribution, data);   // KS, Anderson–Darling, CvM + p-values
```

Tasks:
- [ ] `FitResult<D>` struct: fitted distribution, log-likelihood, AIC/BIC, standard errors, covariance, convergence information, and the method used.
- [ ] Methods: `MLE`, `MethodOfMoments`, `UMVU` (where it exists), `Bayes` (conjugate pairs), and `Quantile`-based (robust initial guess; McCulloch for stable, see 5.1).
- [ ] **Generic MLE fallback:** any distribution exposing `LogLikelihood(params, data)` and parameter bounds is fittable via `optimize.hpp`. Closed-form or specialized solvers override it where they exist (today's 17).
- [ ] Extend fitting to **all Tier 1 distributions** (Appendix B). Closed form where possible, generic MLE otherwise.
- [ ] Standard errors via the observed Fisher information (numeric Hessian of the log-likelihood).
- [ ] `goodness_of_fit`: Kolmogorov–Smirnov (exact distribution via the existing `KolmogorovSmirnovRand`), Anderson–Darling, Cramér–von Mises; χ² for discrete.
- [ ] Keep the old member methods (`X.Fit(data)`) as thin deprecated wrappers until v3.
- [ ] Tests: parameter recovery over a grid of true parameters and sample sizes; CI coverage ≈ nominal (simulation test, nightly only); comparison to SciPy `fit` results on fixed datasets.

### 3.2 Tail estimators and extreme-value distributions (3 days)
- [ ] Tail-index estimators: Hill and Pickands, with a plot-ready Hill function output.
- [ ] Generalized Pareto and Generalized Extreme Value distributions (peaks-over-threshold and block-maxima fitting). These are a natural fit with the existing Pareto, Gumbel and Fréchet.
- [ ] Truncated and censored data support in MLE (likelihood with `F`/`S` terms).

### 3.3 Close #11
- [ ] Close the issue with an explanation: the exact samplers (Bernoulli sum, waiting time, BTRD-style rejection) are already fast for small `np`, and replacing them with a Poisson approximation would add bias for little gain. Remove the unused `POISSON` value from `GENERATOR_ID`.
- [ ] Optional: benchmark Binomial sampling for small `np` and switch to inversion if it's faster.

**Phase 3 done when:** the fitting API covers all Tier 1 distributions (for stable laws, with the accuracy
available at the time); the tail estimators and extreme-value laws are tested against mpmath; and
v2.0.0-beta.1 is tagged.

---

## Phase 4: Documentation, packaging, release (v2.0)

### 4.1 Documentation (2–3 days)
- [ ] **README rewrite:** positioning ("fit distributions and handle heavy tails in C++"), quick start, a comparison table vs `<random>`, Boost.Math and StatsLib (honest about where those are better), and badges.
- [ ] **API reference:** Doxygen with the Doxygen Awesome theme, published to GitHub Pages by CI.
- [ ] **Guides:** fitting, parallel sampling with streams, stable distributions (parameterizations, accuracy table), and accuracy policy.
- [ ] Retire `documentation/document.pdf`. Keep its LaTeX source in `docs/legacy/` as a mathematical reference, and link to it.
- [ ] `examples/`: compiled in CI. At minimum: sampling with `std::mt19937`, parallel streams, fit + model selection, Bayesian update, stable fitting.

### 4.2 Migration guide v1 → v2 (0.5 day)
A table of every breaking change with before/after code: namespace, include paths, `Reseed` const-ness, `Fit` → `fit<>()`, engine changes, and the reproducibility break.

### 4.3 Packaging (1 day)
- [ ] vcpkg port (overlay first, then upstream PR) and a Conan recipe.
- [ ] Single-header amalgamation `randlib_single.hpp`, generated by a script and attached to each GitHub release.
- [ ] Release workflow: tag → GitHub Release with changelog excerpt and amalgamated header.

### 4.4 Release v2.0.0
- [ ] Release candidates for two weeks, then final. Open Phase 5 items are listed in the release notes
  with their current behavior (for example: stable pdf/cdf near α = 1 throw `not_implemented`).
- [ ] Announce on the GitHub repo, r/cpp and the issues of people who asked for features (#18's reporter).

### 4.5 Mathematical distribution reference
- [x] Document every registered concrete distribution name with TeX definitions of its density or
  mass, CDF and direct tail where applicable, generating or characteristic functions, moments,
  support and important special cases.
- [x] Explain the numerical algorithms used for probability functions, quantiles, sampling and
  fitting, including direct-tail and limiting branches that differ from the simplest formula.
- [x] Preserve the useful content of the archived `document.tex` while checking parameterizations
  and implementation claims against the C++20 headers and regression evidence.
- [x] Make the strict documentation target compare the executable distribution registry with the
  mathematical-reference headings, so a newly registered class cannot silently remain undocumented.
- [ ] Add generated plots in a separate visualization phase. Formula and algorithm coverage does not
  depend on figures.

---

## Phase 5: Hard numerics (v2.1)

Research-grade numerical work, deliberately last. Each item is time-boxed; what misses its target keeps
its documented v2.0 behavior and stays in the tracker.

### 5.1 Stable distributions complete (time-boxed: 8 days). Closes #1
Today the pdf and cdf agree with high-precision references to ~10⁻¹⁰ wherever they are allowed. The work is the excluded region and the tails.
- [ ] Add the **S0 parameterization** (Nolan) alongside S1 (the current one), with `StableRand::Parameterization::{S0, S1}`. S0 is continuous in α at α=1, which makes the near-1 region tractable.
- [ ] Pdf and cdf for **α near 1 with β≠0** and for **α=1 with small β**: implement Nolan's integral representation in S0 using the tanh-sinh quadrature from 5.2. Where the integrand is ill-conditioned, switch to series expansions (small or large |x|) or interpolate in α across the gap, with an error bound.
- [ ] Quantile function via TOMS 748 on `F` with good initial guesses.
- [ ] Fitting: McCulloch quantile estimator (fast initial guess), then MLE (3.1). Validate on simulated data.
- [ ] Accuracy target: ≤10⁻⁹ absolute pdf and cdf over α∈[0.1,2], β∈[−1,1], |x|≤10⁶, checked against mpmath reference data.
- [ ] Geometric stable laws with α = 1 and β ≠ 0: their density takes about a minute in the far tails, because
  the mixture needs the stable law at huge arguments (asymptotic expansions of the α = 1 tails), and their
  1.2 reference values are still incomplete.
- [ ] **Time-box rule:** if the near-α=1 region misses the target after 8 days, keep the B7 behavior (sampling
  works; pdf/cdf throw with a clear message) and record what remains.

### 5.2 Characteristic functions of heavy-tailed laws (3–4 days). Closes #6 and #10
Measured in v1.2: the characteristic functions of Pareto, Fréchet, inverse gamma, Fisher F and beta prime
are accurate to about 10⁻³, and the log-normal one to 6×10⁻⁵ at σ²=4 (`tests/data/reference_tolerances.csv`).
- [ ] `math/quadrature.hpp`: **tanh-sinh** for endpoint singularities (stable and log-stable need it) and
  **oscillatory integration** (QAWF-style) for characteristic functions.
- [ ] Re-run the heavy-tailed characteristic functions with the oscillatory quadrature. Change of variables y=log x and splitting at zeros of cos/sin should give most of the gain.
- [ ] Only if the target (≤10⁻⁸) is still missed: implement the Laplace-transform / Lambert-W approach referenced in the issue (needs complex Lambert W, which goes in `math/lambertw.hpp`).

### 5.3 Marcum Q / noncentral χ² accuracy (2–3 days). Closes #4
Measured today: noncentral χ² CDF relative error of 2×10⁻⁴ in the bulk and up to 5×10⁻³ in the lower tail.
- [ ] Implement Gil, Segura & Temme (2014), "Computation of the Marcum Q-function" (ACM TOMS), which the author already cites in #4: series, asymptotic expansions and a quadrature region selected by (μ, x, y).
- [ ] Use it for `NoncentralChiSquaredRand::F`/`S`, the Rice distribution (new, Tier 2, almost free once Marcum Q exists), and the noncentral F if added later.
- [ ] Target ≤10⁻¹² relative error in the bulk and ≤10⁻⁹ in the tails.

### 5.4 Noncentral t (3 days). Closes #7
- [ ] Port `attic/NoncentralTRand` to the v2 templated design.
- [ ] CDF: Lenth's AS 243 / Benton & Krishnamoorthy (2003) series, using the existing regularized incomplete beta. Pdf via the recurrence between the pdf and the cdf at different degrees of freedom.
- [ ] Extreme tails with large noncentrality: implement the integral approach from the paper cited in #7 if time allows, otherwise document the limits.
- [ ] Moments, quantile, sampling (Z+δ)/√(χ²/ν), and tests against mpmath.

### 5.5 Log-stable distribution (2 days, after 5.1). Closes #5
- [ ] `LogStableRand` as the exponential of a stable variable: sampling, pdf, cdf, quantile, median and mode by transformation. Only allowed for β=−1 when moments are needed. Document that for other β no moments exist.
- [ ] **No characteristic function** in v2.0. `CF()` throws `not_implemented` with an explanation. Revisit on demand.

**Phase 5 done when:** #1, #4, #5, #6, #7 and #10 are closed or explicitly moved on under their time-box
rules; accuracy targets are met in CI; and v2.1.0 is tagged.

---

## Phase 6: Stable-count family and stable priors (v2.2)

This phase revisits the experimental [`StableCountRand` commit](https://github.com/astralord/RandLib/commit/c0700988a88767cde9967b41e0cc459b8434cc6a)
against Lihn's revised [Wright-function treatment](https://arxiv.org/abs/2405.04693). The old code is
design evidence, not an implementation to port: its sampler uses an arbitrary hitting-time cutoff,
rejects large stable draws and keeps mutable static state. The v2 implementation starts from the
published density and the Phase 5 stable kernel.

### 6.1 Stable Count distribution (3–5 days)

- [ ] Add `StableCountRand<RealType>` under `include/randlib/continuous/` with exponent
  `0 < alpha < 1`, positive scale and a finite location. Use the location-scale form of
  `N_alpha(x) = L_alpha(1/x) / (Gamma(1/alpha + 1) x)` and state exactly which one-sided stable
  parameterization and scale convention it requires.
- [ ] Implement direct log-density, density, CDF, survival, quantile and all finite moments. Evaluate
  tails directly instead of subtracting from one. Use the exact moment identity
  `E[X^q] = Gamma((q + 1)/alpha) / (Gamma(q + 1) Gamma(1/alpha))` in the standard case, with
  log-gamma arithmetic when the raw moment would overflow.
- [ ] Preserve the `alpha = 1/2` reduction to `Gamma(shape=3/2, scale=4)` as an exact fast path and
  document the degenerate limit as `alpha -> 1`. Do not retain the experimental `[0.1, 1)` or
  near-one exclusions unless new numerical evidence requires a documented bound.
- [ ] Supply a standard-URBG sampler with no hidden state. Inverse-CDF sampling is an acceptable
  correctness-first implementation; a transformed or rejection sampler must have a derivation and
  statistical evidence. Do not use the experimental fixed `100000` threshold or discard large
  one-sided-stable draws.
- [ ] Generate committed 40-digit mpmath reference grids from the Wright/one-sided-stable definition,
  cross-check the half-stable Gamma case and moment identities, and cover normalization, support,
  density/CDF/survival coherence, quantile round trips, reproducibility and seeded sampling tests.
- [ ] Add generic MLE with fixed-parameter variants and seeded recovery tests. Land as Tier 2 first;
  promote to Tier 1 only after the independent reference and sampling evidence meet the promotion
  rule.

### 6.2 Stable-prior semantics (time-boxed: 2–3 days)

- [ ] Formalize the claimed "conjugate prior" before adding an API: identify the random quantity
  (count, count root, scale or reciprocal scale), write the normalized likelihood, and derive the
  posterior hyperparameter update for an arbitrary sample size. Distinguish Bayesian conjugacy from
  the paper's inverse-variable, product-mixture and subordination identities.
- [ ] Add a `fit_bayes` overload for `StableRand` only if the posterior remains in the Stable Count
  family and the supported stable parameterization is explicit. Verify normalization, posterior
  moments and posterior-predictive values against high-precision quadrature, including the
  `alpha = 1/2` special case.
- [ ] If posterior closure cannot be proved within the time box, keep `StableCountRand` as a
  standalone distribution and document the subordination relation without calling it a conjugate
  Bayesian prior. Record the missing derivation as a research issue rather than exposing a
  misleading closed-form update.

The generalized Stable Count family `(alpha, sigma, d, p)` from the 2025 revision is deliberately a
follow-up, not part of the first implementation. It should be proposed separately after the core
Wright-function evaluation and sampler are independently validated.

**Phase 6 done when:** the core Stable Count distribution has reference-backed probability functions,
moments, fitting and a validated standard-URBG sampler; the Bayesian-conjugacy claim is either proved
and implemented or explicitly scoped to a non-Bayesian subordination identity; and v2.2.0 is tagged.

---

## 8. Open issues mapping

| # | Title | Where | Plan |
|---|---|---|---|
| 18 | constexpr not cross-compilable | 1.3, 2.1 | Quick fix in v1.2 (already verified with Clang 18); permanent generated tables in v2 |
| 12 | Unit-testing | 1.6 | Full suite in v1.2, extended in every phase |
| 10 | Integration techniques | 2.5, 5.2 | Gauss–Kronrod in v1.2; tanh-sinh and oscillatory quadrature in Phase 5 |
| 1 | Stable pdf/cdf | 1.4 (B7, B8), 5.1 | Unblock sampling in v1.2; S0 parameterization and near-α=1 work, time-boxed |
| 4 | Special functions (Marcum) | 5.3 | Gil–Segura–Temme algorithm |
| 6 | Log-normal CF | 5.2 | New quadrature first, Lambert-W only if needed |
| 7 | Noncentral t | 5.4 | Port and implement |
| 5 | Log-stable | 5.5 | Everything except CF |
| 11 | Binomial ≈ Poisson | 3.3 | Close as won't-fix with explanation |

---

## 9. Accuracy and testing policy

| | Tier 1 | Tier 2 |
|---|---|---|
| pdf/pmf, cdf, survival (bulk: 10⁻⁶ < F < 1−10⁻⁶) | ≤ 10⁻¹² relative | ≤ 10⁻⁸ relative |
| Tails (F or S down to 10⁻³⁰⁰) | ≤ 10⁻⁹ relative | best effort, documented |
| Quantile round trip | ≤ 10⁻¹⁰ relative | ≤ 10⁻⁷ |
| Moments (closed form) | exact to rounding | exact to rounding |
| Characteristic function | ≤ 10⁻⁸ absolute | best effort |
| Stable (special case) | ≤ 10⁻⁹ absolute pdf/cdf | n/a |
| Goodness-of-fit sampling test | required | required |
| Fitting | required (closed form or generic MLE) | optional |

- Reference data is always from mpmath at ≥30 digits, never from another double-precision library.
- Any loosened tolerance requires a comment in the test that links an issue.
- Survival functions are computed directly, not as `1 − F`, so tails stay accurate.

---

## 10. Risks and open decisions

| Risk / decision | Impact | Mitigation |
|---|---|---|
| **Stable near α=1 is research-grade** | Could stall the release | Moved to Phase 5 (5.1), after v2.0; time-boxed, with clear errors meanwhile |
| **Header-only compile time** with ~55 distributions | Slow builds for users who include everything | Per-distribution headers; `extern template` precompiled option; measure in CI |
| **Reproducibility break** (B1–B6 change sampled values) | Users with stored seeds get different numbers | Call it out in v1.2 release notes; offer `randlib::legacy::JKissEngine` for those who must reproduce old runs (with the old bugs) |
| **C++20 requirement** | Excludes very old toolchains | v1.2 remains C++17 and gets critical fixes for 12 months |
| **Upstream ownership** | Plan needs maintainer buy-in if done from a fork | Phase 0.1. The maintainer has said he will review PRs. Land Phase 1 as small PRs. |
| **Fitting scope creep** (censoring, mixtures, regression) | Phase 3 grows unbounded | v2.0 = i.i.d. univariate fitting + censoring (3.2) only; mixtures listed for v2.x |
| **"Conjugate prior" is underspecified for Stable Count** | A Bayesian API could imply posterior closure that the cited product/subordination identity does not establish | Phase 6.2 requires a normalized likelihood and posterior-update derivation before any `fit_bayes<StableRand>` overload |
| Python bindings (nanobind) | Wider audience, but SciPy covers most of it | Out of scope for v2.0; revisit if stable/fitting proves unique enough |

**Decisions confirmed before Phase 2** (all as recommended):
1. C++20 (recommended) vs staying on C++17 with `std::span` replaced by pointer+size.
2. Keep `*Rand` class names (recommended for v2) vs a snake_case rename (more idiomatic, but more churn).
3. Default engine: Xoshiro256++ (recommended: fast, with jump-ahead) vs PCG64.

---

## Appendix A: Baseline measurements

Measured on the baseline commit with GCC 13.3 and Clang 18. These are the numbers the test suite must improve on or preserve.

**Sampling quality:** passes KS (√n·D < 1.95 at n=2×10⁵) and moment checks for Normal, Gamma(0.3,2), Beta(.5,.5), StudentT(5), LogNormal(0,.25), Stable(1.5,.5), Stable(.7,−.3).

**Reproducibility:** two `NormalRand` instances seeded with 42 produce different sequences if another normal was drawn between them (995,504 of 10⁶ values differ, first at index 4,496).

**Seeding:** seeds 7 and 7+2³² give identical uniform streams. With address randomization disabled, two default-seeded processes in the same second produce identical values.

**Build:**
- GCC 13: one file (`NoncentralTRand.cpp`) fails; with the project's own `-Werror` flags, `BetaPrimeRand.cpp` and `LogisticRand.cpp` fail too.
- Clang 18: 54 files fail on `constexpr` initializers; a 7-line fix makes all of them build, with identical numerical results.

**Stable pdf/cdf (S1) vs high-precision reference:** ≤10⁻¹⁰ over α∈{0.5…1.9}, β∈{−0.7…1}, x∈[−5,10], except:
- the rejected region (B7)
- the β=1, α<1 CDF floor of 3.8×10⁻⁵ (B8)
- tail values of the cdf for α≈1±0.01, β≠0, where differences of ~10⁻³ vs SciPy were seen; this needs mpmath confirmation in Phase 0.4, because SciPy itself is unreliable there

**Noncentral χ² CDF vs exact series:** relative error 3×10⁻⁵ to 2×10⁻⁴ in the bulk; up to 5×10⁻³ in the lower tail (k=50, λ=0.01).

**Log-normal CF vs mpmath:** absolute error ≤3×10⁻¹¹ for σ²≤0.25; 2.5×10⁻⁶ (t=5) and 2.9×10⁻⁵ (t=20) at σ²=1; 2.5×10⁻⁴ to 4.1×10⁻³ at σ²=4. Each evaluation takes 100–600 µs.

---

## Appendix B: Distribution tiers

**Tier 1** (full accuracy targets, fitting required, ~30):
- **Continuous:** Normal, LogNormal, Exponential, Gamma (+ Chi-squared, Erlang), Beta, BetaPrime, StudentT, Cauchy, Lévy, **Stable**, **GeometricStable**, **LogStable** (new), Pareto, Weibull, Fréchet, Gumbel, Laplace, Logistic, InverseGamma, InverseGaussian, Uniform, FisherF, NoncentralChiSquared, **NoncentralT**
- **Discrete:** Poisson, Binomial, NegativeBinomial, Geometric, Bernoulli, UniformDiscrete
- **Bivariate:** NormalInverseGamma (as a conjugate prior)

**Tier 2** (tested for coherence and sampling, best-effort accuracy):
- **Continuous:** ExponentiallyModifiedGaussian, IrwinHall, KolmogorovSmirnov, MarchenkoPastur, Nakagami, Planck, RaisedCosine, Sech, Triangular, WignerSemicircle, Degenerate, Rice (new, via 3.5), StableCount (planned in Phase 6)
- **Discrete:** BetaBinomial, Categorical, HyperGeometric, NegativeHyperGeometric, Logarithmic, Rademacher, Skellam, Yule, Zeta, Zipf
- **Circular:** VonMises, WrappedExponential
- **Singular:** Cantor
- **Bivariate:** BivariateNormal, Trinomial

**Promotion rule:** a Tier 2 distribution becomes Tier 1 when someone needs it and contributes reference data and fitting.

---

## Revision history

- **2026-10-17:** added Phase 6 for a clean Stable Count implementation based on the published
  Wright-function definition and for a separately gated investigation of the claimed stable-law
  conjugate-prior semantics. The 2021 experimental sampler is retained as design evidence only.
- **2026-10-11:** at the maintainer's request, the hard numerical work on heavy-tailed laws (former 2.5
  tanh-sinh and oscillatory quadrature, 3.3–3.7) moved to the new Phase 5, after all engineering and
  refactoring work; v2.0 is released after Phase 4. The v2 branch is `2.x` instead of `main`.
