# RandLib

[![CI](https://github.com/astralord/RandLib/actions/workflows/ci.yml/badge.svg?branch=2.x)](https://github.com/astralord/RandLib/actions/workflows/ci.yml)
[![Documentation](https://github.com/astralord/RandLib/actions/workflows/docs.yml/badge.svg?branch=2.x)](https://astralord.github.io/RandLib/)
[![codecov](https://codecov.io/gh/astralord/RandLib/branch/2.x/graph/badge.svg)](https://codecov.io/gh/astralord/RandLib)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C.svg)](https://en.cppreference.com/w/cpp/20)

Fit probability distributions and handle heavy tails in modern C++. RandLib is a dependency-free,
header-only C++20 library for sampling, probability functions, moments, parameter estimation, model
selection, goodness-of-fit testing, extreme-value analysis, and stable laws.

The `2.x` branch is the next major version. The source-compatible C++17 maintenance release remains
on `1.x`; see the [v1 to v2 migration guide](docs/migration-v2.md).

## Quick start

```cpp
#include <randlib/randlib.hpp>

#include <iostream>
#include <random>
#include <vector>

int main()
{
    std::mt19937_64 engine(42);
    const randlib::GammaRand<double> model(2.5, 0.5);
    std::vector<double> sample(1'000);
    model.Sample(engine, sample);

    const auto fitted = randlib::fit<randlib::GammaRand<double>>(sample);
    std::cout << fitted.distribution.Name() << " AIC=" << fitted.aic << '\n';
}
```

Every sampler accepts any standard uniform random bit generator. Fitting returns a new distribution;
it never mutates the model used to create the data.

## When to use RandLib

| Library | Sampling | PDF/CDF/quantiles | Fitting and diagnostics | Heavy-tail focus | Prefer it when… |
|---|---:|---:|---:|---:|---|
| C++ `<random>` | Yes | No | No | No | You only need the standard set of portable samplers. |
| Boost.Math | Separate Boost.Random library | Excellent, including multiprecision and configurable policies | No general fitting API | Some heavy-tailed laws | Accuracy, multiprecision, or mature special functions matter most. |
| StatsLib | Yes | Yes, with optional matrix APIs | No general fitting API | Limited | You prefer an R-like functional API or Eigen/Armadillo/Blaze integration. |
| RandLib | Yes | Yes | MLE, moments, conjugate Bayes, model selection, goodness-of-fit, censoring and truncation | Stable, Pareto, GPD/GEV and tail estimators | One object-oriented library should cover simulation through fitted heavy-tail models. |

These projects are complementary. See the standard [`<random>` inventory](https://en.cppreference.com/w/cpp/numeric/random),
[Boost.Math distribution reference](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/dist_ref/dists.html),
and [StatsLib feature list](https://github.com/kthohr/stats) for their current scope.

## Use it

Add `include/` to your compiler's include path, use the generated single header, or install with
CMake 3.21 or newer:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /your/prefix
```

```cmake
find_package(RandLib 2 REQUIRED)
target_link_libraries(your_target PRIVATE RandLib::RandLib)
```

The project also supplies a vcpkg overlay port and a Conan 2 recipe. Build controls include:

- `RANDLIB_BUILD_TESTS`, `RANDLIB_BUILD_EXAMPLES`, `RANDLIB_BUILD_BENCHMARKS`, and
  `RANDLIB_BUILD_DOCS`;
- `RANDLIB_WARNINGS_AS_ERRORS` for project-owned targets;
- `RANDLIB_PRECOMPILED` to compile common explicit instantiations once; and
- `RANDLIB_PORTABLE_SPECIAL_FUNCTIONS` to force RandLib's Bessel, beta, and zeta implementations.

## Learn more

- [API reference](https://astralord.github.io/RandLib/)
- [Distribution and numerical-method reference](docs/reference/index.md)
- [Fitting distributions](docs/fitting.md)
- [Parallel sampling and streams](docs/parallel-sampling.md)
- [Stable distributions](docs/stable-distributions.md)
- [Stable Count distribution](docs/reference/distributions/stable-count.md)
- [Accuracy and reference-data policy](docs/accuracy-policy.md)
- [v1 to v2 migration](docs/migration-v2.md)
- [Packaging and the single header](docs/packaging.md)
- [Compiled examples](examples/)
- [Legacy mathematical notes](docs/legacy/README.md)

The [changelog](CHANGELOG.md) records behavior changes and known limitations. Numerical claims are
backed by committed high-precision reference data; documented exceptions are tracked in
`tests/data/reference_tolerances.csv`.

## Development

```sh
cmake --preset dev
cmake --build --preset dev --parallel
ctest --preset dev-unit
ctest --preset dev-statistical
```

Supported CI compilers include GCC 11–14, Clang 15–18, AppleClang, and MSVC 2022. The library is MIT
licensed.
