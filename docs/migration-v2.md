# Migrating from RandLib 1.x to 2.x

RandLib 2 is a C++20, header-only library. The probability model is familiar, but namespace, layout,
engines, fitting, and several type-level APIs intentionally changed. Migrate one distribution at a
time and regenerate any seed-based fixtures.

| Change | 1.x | 2.x |
|---|---|---|
| Language | C++17 | C++20 |
| Headers | `distributions/.../*.h` | `randlib/.../*.hpp` or `randlib/randlib.hpp` |
| Namespace | Mostly global | `randlib` |
| Linkage | Compiled library | Header-only `RandLib::RandLib`; optional explicit-instantiation archive |
| Sampling engine | Engine stored inside each distribution | Any caller-owned standard URBG passed to `operator()` or `Sample` |
| Reseeding | `const` member changed hidden state | Non-const compatibility member; prefer reseeding the external engine |
| Bulk arguments | Primarily `std::vector` | `std::span` plus contiguous-range conveniences |
| Fitting | Mutating `X.Fit(sample)` | Non-mutating `fit<Distribution>(sample)` returning `FitResult` |
| Errors | General standard exceptions | `randlib::fit_error` with `FitErrorCode` for fitting failures |
| Support category | Unscoped values | `enum class randlib::Support` |
| Inheritance | Concrete classes could be subclassed | Concrete distribution classes are `final` |
| Engine selection | `RANDLIB_GENERATOR` and engine preprocessor macros | Explicit engine type; `DefaultEngine` is Xoshiro256++ |
| Reproducibility | v1.1 algorithms | Corrected streams; no bit-for-bit compatibility promise |

## Includes and namespace

Before:

```cpp
#include "distributions/univariate/continuous/NormalRand.h"

NormalRand<double> normal(0, 1);
```

After:

```cpp
#include <randlib/continuous/normal.hpp>

randlib::NormalRand<double> normal(0, 1);
```

The umbrella `<randlib/randlib.hpp>` is convenient for applications. Libraries sensitive to compile
time can include per-distribution headers.

## Engines and sampling

Before:

```cpp
NormalRand<double> normal(0, 1, 42);
double x = normal.Variate();
normal.Reseed(7);
```

After:

```cpp
randlib::NormalRand<double> normal(0, 1);
std::mt19937_64 engine(42);
double x = normal(engine);
engine.seed(7);
```

Bulk output accepts a span-compatible contiguous range:

```cpp
std::vector<double> values(1'000);
normal.Sample(engine, values);
```

`make_streams(seed, count)` creates deterministic per-worker engines. The CMake and preprocessor
`RANDLIB_GENERATOR` selection mechanism is removed. The legacy JKiss/JLKiss engines remain under
`randlib::legacy` only for investigations that require old engine types.

## Fitting

Before:

```cpp
NormalRand<double> normal;
normal.Fit(values);
```

After:

```cpp
auto result = randlib::fit<randlib::NormalRand<double>>(values);
const auto normal = result.distribution;
```

Use `fit_moments`, `fit_bayes`, `fit_best`, and `goodness_of_fit` for named workflows. Legacy member
functions remain deprecated until v3, but new code should not depend on mutation. Catch `fit_error`
when a fitting failure needs to be distinguished from invalid construction parameters.

## Other source changes

- `Support` values require qualification, for example `randlib::Support::Finite`.
- Concrete distributions are `final`; compose a distribution instead of inheriting from it.
- Factories replace constructors whose boolean flags were ambiguous.
- Internal numerical helpers moved under `randlib::math::detail` and are not API.
- Bernoulli sampling no longer shares cached bits, and discrete-uniform sampling uses unbiased
  rejection. These fixes, together with full-mantissa uniforms and engine corrections, change draws
  for the same seed.
- `Reseed` is non-const because mutation is visible. It only affects the compatibility engine owned
  by a distribution; explicit external engines are the preferred v2 design.

Do not update a reproducibility fixture by accepting new values blindly. Where the fixture represents
a probability function rather than a sample stream, follow the [accuracy evidence
policy](accuracy-policy.md).
