# Stable distributions

`StableRand<Real>` models an alpha-stable law with characteristic exponent alpha in `[0.1, 2]`,
skewness beta in `[-1, 1]`, positive scale gamma, and location mu. Nolan's S1 parameterization is the
source-compatible default; S0 is available explicitly. Important special cases are normal (`alpha=2`), Cauchy
(`alpha=1, beta=0`), and Levy (`alpha=0.5, |beta|=1`).

```cpp
#include <randlib/randlib.hpp>

#include <random>

std::mt19937_64 engine(42);
using Stable = randlib::StableRand<double>;
Stable stable(1.5, -0.3, 2.0, 1.0, Stable::Parameterization::S0);
double draw = stable(engine);
double density = stable.f(0.0);
double probability = stable.F(0.0);
```

## Parameterization matters

S1's location is discontinuous at `alpha=1` when `beta != 0`; S0 is continuous in all four
parameters. The constructor, `GetLocation()`, and location setter use the selected convention. S0
and S1 locations for the same law are related by

```text
mu1 = mu0 - beta * gamma * tan(pi * alpha / 2)       (alpha != 1)
mu1 = mu0 - (2 / pi) * beta * gamma * log(gamma)    (alpha == 1)
```

`SetParameterization()` changes only the convention used to report the location; it does not change
the law. `SetScale()` preserves the location in the currently selected convention.

## Accuracy in v2.0

| Operation/region | Current behavior |
|---|---|
| Sampling and characteristic function | Available over the complete valid parameter range. |
| PDF/CDF away from skewed `alpha` near 1 | Compared with 40-digit mpmath references; target is `1e-9` absolute error. |
| PDF/CDF when `|alpha-1| < 0.01` and `beta != 0` | Throws `randlib::not_implemented`; sampling remains available. |
| Direct survival probabilities | Used instead of `1-F` to avoid cancellation in supported tails. |
| Far-tail quantiles | Seeded from the stable power-law asymptote before root finding. |
| Fitting | Scheduled with S0 support and the near-`alpha=1` numerical work. |

The executable tolerance registry in `tests/data/reference_tolerances.csv` is authoritative. Stable
reference generation is deliberately expensive; see [the accuracy policy](accuracy-policy.md) before
changing a result or tolerance.

Geometric stable laws use an exponential mixture of stable laws. Skewed `alpha=1` density and CDF
evaluation can be very slow in the far tails because the mixture probes the stable law at enormous
arguments. This is a documented numerical limitation, not a sampling limitation.
