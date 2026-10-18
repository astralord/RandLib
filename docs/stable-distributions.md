# Stable distributions

`StableRand<Real>` models an alpha-stable law with characteristic exponent alpha in `(0, 2]`,
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

## Accuracy

| Operation/region | Current behavior |
|---|---|
| Sampling and characteristic function | Available over the complete valid parameter range. |
| PDF/CDF | Compared with 40-digit mpmath references, including skewed inputs near `alpha=1`; target is `1e-9` absolute error. |
| Direct survival probabilities | Used instead of `1-F` to avoid cancellation in supported tails. |
| Far-tail quantiles | Seeded from the stable power-law asymptote before root finding. |
| Fitting | McCulloch quantile fitting and constrained MLE return S0 parameters. |

The executable tolerance registry in `tests/data/reference_tolerances.csv` is authoritative. Stable
reference generation is deliberately expensive; see [the accuracy policy](accuracy-policy.md) before
changing a result or tolerance.

Geometric stable laws use an exponential mixture of stable laws. Alpha-one tail expansions keep
that mixture practical even when it probes stable probabilities at enormous arguments.

## Log-stable distributions

`LogStableRand<Real>` models `Y = exp(X)` for a stable `X` and accepts the same five constructor
arguments, including the S0/S1 selection. Density, probability, sampling, quantiles and the median
are evaluated through that monotone transform. The mode is not generally `exp(X.Mode())`: the
Jacobian contributes a factor `1 / Y`. For `alpha < 2` and `beta < 1`, the polynomial left tail
makes the transformed density diverge at zero, so the mode is zero. Totally right-skewed cases use
a log-domain numerical maximum, and the Gaussian case has its closed form.

```cpp
using LogStable = randlib::LogStableRand<double>;
LogStable loss(1.5, -1.0, 0.5, 0.25, LogStable::Parameterization::S0);
double draw = loss(engine);
double upper_tail = loss.S(10.0);
auto estimate = randlib::fit<LogStable>(sample, randlib::Method::Quantile);
```

Positive moments of a non-Gaussian log-stable law exist only for `beta = -1`. Writing `mu1` for the
S1 location, their logarithms are

\f[
 \log \mathbb{E}[Y^r]=
 \begin{cases}
  r\mu_1-(\gamma r)^\alpha/\cos(\pi\alpha/2),&\alpha\ne1,\\
  r\mu_1+(2\gamma/\pi)r\log r,&\alpha=1.
 \end{cases}
\f]

The Gaussian `alpha = 2` case has moments for every beta because beta does not change that law. For
other non-Gaussian skewnesses, `Mean()` and `Variance()` return infinity and standardized higher
moments return NaN. `CF(0)` is one as for every distribution; `CF(t)` for nonzero `t` throws
`randlib::not_implemented` because no production-quality log-stable characteristic-function
algorithm is included in v2.1.

As with any exponential transform, finite-precision sampling can underflow to zero or overflow to
infinity when the stable draw lies outside the logarithmic range of `Real`. This is observable for
very heavy tails and is not a point mass in the mathematical distribution.
