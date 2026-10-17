# Mathematical reference {#mathematical_reference}

This reference describes the probability laws implemented by RandLib and the algorithms used to
evaluate and sample them. It is the current counterpart of the archived 1.x `document.tex`: the
legacy source remains useful historical material, but statements here are checked against the v2
headers and tests.

RandLib currently exposes 71 concrete distribution classes: 49 continuous, 16 discrete, three
bivariate, two circular and one singular. Related classes that share one mathematical family are
documented together. Internal base classes are discussed only when their parameterization or
algorithm affects a public law.

## How to read an entry

Each distribution entry records, when applicable:

- the RandLib class and constructor parameterization;
- support, probability density or mass, CDF and direct survival function;
- characteristic function or probability-generating function;
- mean, variance and principal higher moments, including existence conditions;
- relations to other implemented laws;
- the sampling algorithm and caller-owned random-engine behavior;
- the numerical methods behind probability functions and quantiles;
- fitting methods and important limitations.

The formulas use TeX notation and are rendered by MathJax in the generated Doxygen site. Source
Markdown retains the TeX directly, so the mathematical content does not depend on generated image
files.

## Notation and API correspondence

For a continuous random variable \f$X\f$,

\f[
 f(x)=\frac{d}{dx}\Pr[X\le x],\qquad
 F(x)=\Pr[X\le x],\qquad
 S(x)=\Pr[X>x],\qquad
 \phi(t)=\operatorname{E}[e^{itX}].
\f]

For a discrete random variable, \f$p(k)=\Pr[X=k]\f$ replaces the density. RandLib uses `f` for a
continuous density, `P` for a discrete mass, `F` for the CDF, `S` for the survival function, `Variate`
or `operator()` for one draw, and `CF` for the characteristic function. `Quantile(p)` is the lower
quantile

\f[
 Q(p)=\inf\{x:F(x)\ge p\}.
\f]

The API deliberately evaluates upper tails directly when subtraction \f$1-F(x)\f$ would lose
precision. Formulas in this reference may state \f$S=1-F\f$ as a mathematical identity even when
the implementation uses a separate series, special function or asymptotic expansion.

Unless an entry says otherwise, `Variance()` returns the second central moment, `Skewness()` the
standardized third central moment and `ExcessKurtosis()` the standardized fourth central moment
minus three. Infinity denotes a divergent moment; NaN denotes a standardized moment that is not
defined.

## Distribution families

- [Continuous distributions](continuous.md)
- [Discrete distributions](discrete.md)
- [Bivariate distributions](bivariate.md)
- [Circular and singular distributions](other.md)
- [Shared numerical methods](numerical-methods.md)

The [stable-law guide](../stable-distributions.md) gives additional detail about Nolan S0/S1
parameterizations and the supported accuracy regions. The [fitting guide](../fitting.md) documents
the common non-mutating fitting API, censoring, truncation and model diagnostics.

## Source and evidence policy

The implementation is authoritative for RandLib parameter names and conventions; a standard
mathematical source is authoritative for the law itself. When the archived notes disagree with v2,
the current entry explains the v2 behavior rather than silently copying the old text. Numerical
claims follow the [accuracy and reference-data policy](../accuracy-policy.md) and are backed by exact
identities, committed high-precision data or focused regression tests.

Plots from the archived PDF are reused only when their parameterization still matches the current
law. The formulas and algorithm descriptions are complete without those plots; a separate
visualization framework can add generated figures later.
