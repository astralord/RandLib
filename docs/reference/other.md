# Circular and singular distributions {#other_distribution_reference}

This page covers the circular `VonMisesRand` and `WrappedExponentialRand` laws and the singular
`CantorRand` law. Their support and probability operations do not fit the ordinary real-line density
template, so each entry explicitly states its measure and CDF conventions.

## Circular distributions

For circular laws, ordinary linear moments depend on the chosen cut of the circle and RandLib
returns NaN from `Mean()`, `Variance()`, `Skewness()` and `ExcessKurtosis()`. `CircularMean()` is an
angle and `CircularVariance()` is \f$1-|\operatorname{E}e^{iX}|\f$. Densities and CDFs use the
single representative interval reported by `MinValue()` and `MaxValue()`.

### Von Mises distribution (`VonMisesRand`)

`VonMisesRand(mu, kappa)` requires \f$\kappa>0\f$ and represents angles on
\f$[\mu-\pi,\mu+\pi]\f$. Its density is

\f[
 f(x)=\frac{\exp\{\kappa\cos(x-\mu)\}}{2\pi I_0(\kappa)},
\f]

where \f$I_\nu\f$ is the modified Bessel function of the first kind. The CDF has no elementary
closed form. Below the implementation crossover \f$\kappa=10.5\f$, RandLib evaluates a backwards
trigonometric recurrence whose order is selected from \f$\kappa\f$. At larger concentration it
uses a normal/`erfc` approximation designed for the concentrated regime. The large-concentration
survival function has its own complementary-`erfc` branch; below the crossover it is \f$1-F\f$.

At integer harmonic \f$m\f$, the trigonometric moments are

\f[
 \operatorname{E}(e^{imX})=
 e^{im\mu}\frac{I_m(\kappa)}{I_0(\kappa)}.
\f]

This Bessel ratio is not the ordinary characteristic function at noninteger arguments. Accordingly,
`CF(t)` uses it only when \f$t\f$ is an integer; otherwise the shared continuous-distribution path
numerically integrates \f$e^{itx}f(x)\f$ over the representative interval.

The circular summaries are

\f[
 \operatorname{CMean}(X)=\mu,
 \qquad
 \operatorname{CVar}(X)=1-\frac{I_1(\kappa)}{I_0(\kappa)}.
\f]

The linear representative's median and mode are both \f$\mu\f$. Sampling uses the
ratio-of-uniforms rejection method with a concentration-dependent envelope scale and returns an
angle in \f$[\mu-\pi,\mu+\pi]\f$. There is currently no fitting specialization or dedicated
entropy member.

### Wrapped exponential distribution (`WrappedExponentialRand`)

`WrappedExponentialRand(lambda)` uses rate \f$\lambda>0\f$ and represents the exponential law
modulo \f$2\pi\f$ on \f$[0,2\pi]\f$:

\f[
 f(x)=\frac{\lambda e^{-\lambda x}}{1-e^{-2\pi\lambda}},
\f]

\f[
 F(x)=\frac{1-e^{-\lambda x}}{1-e^{-2\pi\lambda}},
 \qquad
 S(x)=\frac{e^{-\lambda x}-e^{-2\pi\lambda}}
 {1-e^{-2\pi\lambda}}.
\f]

`log1mexp` and `expm1` retain accuracy when \f$\lambda\f$ or \f$x\f$ is small, and the survival
function is evaluated directly. For arbitrary real \f$t\f$, the characteristic function is

\f[
 \phi(t)=
 \frac{\lambda\{1-e^{-2\pi(\lambda-it)}\}}
 {(\lambda-it)(1-e^{-2\pi\lambda})}.
\f]

At integer \f$t\f$, the wrap factor cancels and this reduces to
\f$\lambda/(\lambda-it)\f$. The circular summaries are

\f[
 \operatorname{CMean}(X)=\frac\pi2-\arctan\lambda,
 \qquad
 \operatorname{CVar}(X)=1-\frac{\lambda}{\sqrt{1+\lambda^2}}.
\f]

The representative mode is zero and its median is

\f[
 \frac{\log2-\log(1+e^{-2\pi\lambda})}{\lambda}.
\f]

Both lower- and upper-tail quantiles use closed logarithmic forms. Sampling applies the lower
quantile to one uniform variate, so it is exact and rejection-free. There is currently no fitting
specialization or dedicated entropy member.

## Singular distributions

### Cantor distribution (`CantorRand`)

The Cantor law is the distribution of

\f[
 X=2\sum_{j=1}^{\infty}\frac{B_j}{3^j},
 \qquad B_j\mathrel{\mathop{\sim}^{\mathrm{iid}}}\operatorname{Bernoulli}(1/2).
\f]

It is supported on the Cantor set in \f$[0,1]\f$: its CDF is continuous, is constant on every
removed middle-third interval, and has no Lebesgue density. `F(x)` follows the ternary interval
containing \f$x\f$, updating the corresponding binary CDF value until the interval is below the
shared numerical tolerance; it returns immediately when \f$x\f$ falls in a removed interval.

The principal summaries are

\f[
 \operatorname{E}X=\frac12,
 \qquad
 \operatorname{Var}X=\frac18,
 \qquad
 \operatorname{Skew}X=0,
 \qquad
 \operatorname{ExKurt}X=-\frac85.
\f]

The lower median is \f$1/3\f$. A density-based mode is not defined, so the singular base returns
NaN from `Mode()`, `Hazard()`, likelihood methods, and generic expectation integration.

The characteristic function has the product representation

\f[
 \phi(t)=e^{it/2}\prod_{j=1}^{\infty}\cos\!\left(\frac{t}{3^j}\right).
\f]

RandLib evaluates the first 30 factors. Sampling likewise draws 30 unbiased Bernoulli digits and
returns the corresponding truncated ternary expansion; the omitted remainder is below double's
practical resolution. Quantiles invert the Cantor CDF with the shared bracketed root finder. There
is currently no fitting specialization or entropy member.
