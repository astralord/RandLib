# Stable Count distribution {#stable_count_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Stable Count distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Stable families</span>
  <span role="listitem"><strong>C++ type</strong><code>StableCountRand</code></span>
</div>

## Definition and properties

`StableCountRand<Real>(exponent, scale, location)` implements the core Stable Count law for
\f$0<\alpha<1\f$, \f$\sigma>0\f$, and finite \f$\mu\f$. If \f$X\f$ is standard Stable Count, the
returned variate is \f$Y=\mu+\sigma X\f$ and has support \f$(\mu,\infty)\f$.

Let \f$L_\alpha\f$ denote the positive \f$\alpha\f$-stable density normalized by

\f[
 \mathbb{E}\!\left[e^{-sL_\alpha}\right]=e^{-s^\alpha}.
\f]

The standard Stable Count density is

\f[
 N_\alpha(x)=
 \frac{L_\alpha(1/x)}
      {\Gamma(1/\alpha+1)\,x},
 \qquad x>0.
\f]

For \f$z=(y-\mu)/\sigma>0\f$, the location-scale density, CDF, and survival function are

\f[
 f_Y(y)=\frac{N_\alpha(z)}{\sigma},
\f]

\f[
 F_Y(y)=
 \frac{1}{\Gamma(1/\alpha+1)}
 \int_{1/z}^{\infty}\frac{L_\alpha(u)}{u}\,du,
 \qquad
 S_Y(y)=
 \frac{1}{\Gamma(1/\alpha+1)}
 \int_0^{1/z}\frac{L_\alpha(u)}{u}\,du.
\f]

Below the support, \f$f_Y(y)=F_Y(y)=0\f$ and \f$S_Y(y)=1\f$. At the finite endpoint the density is
zero. The characteristic function has no simpler general closed form; RandLib's convention is
captured by

\f[
 \phi_Y(t)=e^{i\mu t}
 \int_0^\infty e^{i\sigma t x}N_\alpha(x)\,dx.
\f]

An equivalent Wright-function representation is

\f[
 N_\alpha(x)=
 \frac{W_{-\alpha,0}(-x^\alpha)}
      {\Gamma(1/\alpha+1)}.
\f]

The positive stable scale conversion is essential. For \f$\alpha\ne1\f$, RandLib's S1 stable
characteristic function is

\f[
 \exp\!\left\{
  -\gamma^\alpha |t|^\alpha
  \left[1-i\beta\,\operatorname{sgn}(t)\tan\frac{\pi\alpha}{2}\right]
  +i\mu t
 \right\}.
\f]

Taking \f$\beta=1\f$, \f$\mu=0\f$, and

\f[
 \gamma=\left(\cos\frac{\pi\alpha}{2}\right)^{1/\alpha}
\f]

gives the required Laplace transform \f$e^{-s^\alpha}\f$. The implementation uses this scale, not
unit S1 scale, when it evaluates the reciprocal-stable definition.

At \f$\alpha=1/2\f$, the law reduces exactly to
\f$\operatorname{Gamma}(\text{shape}=3/2,\text{scale}=4)\f$. As \f$\alpha\to1^-\f$, the standard
law converges weakly to a point mass at one, so the nonstandard limit is a point mass at
\f$\mu+\sigma\f$. The continuous endpoint \f$\alpha=1\f$ is not accepted by the constructor.

## Numerical methods

The density uses the Wright series where its cancellation estimate is safe:

\f[
 W_{-\alpha,0}(-x^\alpha)
 =\sum_{k=1}^{\infty}
 \frac{(-x^\alpha)^k}{k!\,\Gamma(-\alpha k)}.
\f]

Reflection of the gamma function, compensated summation, and a finite iteration guard stabilize
that path. Regions with strong cancellation fall back to the reciprocal positive-stable density.
The CDF uses the integrated Wright series where safe and otherwise evaluates the finite-tail
integrals above.

Small upper-tail probabilities are integrated directly over \f$[0,1/z]\f$ instead of subtracting a
rounded CDF from one. Generic TOMS 748 root finding inverts the appropriate direct tail for
`Quantile` and `Quantile1m`. The \f$\alpha=1/2\f$ path dispatches to the exact Gamma formulas.

The near-one upper tail can make the Wright series strongly cancelling. Those inputs deliberately
take the slower stable-density quadrature path. This prioritizes accuracy over the speed of an
unverified asymptotic shortcut.

## Moments and shape

Every standard raw moment of order \f$q>-1\f$ exists and is

\f[
 m_q=\mathbb{E}[X^q]
 =\frac{\Gamma((q+1)/\alpha)}
        {\Gamma(q+1)\Gamma(1/\alpha)}.
\f]

Writing \f$m_j=m_q|_{q=j}\f$ and \f$v=m_2-m_1^2\f$, the reported summaries are

\f[
 \mathbb{E}[Y]=\mu+\sigma m_1,
 \qquad
 \operatorname{Var}[Y]=\sigma^2v,
\f]

\f[
 \operatorname{Skew}[Y]=
 \frac{m_3-3m_1m_2+2m_1^3}{v^{3/2}},
\f]

\f[
 \operatorname{ExKurt}[Y]=
 \frac{m_4-4m_1m_3+6m_1^2m_2-3m_1^4}{v^2}-3.
\f]

RandLib evaluates general \f$m_q\f$ through log-gamma differences to avoid intermediate overflow.
The first four half-Gamma raw moments use their exact integer values so derived shape summaries are
portable when `long double` is binary64. The class does not expose an unproved closed-form mode.

## Sampling

Sampling is correctness-first inverse-CDF sampling from one caller-owned standard random engine.
The exact half case uses the Gamma decomposition. There is no renewal cutoff, discarded large
stable draw, mutable static state, or library-owned engine, so seeded standard engines remain
reproducible.

## Fitting

`fit<StableCountRand<double>>(sample)` performs constrained maximum likelihood over exponent, scale,
and location. `FitOptions::fixed` accepts `"exponent"`, `"scale"`, and `"location"`. Initialization
uses exponent \f$1/2\f$, a location just below the sample minimum, and moment-scaled spread. Empty,
non-finite, or constant samples fail cleanly.
