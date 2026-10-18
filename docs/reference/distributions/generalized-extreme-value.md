# Generalized extreme-value distribution {#generalized_extreme_value_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Generalized extreme-value distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>GeneralizedExtremeValueRand</code></span>
</div>

## Definition and properties

`GeneralizedExtremeValueRand(shape, scale, location)` uses finite \f$\xi,\mu\f$ and
\f$\sigma>0\f$. Define

\f[
 z=\frac{x-\mu}{\sigma},\qquad T(x)=1+\xi z.
\f]

For \f$\xi\ne0\f$ and \f$T(x)>0\f$,

\f[
 F(x)=e^{-T(x)^{-1/\xi}},\qquad
 f(x)=\frac1\sigma T(x)^{-1-1/\xi}e^{-T(x)^{-1/\xi}}.
\f]

The support begins at \f$\mu-\sigma/\xi\f$ when \f$\xi>0\f$ and ends there when
\f$\xi<0\f$. The continuous \f$\xi=0\f$ limit is the Gumbel law above. `log1p` and `expm1`
forms preserve this parameterization near zero, and survival probabilities are evaluated directly.
If \f$E\sim\operatorname{Exponential}(1)\f$,

\f[
 X=
 \begin{cases}
  \mu-\sigma\log E,&\xi=0,\\
  \mu+\dfrac{\sigma}{\xi}(E^{-\xi}-1),&\xi\ne0,
 \end{cases}
\f]

which supplies both the sampler and the inverse CDF after replacing \f$E\f$ by
\f$-\log p\f$. No dedicated closed-form characteristic function is used; the common expectation
quadrature evaluates \f$\phi(t)=\int e^{itx}f(x)\,dx\f$.

## Moments and shape

For \f$\xi\ne0\f$, write \f$g_r=\Gamma(1-r\xi)\f$. The finite first two moments are

\f[
 \mathbb{E}[X]=\mu+\sigma\frac{g_1-1}{\xi}\quad(\xi<1),
\qquad
 \operatorname{Var}[X]=\frac{\sigma^2}{\xi^2}(g_2-g_1^2)
 \quad(\xi<1/2).
\f]

The Gumbel limits are used exactly at \f$\xi=0\f$. For \f$\xi>-1\f$ the mode is

\f[
 \operatorname{Mode}[X]=\mu+\frac{\sigma}{\xi}
 \left[(1+\xi)^{-\xi}-1\right],
\f]

while for \f$\xi\leq-1\f$ the API returns the finite upper endpoint. Skewness exists for
\f$\xi<1/3\f$, excess kurtosis for \f$\xi<1/4\f$, and both are formed from
\f$g_1,\ldots,g_4\f$ by central-moment identities. The entropy is

\f[
 h[X]=1+\log\sigma+(1+\xi)\gamma_E.
\f]

## Fitting

`fit<GeneralizedExtremeValueRand>` uses constrained three-parameter MLE initialized at the Gumbel
moment estimate.
