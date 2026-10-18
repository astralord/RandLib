# Hyperbolic secant distribution {#hyperbolic_secant_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Hyperbolic secant distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Exponential-tail location families</span>
  <span role="listitem"><strong>C++ type</strong><code>SechRand</code></span>
</div>

## Definition and properties

`SechRand` is parameter-free and standardized to unit variance:

\f[
 f(x)=\frac1{2\cosh(\pi x/2)},\qquad
 F(x)=\frac2\pi\arctan(e^{\pi x/2}),\qquad
 S(x)=\frac2\pi\arctan(e^{-\pi x/2}),
\f]

\f[
 Q(p)=\frac2\pi\log\tan\frac{\pi p}{2},\qquad
 \phi(t)=\frac1{\cosh t}.
\f]

The log density uses softplus rather than taking the logarithm of a potentially overflowing
hyperbolic cosine.

## Sampling

Sampling uses
\f$(2/\pi)\log|C|\f$ for a standard Cauchy variate \f$C\f$. Finally,

\f[
 \mathbb{E}[X]=\operatorname{Median}[X]=\operatorname{Mode}[X]=0,\qquad
 \operatorname{Var}[X]=1,\qquad
 \operatorname{Skew}[X]=0,\qquad
 \operatorname{ExKurt}[X]=2,
\f]

\f[
 h[X]=\frac{4G}{\pi},
\f]

where \f$G\f$ is Catalan's constant. There are no parameters to fit.
