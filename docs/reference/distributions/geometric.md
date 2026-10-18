# Geometric distribution {#geometric_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Geometric distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Geometric and logarithmic-series families</span>
  <span role="listitem"><strong>C++ type</strong><code>GeometricRand</code></span>
</div>

## Definition and properties

`GeometricRand(p)` counts failures before the first success. With \f$0<p\leq1\f$, \f$q=1-p\f$,
and \f$k\in\{0,1,\ldots\}\f$,

\f[
 \Pr[X=k]=p q^k,
 \qquad
 F(k)=1-q^{k+1},
 \qquad
 S(k)=q^{k+1}.
\f]

## Numerical evaluation and quantiles

The implementation evaluates \f$F\f$ with `expm1` and evaluates \f$S\f$ directly, so neither short
tail is formed by subtracting a rounded probability from one. The probability-generating and
characteristic functions are

\f[
 G(z)=\frac{p}{1-qz},
 \qquad
 \phi(t)=\frac{p}{1-qe^{it}}.
\f]

Its principal summaries are

\f[
 \mathbb{E}[X]=\frac qp,
 \qquad
 \operatorname{Var}[X]=\frac q{p^2},
 \qquad
 \operatorname{Skew}[X]=\frac{1+q}{\sqrt q},
 \qquad
 \operatorname{ExKurt}[X]=6+\frac{p^2}{q},
\f]

with the standardized moments undefined at the deterministic endpoint \f$p=1\f$. The mode is zero
and the lower median is

\f[
 \left\lceil\frac{-\log2}{\log q}\right\rceil-1,
\f]

interpreted as zero when \f$p=1\f$. In natural units its Shannon entropy is

\f[
 H(X)=-\log p-\frac qp\log q
     =\frac{-p\log p-q\log q}{p},
\f]

with the usual \f$0\log0=0\f$ limit.

## Sampling

Sampling uses the identity

\f[
 X=\left\lfloor\frac{E}{-\log q}\right\rfloor,
 \qquad E\sim\operatorname{Exp}(1),
\f]

for small \f$p\f$. Otherwise a constructed distribution inverts a cached 16-entry CDF table,
recursively shifting by 16 when a uniform falls in the remaining tail. The static convenience
sampler uses the exponential transform below a slightly lower crossover and ordinary sequential
CDF inversion above it. All three paths sample the exact law.

## Fitting

Both `fit<GeometricRand>` and `fit_moments<GeometricRand>` use the closed-form estimate

\f[
 \widehat p=\frac1{1+\bar x}.
\f]

The inherited legacy beta-prior member gives
\f$\operatorname{Beta}(a+m,b+\sum_jx_j)\f$ for \f$m\f$ observations, but the unified
`fit_bayes<GeometricRand>` specialization is not currently provided.
