# Bernoulli distribution {#bernoulli_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Bernoulli distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Binomial family</span>
  <span role="listitem"><strong>C++ type</strong><code>BernoulliRand</code></span>
</div>

## Definition and properties

`BernoulliRand(probability)` accepts \f$0\leq p\leq1\f$ and has support \f$\{0,1\}\f$:

\f[
 \Pr[X=k]=
 \begin{cases}
  1-p,&k=0,\\
  p,&k=1,\\
  0,&\text{otherwise}.
 \end{cases}
\f]

Consequently,

\f[
 F(k)=
 \begin{cases}
  0,&k<0,\\
  1-p,&0\leq k<1,\\
  1,&k\geq1,
 \end{cases}
 \qquad
 S(k)=
 \begin{cases}
  1,&k<0,\\
  p,&0\leq k<1,\\
  0,&k\geq1.
 \end{cases}
\f]

Its probability-generating and characteristic functions are

\f[
 G(z)=1-p+pz,\qquad
 \phi(t)=1-p+pe^{it}.
\f]

## Moments and shape

The summaries are

\f[
 \mathbb{E}[X]=p,\qquad
 \operatorname{Var}[X]=p(1-p),
\f]

\f[
 \operatorname{Skew}[X]=\frac{1-2p}{\sqrt{p(1-p)}},
 \qquad
 \operatorname{ExKurt}[X]=\frac{1}{p(1-p)}-6,
\f]

for non-degenerate \f$p\f$. The API median is zero when \f$p\leq1/2\f$ and one otherwise. Its
single reported mode is zero for \f$p<1/2\f$ and one for \f$p\geq1/2\f$; mathematically both
points are modes at equality. The entropy is

\f[
 H(X)=-p\log p-(1-p)\log(1-p),
\f]

with \f$0\log0=0\f$, so the deterministic endpoints return zero.

## Sampling

## Fitting

Sampling compares one uniform variate with \f$p\f$. Bulk sampling at \f$p=1/2\f$ instead consumes
one unbiased random bit per result. `fit<BernoulliRand>` supports closed-form MLE and
method-of-moments fitting, both with \f$\widehat p=\bar x\f$.
