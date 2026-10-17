# Cantor distribution {#cantor_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="other_distribution_reference.html">Circular and singular distributions</a>
  <span aria-hidden="true">/</span>
  <span>Cantor distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Circular or singular</span>
  <span role="listitem"><strong>Family</strong>Singular distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>CantorRand</code></span>
</div>

## Definition and properties

The Cantor law is the distribution of

\f[
 X=2\sum_{j=1}^{\infty}\frac{B_j}{3^j},
 \qquad B_j\mathrel{\mathop{\sim}^{\mathrm{iid}}}\operatorname{Bernoulli}(1/2).
\f]

It is supported on the Cantor set in \f$[0,1]\f$: its CDF is continuous, is constant on every
removed middle-third interval, and has no Lebesgue density. `F(x)` follows the ternary interval
containing \f$x\f$, updating the corresponding binary CDF value until the interval is below the
shared numerical tolerance; it returns immediately when \f$x\f$ falls in a removed interval.

## Moments and shape

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

## Numerical evaluation and quantiles

RandLib evaluates the first 30 factors.

## Sampling

Sampling likewise draws 30 unbiased Bernoulli digits and
returns the corresponding truncated ternary expansion; the omitted remainder is below double's
practical resolution. Quantiles invert the Cantor CDF with the shared bracketed root finder. There
is currently no fitting specialization or entropy member.
