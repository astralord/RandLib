# Discrete-uniform distribution {#discrete_uniform_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Discrete-uniform distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Finite categorical laws</span>
  <span role="listitem"><strong>C++ type</strong><code>UniformDiscreteRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/uniform_discrete_pmf.png" loading="lazy" alt="Representative Discrete-uniform distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/uniform_discrete_cdf.png" loading="lazy" alt="Representative Discrete-uniform distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`UniformDiscreteRand(minimum, maximum)` uses integer bounds \f$a\leq b\f$. Set
\f$n=b-a+1\f$. On \f$\{a,\ldots,b\}\f$,

\f[
 P(X=k)=\frac1n,
 \qquad
 F(k)=\frac{k-a+1}{n}\quad(a\leq k\leq b).
\f]

The probability-generating function and characteristic function are

\f[
 G(z)=\frac{z^a(1-z^n)}{n(1-z)},
\f]

\f[
 \phi(t)=
 e^{i(a+b)t/2}\frac{\sin(nt/2)}{n\sin(t/2)}.
\f]

At integer multiples of \f$2\pi\f$, the continuous limit is exactly one. RandLib first reduces the
angle modulo \f$2\pi\f$, returns that exact limit at a lattice period, and otherwise evaluates the
centered sinc form. This avoids both an indeterminate \f$0/0\f$ and precision loss from subtracting nearby
trigonometric values.

For \f$n>1\f$,

\f[
 \operatorname{E}X=\frac{a+b}{2},
 \qquad
 \operatorname{Var}X=\frac{n^2-1}{12},
\f]

\f[
 \operatorname{Skew}X=0,
 \qquad
 \operatorname{ExKurt}X=
 -\frac{6(n^2+1)}{5(n^2-1)}.
\f]

The entropy is \f$\log n\f$. Every support point is a mathematical mode; the API reports an
integer midpoint, and reports the lower integer midpoint as its median when two central points
exist. The \f$n=1\f$ law is degenerate, with zero variance and undefined standardized fourth
moment.

## Sampling

## Fitting

Sampling uses rejection-based `uniformBelow`, so reducing raw engine output modulo \f$n\f$ cannot
introduce bias. `fit<UniformDiscreteRand>` provides the closed-form MLE

\f[
 \widehat a=\min_jx_j,\qquad
 \widehat b=\max_jx_j.
\f]
