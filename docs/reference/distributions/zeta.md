# Zeta distribution {#zeta_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Zeta distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Discrete power-law families</span>
  <span role="listitem"><strong>C++ type</strong><code>ZetaRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/zeta_pmf.png" loading="lazy" alt="Representative Zeta distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/zeta_cdf.png" loading="lazy" alt="Representative Zeta distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

For exponent \f$s>1\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 P(X=k)=\frac{k^{-s}}{\zeta(s)},
\f]

where \f$\zeta\f$ is the Riemann zeta function. With the generalized harmonic number
\f$H_{k,s}=\sum_{j=1}^k j^{-s}\f$,

\f[
 F(k)=\frac{H_{k,s}}{\zeta(s)},
 \qquad
 S(k)=\frac{\zeta(s,k+1)}{\zeta(s)},
\f]

where the second numerator is the Hurwitz zeta tail. Up to \f$k=10^5\f$, the implementation forms
\f$F\f$ from the harmonic sum and \f$S=1-F\f$. Beyond that threshold it evaluates \f$S\f$
directly with a ten-correction Euler--Maclaurin expansion and obtains \f$F=1-S\f$.

The generating functions are

\f[
 G(z)=\frac{\operatorname{Li}_s(z)}{\zeta(s)},
 \qquad
 \phi(t)=\frac{\operatorname{Li}_s(e^{it})}{\zeta(s)}.
\f]

RandLib's generic discrete CF evaluator sums outward from the median and applies a short Euler
transform to the oscillatory remainder; it does not call a complex polylogarithm. More generally,

\f[
 \operatorname{E}(X^m)=
 \begin{cases}
 \dfrac{\zeta(s-m)}{\zeta(s)},&s>m+1,\\
 +\infty,&s\leq m+1.
 \end{cases}
\f]

Thus

\f[
 \operatorname{E}X=\frac{\zeta(s-1)}{\zeta(s)}\quad(s>2),
\f]

\f[
 \operatorname{Var}X=
 \frac{\zeta(s-2)}{\zeta(s)}-
 \left(\frac{\zeta(s-1)}{\zeta(s)}\right)^2
 \quad(s>3),
\f]

## Moments and shape

and the absent moments are reported as infinity. Skewness and excess kurtosis are assembled from
the first four raw zeta-ratio moments and therefore require \f$s>4\f$ and \f$s>5\f$,
respectively. The mode is one; the median is found by the generic discrete quantile search.

## Sampling

Sampling uses Devroye's rejection method with a rounded Pareto proposal of shape \f$s-1\f$. For
moderate proposals the acceptance ratio is evaluated directly; above \f$10^4\f$ it switches to
`log1p` and `expm1` to avoid cancellation. Accepted values beyond the result type saturate at its
maximum. There is currently no fitting specialization or dedicated entropy member.
