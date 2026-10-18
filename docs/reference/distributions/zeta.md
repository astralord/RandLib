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
 \Pr[X=k]=\frac{k^{-s}}{\zeta(s)},
\f]

where \f$\zeta\f$ is the Riemann zeta function. With the generalized harmonic number
\f$H_{k,s}=\sum_{j=1}^k j^{-s}\f$,

\f[
 F(k)=\frac{H_{k,s}}{\zeta(s)},
 \qquad
 S(k)=\frac{\zeta(s,k+1)}{\zeta(s)},
\f]

where the second numerator is the Hurwitz zeta tail. The implementation evaluates \f$S\f$ directly:
terms below \f$n=16\f$ are summed, and the rest of the tail uses a ten-correction Euler--Maclaurin
expansion. For \f$k<16\f$ it forms \f$F\f$ from the harmonic sum and otherwise obtains \f$F=1-S\f$, so
neither tail is evaluated by subtraction and no evaluation sums more than sixteen terms.

The generating functions are

\f[
 G(z)=\frac{\operatorname{Li}_s(z)}{\zeta(s)},
 \qquad
 \phi(t)=\frac{\operatorname{Li}_s(e^{it})}{\zeta(s)}.
\f]

RandLib's generic discrete CF evaluator sums outward from the median and applies a short Euler
transform to the oscillatory remainder; it does not call a complex polylogarithm.

## Moments and shape

The raw moments are

\f[
 \mathbb{E}[X^m]=
 \begin{cases}
 \dfrac{\zeta(s-m)}{\zeta(s)},&s>m+1,\\
 +\infty,&s\leq m+1.
 \end{cases}
\f]

Thus

\f[
 \mathbb{E}[X]=\frac{\zeta(s-1)}{\zeta(s)}\quad(s>2),
\f]

\f[
 \operatorname{Var}[X]=
 \frac{\zeta(s-2)}{\zeta(s)}-
 \left(\frac{\zeta(s-1)}{\zeta(s)}\right)^2
 \quad(s>3),
\f]

and the absent moments are reported as infinity. Skewness and excess kurtosis are assembled from
the first four raw zeta-ratio moments and therefore require \f$s>4\f$ and \f$s>5\f$,
respectively. The mode is one; the median is found by the generic discrete quantile search.

## Sampling

Sampling uses Devroye's rejection method with a rounded Pareto proposal of shape \f$s-1\f$. For
moderate proposals the acceptance ratio is evaluated directly; above \f$10^4\f$ it switches to
`log1p` and `expm1` to avoid cancellation. Accepted values beyond the result type saturate at its
maximum. There is currently no fitting specialization or dedicated entropy member.
