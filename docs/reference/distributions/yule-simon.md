# Yule--Simon distribution {#yule_simon_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Yule--Simon distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Discrete power-law families</span>
  <span role="listitem"><strong>C++ type</strong><code>YuleRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/yule_pmf.png" loading="lazy" alt="Representative Yule--Simon distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/yule_cdf.png" loading="lazy" alt="Representative Yule--Simon distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

For shape \f$\rho>0\f$ and \f$k\in\{1,2,\ldots\}\f$,

\f[
 \Pr[X=k]=\rho B(k,\rho+1)
 =\rho\frac{\Gamma(\rho+1)\Gamma(k)}{\Gamma(k+\rho+1)}.
\f]

## Numerical evaluation and quantiles

RandLib evaluates this expression in the log domain. Its direct tails are

\f[
 S(k)=kB(k,\rho+1),
 \qquad
 F(k)=1-kB(k,\rho+1),
\f]

with the CDF formed by `expm1` so a small lower tail is retained. The probability-generating and
characteristic functions are

\f[
 G(z)=\frac{\rho z}{\rho+1}
 {}_2F_1(1,1;\rho+2;z),
 \qquad
 \phi(t)=G(e^{it}).
\f]

The generic discrete CF evaluator uses probability summation with an accelerated oscillatory tail,
rather than evaluating the hypergeometric function. The principal summaries are

\f[
 \mathbb{E}[X]=
 \begin{cases}\dfrac{\rho}{\rho-1},&\rho>1,\\+\infty,&\rho\leq1,\end{cases}
\f]

\f[
 \operatorname{Var}[X]=
 \begin{cases}
 \dfrac{\rho^2}{(\rho-1)^2(\rho-2)},&\rho>2,\\
 +\infty,&\rho\leq2,
 \end{cases}
\f]

\f[
 \operatorname{Skew}[X]=
 \frac{(\rho+1)^2\sqrt{\rho-2}}{\rho(\rho-3)}
 \quad(\rho>3),
\f]

\f[
 \operatorname{ExKurt}[X]=
 \rho+3+
 \frac{\rho(11\rho^2-49)-22}{\rho(\rho-4)(\rho-3)}
 \quad(\rho>4).
\f]

## Moments and shape

Unavailable moments are reported as infinity. The mode is one, and the median uses the generic
discrete quantile search.

## Sampling

Sampling uses the exact mixture

\f[
 Y\sim\operatorname{Pareto}(\rho,1),
 \qquad
 X=1+\operatorname{Geometric}(1/Y).
\f]

Both constructed and static samplers implement that transformation. There is currently no fitting
specialization or dedicated entropy member.
