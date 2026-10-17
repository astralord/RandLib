# Pareto distribution {#pareto_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Pareto distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Extreme-value and power-tail families</span>
  <span role="listitem"><strong>C++ type</strong><code>ParetoRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/pareto_pdf.png" loading="lazy" alt="Representative Pareto distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/pareto_cdf.png" loading="lazy" alt="Representative Pareto distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`ParetoRand(shape, scale)` uses \f$\alpha,\sigma>0\f$ and support \f$[\sigma,\infty)\f$:

\f[
 f(x)=\frac{\alpha\sigma^\alpha}{x^{\alpha+1}},\qquad
 F(x)=1-\left(\frac{\sigma}{x}\right)^\alpha,\qquad
 S(x)=\left(\frac{\sigma}{x}\right)^\alpha.
\f]

The lower probability is evaluated with `expm1`, while the survival probability is evaluated
directly. Quantiles and the sampling identity are

\f[
 Q(p)=\sigma(1-p)^{-1/\alpha},\qquad
 X=\sigma e^{E/\alpha},\quad E\sim\operatorname{Exponential}(1).
\f]

Specialized reciprocal and reciprocal-square-root transforms are used when \f$\alpha\f$ is one or
two. The characteristic function can be written using an upper incomplete gamma function,

\f[
 \phi(t)=\alpha(-i\sigma t)^\alpha\Gamma(-\alpha,-i\sigma t),
\f]

but branch conventions make that representation awkward numerically. RandLib instead applies its
semi-infinite oscillatory quadrature directly to \f$e^{itx}f(x)\f$.

The finite summaries are

\f[
 \operatorname{E}X=\frac{\alpha\sigma}{\alpha-1}\quad(\alpha>1),
\qquad
 \operatorname{Var}X=\frac{\alpha\sigma^2}
 {\left(\alpha-1\right)^2(\alpha-2)}\quad(\alpha>2),
\f]

\f[
 \operatorname{Median}X=\sigma2^{1/\alpha},\qquad
 \operatorname{Mode}X=\sigma,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2(\alpha+1)}{\alpha-3}\sqrt{\frac{\alpha-2}{\alpha}}
 \quad(\alpha>3),
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{6(\alpha^3+\alpha^2-6\alpha-2)}
 {\alpha(\alpha-3)(\alpha-4)}\quad(\alpha>4),
 \qquad
 h(X)=\log\frac{\sigma}{\alpha}+1+\frac1\alpha.
\f]

Methods whose required moment does not exist report infinity; the variance is undefined and reports
NaN when even the mean does not exist.

## Fitting

`fit<ParetoRand>` supports MLE and UMVU estimators. Joint MLE sets

\f[
 \hat\sigma=x_{(1)},\qquad
 \hat\alpha=\frac{n}{\sum_i\log(x_i/x_{(1)})}.
\f]

The implementation also retains the historical unbiased corrections for joint or separately fixed
parameters. With \f$\sigma\f$ fixed, a Gamma \f$(a,b)\f$ shape--rate prior for \f$\alpha\f$ is
conjugate:

\f[
 \alpha\mid x\sim\operatorname{Gamma}\left(
 a+n,\ b+\sum_i\log(x_i/\sigma)\right).
\f]

The free `fit_bayes` API exposes this conditional update. Hill and Pickands estimators and
peaks-over-threshold diagnostics are separately available from `randlib/fit/tail.hpp`.
