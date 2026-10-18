# Noncentral chi-squared distribution {#noncentral_chi_squared_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Noncentral chi-squared distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Noncentral Gaussian norms</span>
  <span role="listitem"><strong>C++ type</strong><code>NoncentralChiSquaredRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/noncentral_chi-squared_pdf.png" loading="lazy" alt="Representative Noncentral chi-squared distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/noncentral_chi-squared_cdf.png" loading="lazy" alt="Representative Noncentral chi-squared distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`NoncentralChiSquaredRand(degree, noncentrality)` uses \f$k>0\f$ and \f$\lambda\geq0\f$.
For \f$x>0\f$ and \f$\lambda>0\f$ its density is

\f[
 f(x)=\frac12 e^{-(x+\lambda)/2}
 \left(\frac{x}{\lambda}\right)^{k/4-1/2}
 I_{k/2-1}(\sqrt{\lambda x}).
\f]

In the central limit \f$\lambda=0\f$ this becomes the Gamma law with shape \f$k/2\f$ and rate
\f$1/2\f$:

\f[
 f(x)=\frac{x^{k/2-1}e^{-x/2}}{2^{k/2}\Gamma(k/2)}.
\f]

The distribution and survival functions use the generalized Marcum functions,

\f[
 F(x)=P_{k/2}(\sqrt\lambda,\sqrt x),\qquad
 S(x)=Q_{k/2}(\sqrt\lambda,\sqrt x).
\f]

## Numerical evaluation and quantiles

RandLib evaluates the smaller tail directly through a Poisson mixture centered at the combined
Poisson/incomplete-gamma saddle, the lower one below the mean \f$k/2+\lambda/2\f$ of the mixture and
the upper one above it; the other tail is its complement, so \f$F+S=1\f$ to rounding. The summation
proceeds outward in both directions with stable recurrences, avoiding both a long walk from index
zero and the cancellation in \f$1-F\f$. Below \f$\lambda/2=30\f$ the lower tail uses the direct
Marcum series. At
\f$\lambda=0\f$, the exact regularized incomplete-gamma limits avoid sending a singular zero
argument through the Bessel or Marcum formulas.

The characteristic function is

\f[
 \phi(t)=(1-2it)^{-k/2}
 \exp\!\left(\frac{i\lambda t}{1-2it}\right).
\f]

## Sampling

Sampling uses one of two equivalent decompositions. When \f$k\geq1\f$,

\f[
 X=(Z+\sqrt\lambda)^2+G,\qquad
 Z\sim N(0,1),\quad G\sim\chi^2_{k-1},
\f]

whereas for \f$k<1\f$,

\f[
 J\sim\operatorname{Poisson}(\lambda/2),\qquad
 X\mid J\sim\operatorname{Gamma}(k/2+J,\text{scale }2).
\f]

The exact central gamma sampler is used at \f$\lambda=0\f$. The summaries are

\f[
 \mathbb{E}[X]=k+\lambda,\qquad
 \operatorname{Var}[X]=2(k+2\lambda),
\f]

\f[
 \operatorname{Skew}[X]=
 \frac{2\sqrt2(k+3\lambda)}{(k+2\lambda)^{3/2}},
 \qquad
 \operatorname{ExKurt}[X]=
 \frac{12(k+4\lambda)}{(k+2\lambda)^2}.
\f]

For the central law the mode is \f$\max(k-2,0)\f$. In the noncentral case, the implementation
compares the boundary density with the positive mode found by numerical minimization of the negative
log density. Medians and general quantiles use TOMS 748 over the direct tails.

## Fitting

`fit<NoncentralChiSquaredRand>` performs constrained two-parameter MLE. It initializes from
\f$\bar x=k+\lambda\f$ and \f$s_n^2=2(k+2\lambda)\f$, clipping both starting parameters away from
zero. The class does not currently expose a dedicated entropy member.
