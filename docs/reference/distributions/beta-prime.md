# Beta-prime distribution {#beta_prime_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Beta-prime distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Transformed beta and gamma laws</span>
  <span role="listitem"><strong>C++ type</strong><code>BetaPrimeRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/beta_prime_pdf.png" loading="lazy" alt="Representative Beta-prime distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/beta_prime_cdf.png" loading="lazy" alt="Representative Beta-prime distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`BetaPrimeRand(alpha, beta)` accepts \f$\alpha,\beta>0\f$ and has support \f$[0,\infty)\f$. Its
density and tail probabilities are

\f[
 f(x)=\frac{x^{\alpha-1}(1+x)^{-\alpha-\beta}}{B(\alpha,\beta)},
 \qquad
 F(x)=I_{x/(1+x)}(\alpha,\beta),
\f]

\f[
 S(x)=I_{1/(1+x)}(\beta,\alpha).
\f]

The direct survival expression is important numerically: for large \f$x\f$, forming
\f$x/(1+x)\f$ can round to one. RandLib instead passes \f$1/(1+x)\f$ and separately computed
logarithms to the incomplete-beta implementation.

The transformation

\f[
 Y=\frac{X}{1+X}\sim\operatorname{Beta}(\alpha,\beta),
 \qquad
 X=\frac{Y}{1-Y}
\f]

drives probability evaluation and sampling. A Beta quantile transformed by \f$y/(1-y)\f$ supplies
the initial guess for a final TOMS 748 inversion; that refinement avoids losing upper-tail accuracy
in the transformation itself.

## Sampling

Sampling applies the same transform to the selected Beta sampler.

The mathematical characteristic function can be written with Tricomi's confluent hypergeometric
function,

\f[
 \phi(t)=\frac{\Gamma(\alpha+\beta)}{\Gamma(\beta)}
 U(\alpha,1-\beta,-it).
\f]

RandLib does not evaluate \f$U\f$; it uses the shared half-period Fourier integration and Wynn
epsilon acceleration on the density over \f$[0,\infty)\f$.

## Moments and shape

The principal moments, with their existence conditions, are

\f[
 \mathbb{E}[X]=\frac{\alpha}{\beta-1}\quad(\beta>1),
\qquad
 \operatorname{Var}[X]=\frac{\alpha(\alpha+\beta-1)}
 { (\beta-2)(\beta-1)^2}\quad(\beta>2),
\f]

\f[
 \operatorname{Mode}[X]=\max\!\left\{\frac{\alpha-1}{\beta+1},0\right\},
\f]

\f[
 \operatorname{Skew}[X]=
 \frac{2(2\alpha+\beta-1)}{\beta-3}
 \sqrt{\frac{\beta-2}{\alpha(\alpha+\beta-1)}}
 \quad(\beta>3),
\f]

\f[
 \operatorname{ExKurt}[X]=
 \frac{6}{(\beta-3)(\beta-4)}
 \left[
 \frac{(\beta-1)^2(\beta-2)}{\alpha(\alpha+\beta-1)}+5\beta-11
 \right]
 \quad(\beta>4).
\f]

RandLib returns infinity when a listed moment does not exist. The median is one when
\f$\alpha=\beta\f$ and is otherwise found through the quantile path.

## Fitting

`fit<BetaPrimeRand>` performs MLE by transforming every observation to
\f$Y_i=X_i/(1+X_i)\f$ and applying the Beta score-equation solver to
\f$\overline{\log Y}\f$ and \f$\overline{\log(1-Y)}\f$.
