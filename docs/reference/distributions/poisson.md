# Poisson distribution {#poisson_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Poisson distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Poisson family</span>
  <span role="listitem"><strong>C++ type</strong><code>PoissonRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/poisson_pmf.png" loading="lazy" alt="Representative Poisson distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/poisson_cdf.png" loading="lazy" alt="Representative Poisson distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`PoissonRand(rate)` uses \f$\lambda>0\f$ and support
\f$\{0,1,\ldots\}\f$:

\f[
 P(X=k)=e^{-\lambda}\frac{\lambda^k}{k!}.
\f]

The interior log PMF uses Loader's saddle-point form,

\f[
 \log P(X=k)=
 -\operatorname{stirlerr}(k)-\operatorname{bd0}(k,\lambda)
 -\frac12\log(2\pi k),
\f]

instead of subtracting three large terms in
\f$k\log\lambda-\lambda-\log(k!)\f$. The direct tails are regularized incomplete-gamma ratios:

\f[
 F(k)=Q(k+1,\lambda),
 \qquad
 S(k)=P(k+1,\lambda),
 \qquad k\geq0.
\f]

Here \f$P(a,x)=\gamma(a,x)/\Gamma(a)\f$ and
\f$Q(a,x)=\Gamma(a,x)/\Gamma(a)\f$.

The probability-generating and characteristic functions are

\f[
 G(z)=e^{\lambda(z-1)},
 \qquad
 \phi(t)=e^{\lambda(e^{it}-1)}.
\f]

## Moments and shape

All cumulants equal \f$\lambda\f$, giving

\f[
 \operatorname{E}X=\operatorname{Var}X=\lambda,
 \qquad
 \operatorname{Skew}X=\lambda^{-1/2},
 \qquad
 \operatorname{ExKurt}X=\lambda^{-1}.
\f]

The mode is \f$\lfloor\lambda\rfloor\f$; for positive integer \f$\lambda\f$,
\f$\lambda-1\f$ is also a mode and the API reports \f$\lambda\f$. `Median()` returns the
standard accurate approximation

\f[
 \max\!\left(0,\left\lfloor\lambda+\frac13-\frac{0.02}{\lambda}\right\rfloor\right)
\f]

rather than running a probability search.

## Sampling

Sampling uses inversion around \f$\lfloor\lambda\rfloor\f$ for \f$\lambda<10\f$. For larger rates
it uses a multi-region normal/exponential rejection envelope and accepts against the exact Poisson
log-PMF ratio. Constructed objects cache the envelope constants; the static convenience function
builds the same constants for its supplied rate and engine, so both entry points sample the exact
law.

## Fitting

`fit<PoissonRand>` uses the shared closed-form MLE and method-of-moments estimate

\f[
 \widehat\lambda=\bar x.
\f]

For a shape--rate prior \f$\lambda\sim\operatorname{Gamma}(a,b)\f$, `fit_bayes` returns

\f[
 \lambda\mid x_1,\ldots,x_m\sim
 \operatorname{Gamma}\!\left(a+\sum_jx_j,\ b+m\right).
\f]

The class exposes exponential-family adjusted cross-entropy terms for KL calculations, but no
dedicated Shannon-entropy member.
