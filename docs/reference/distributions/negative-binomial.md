# Negative-binomial distribution {#negative_binomial_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Negative-binomial distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Geometric and logarithmic-series families</span>
  <span role="listitem"><strong>C++ type</strong><code>NegativeBinomialRand</code>, <code>PascalRand</code>, <code>PolyaRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/negative-binomial_pmf.png" loading="lazy" alt="Representative Negative-binomial distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/negative-binomial_cdf.png" loading="lazy" alt="Representative Negative-binomial distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

RandLib counts failures before the \f$r\f$th success. `PascalRand<Int>` is the integer-shape alias;
`PolyaRand<Int>` admits any real \f$r>0\f$. For \f$0<p\leq1\f$, \f$q=1-p\f$, and
\f$k\in\{0,1,\ldots\}\f$,

\f[
 P(X=k)=
 \frac{\Gamma(k+r)}{\Gamma(r)\,\Gamma(k+1)}p^r q^k.
\f]

The two tails are evaluated as regularized incomplete-beta functions,

\f[
 F(k)=I_p(r,k+1),
 \qquad
 S(k)=I_q(k+1,r),
\f]

and the log PMF uses `lgamma` and `lfact`. Its generating functions are

\f[
 G(z)=\left(\frac{p}{1-qz}\right)^r,
 \qquad
 \phi(t)=\left(\frac{p}{1-qe^{it}}\right)^r.
\f]

## Moments and shape

The principal summaries are

\f[
 \operatorname{E}X=\frac{rq}{p},
 \qquad
 \operatorname{Var}X=\frac{rq}{p^2},
\f]

\f[
 \operatorname{Skew}X=\frac{1+q}{\sqrt{rq}},
 \qquad
 \operatorname{ExKurt}X=\frac6r+\frac{p^2}{rq}.
\f]

The reported mode is

\f[
 \max\!\left(0,\left\lfloor\frac{(r-1)q}{p}\right\rfloor\right),
\f]

and the median is found by the generic discrete quantile search.

For integer \f$r\f$, sampling sums \f$r\f$ geometric variates, using the table or exponential
method described above. For noninteger \f$r\f$ and the remaining parameter region it uses the exact
gamma--Poisson mixture

\f[
 \Lambda\sim\operatorname{Gamma}\!\left(r,\frac p q\right),
 \qquad
 X\mid\Lambda\sim\operatorname{Poisson}(\Lambda),
\f]

where the gamma parameter shown is a rate. Overflowing geometric sums saturate at the integer
result type's maximum.

## Fitting

`fit<PolyaRand>` starts from the method-of-moments shape

\f[
 r_0=\frac{\bar x^2}{s^2-\bar x},
 \qquad
 p_0=\frac{r_0}{r_0+\bar x},
\f]

then solves the likelihood score for \f$\log r\f$ with guarded Newton steps and sets
\f$\widehat p=\widehat r/(\widehat r+\bar x)\f$. It requires \f$s^2>\bar x\f$.
With fixed \f$r\f$ and prior \f$p\sim\operatorname{Beta}(a,b)\f$, `fit_bayes` returns

\f[
 p\mid x_1,\ldots,x_m
 \sim\operatorname{Beta}\!\left(a+rm,\ b+\sum_jx_j\right).
\f]

There is no dedicated Shannon-entropy member.
