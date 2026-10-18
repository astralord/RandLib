# Binomial distribution {#binomial_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Binomial distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Binomial family</span>
  <span role="listitem"><strong>C++ type</strong><code>BinomialRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/binomial_pmf.png" loading="lazy" alt="Representative Binomial distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/binomial_cdf.png" loading="lazy" alt="Representative Binomial distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`BinomialRand(number, probability)` uses integer \f$n>0\f$ and \f$0\leq p\leq1\f$. For
\f$k=0,\ldots,n\f$,

\f[
 \Pr[X=k]={n\choose k}p^k(1-p)^{n-k}.
\f]

## Numerical evaluation and quantiles

RandLib evaluates interior log probabilities with a saddle-point expression based on Stirling
errors and the deviance term `bd0`; this avoids cancellation in the direct logarithm of the
binomial coefficient for large \f$n\f$. The two tails are

\f[
 F(k)=I_{1-p}(n-k,k+1),
 \qquad
 S(k)=I_p(k+1,n-k),
\f]

for \f$0\leq k<n\f$, with exact support-end values outside that range. Both use direct regularized
incomplete-beta evaluations.

The generating functions are

\f[
 G(z)=(1-p+pz)^n,\qquad
 \phi(t)=(1-p+pe^{it})^n.
\f]

Writing \f$q=1-p\f$,

\f[
 \mathbb{E}[X]=np,\qquad
 \operatorname{Var}[X]=npq,
\f]

\f[
 \operatorname{Skew}[X]=\frac{q-p}{\sqrt{npq}},
 \qquad
 \operatorname{ExKurt}[X]=\frac{1/(pq)-6}{n}.
\f]

The standardized moments apply for \f$0<p<1\f$.
The mode returned by the API is
\f$\min\{n,\lfloor(n+1)p\rfloor\}\f$; when \f$(n+1)p\f$ is an integer in the interior,
the preceding integer is also a mathematical mode. The median is selected from
\f$\lfloor np\rfloor\f$ and \f$\lceil np\rceil\f$ by evaluating the direct CDF.

The cached sampler chooses among three algorithms:

- direct Bernoulli summation for tiny \f$n\f$ or moderate \f$n\f$ close to \f$p=1/2\f$;
- a geometric waiting-time construction when \f$n\min(p,q)\f$ is small;
- the Devroye--Naderlsamanl rejection method otherwise.

The rejection branch first samples at the nearby lattice probability
\f$p_0=\lfloor n\min(p,q)\rfloor/n\f$, adds an independent residual binomial by waiting times when
needed, and reflects the result for \f$p>1/2\f$. Exact \f$p=0\f$ and \f$p=1\f$ bulk samples are
filled without consuming the engine.

The static `BinomialRand::Variate(n, p, g)` builds no tables, because callers such as the
beta-binomial, negative-hypergeometric and trinomial samplers change \f$p\f$ on every draw. For
\f$n\geq10\f$ it uses sequential inversion from zero when \f$n\min(p,q)<10\f$ and Hörmann's
transformed rejection with squeeze (BTRS, 1993) otherwise, so its expected work is bounded in
\f$n\f$; the acceptance test evaluates \f$\log(P(k)/P(\text{mode}))\f$ through Stirling's series
instead of differences of large log-factorials.

Sampling stays exact for small \f$p\f$ instead of switching to a Poisson approximation. The
total-variation distance of \f$\operatorname{Poisson}(np)\f$ from the binomial law is at most
\f$p(1-e^{-np})\f$ (Barbour and Hall, 1984), so an approximation indistinguishable in double
precision would require \f$p\lesssim 2\cdot10^{-13}\f$ even for \f$n=2^{31}\f$; there the
cached waiting-time sampler already draws a single geometric variate. Measured for \f$np\f$ from 0.1 to
\f$10^4\f$, the exact samplers cost 1.1 to 2.8 times a Poisson variate of the same mean.

## Fitting

`fit<BinomialRand>` currently requires \f$n\f$ as a fixed parameter and returns the closed-form
MLE

\f[
 \widehat p=\frac{1}{mn}\sum_{j=1}^m x_j.
\f]

With a \f$\operatorname{Beta}(a,b)\f$ prior, the conditional conjugate helper returns

\f[
 p\mid x_1,\ldots,x_m\sim
 \operatorname{Beta}\!\left(a+\sum_jx_j,\ b+mn-\sum_jx_j\right).
\f]
