# Skellam distribution {#skellam_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Skellam distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Poisson family</span>
  <span role="listitem"><strong>C++ type</strong><code>SkellamRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/skellam_pmf.png" loading="lazy" alt="Representative Skellam distribution mass curves"><span class="distribution-plot__caption">Representative mass curves.</span></div><div class="distribution-plot"><img src="images/skellam_cdf.png" loading="lazy" alt="Representative Skellam distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

If \f$Y\sim\operatorname{Poisson}(\mu_1)\f$ and
\f$Z\sim\operatorname{Poisson}(\mu_2)\f$ are independent, then
\f$X=Y-Z\f$ has the Skellam law. RandLib requires \f$\mu_1,\mu_2>0\f$ and supports every integer.
Its PMF is

\f[
 P(X=k)=
 e^{-(\mu_1+\mu_2)}
 \left(\frac{\mu_1}{\mu_2}\right)^{k/2}
 I_{|k|}\!\left(2\sqrt{\mu_1\mu_2}\right).
\f]

`logP` evaluates this with the scaled log-Bessel kernel. With
\f$Q_M(a,b)\f$ the generalized Marcum Q function and
\f$P_M(a,b)=1-Q_M(a,b)\f$, the two direct tails are

\f[
 F(k)=
 \begin{cases}
 P_{-k}(\sqrt{2\mu_1},\sqrt{2\mu_2}),&k<0,\\
 Q_{k+1}(\sqrt{2\mu_2},\sqrt{2\mu_1}),&k\geq0,
 \end{cases}
\f]

\f[
 S(k)=
 \begin{cases}
 Q_{-k}(\sqrt{2\mu_1},\sqrt{2\mu_2}),&k<0,\\
 P_{k+1}(\sqrt{2\mu_2},\sqrt{2\mu_1}),&k\geq0.
 \end{cases}
\f]

Thus neither tail is obtained by subtracting the other from one.

The generating functions are

\f[
 G(z)=
 \exp\!\left[\mu_1(z-1)+\mu_2(z^{-1}-1)\right],
\f]

\f[
 \phi(t)=
 \exp\!\left[
  (\mu_1+\mu_2)(\cos t-1)
  +i(\mu_1-\mu_2)\sin t
 \right].
\f]

## Moments and shape

The principal summaries are

\f[
 \operatorname{E}X=\mu_1-\mu_2,
 \qquad
 \operatorname{Var}X=\mu_1+\mu_2,
\f]

\f[
 \operatorname{Skew}X=
 \frac{\mu_1-\mu_2}{(\mu_1+\mu_2)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=\frac1{\mu_1+\mu_2}.
\f]

The median is found by discrete probability search initialized at the mean. The mode search starts
at that median and walks through adjacent log probabilities until it reaches a maximum; merely
truncating \f$\mu_1-\mu_2\f$ is not generally correct.

## Sampling

## Fitting

Sampling subtracts independent draws from two cached `PoissonRand` objects. There is currently no
`fit<SkellamRand>` specialization or dedicated entropy member.
