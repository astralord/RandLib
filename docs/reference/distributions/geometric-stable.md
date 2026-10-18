# Geometric-stable distribution {#geometric_stable_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Geometric-stable distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Stable families</span>
  <span role="listitem"><strong>C++ type</strong><code>GeometricStableRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/geometric_stable_pdf.png" loading="lazy" alt="Representative Geometric-stable distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/geometric_stable_cdf.png" loading="lazy" alt="Representative Geometric-stable distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`GeometricStableRand(alpha, beta, gamma, mu)` is most compactly defined by exponential
subordination. Let \f$W\sim\operatorname{Exp}(1)\f$, and conditional on \f$W=w\f$ let

\f[
 X\mid W=w\sim S_1\!\left(\alpha,\beta,\gamma w^{1/\alpha},\mu w\right).
\f]

For \f$\alpha=1\f$, the S1 scale law includes its usual logarithmic location correction. Define

\f[
 \psi(t)=
 \begin{cases}
 |\gamma t|^\alpha
 \left[1-i\beta\operatorname{sgn}(t)\tan(\pi\alpha/2)\right]-i\mu t,
 &\alpha\ne1,\\[3pt]
 \gamma|t|
 \left[1+i\beta\dfrac{2}{\pi}\operatorname{sgn}(t)\log|t|\right]-i\mu t,
 &\alpha=1.
 \end{cases}
\f]

Since \f$\mathbb{E}[e^{-W\psi}]=1/(1+\psi)\f$, the characteristic function is

\f[
 \phi(t)=\frac1{1+\psi(t)}.
\f]

The density and CDF are the exponential mixtures

\f[
 f(x)=\int_0^\infty e^{-w}
 f_{S_1(\alpha,\beta,\gamma w^{1/\alpha},\mu w)}(x)\,dw,
\f]

\f[
 F(x)=\int_0^\infty e^{-w}
 F_{S_1(\alpha,\beta,\gamma w^{1/\alpha},\mu w)}(x)\,dw.
\f]

RandLib normally evaluates these in \f$u=\log w\f$, because the conditional stable features retain
roughly constant width on that scale even in extreme tails. The omitted lower head is added
analytically with the conditional function frozen at the cutoff, while \f$w>e^7\f$ is negligible
under the exponential weight. The \f$\alpha=1\f$ path uses a separate bounded transformation so
that tiny mixing variables do not repeatedly force expensive stable far-tail evaluations.

Closed or specialized paths cover \f$\alpha=2\f$ (asymmetric Laplace), symmetric
\f$\alpha=1\f$ (a Cauchy mixture), and \f$\alpha=1/2\f$ cases. At \f$\alpha=2\f$, define

\f[
 \kappa=\frac{\sqrt{\mu^2+4\gamma^2}-\mu}{2\gamma}.
\f]

Then the density is the asymmetric-Laplace form

\f[
 f(x)=\frac{\kappa}{\gamma(1+\kappa^2)}
 \begin{cases}
  e^{x/(\kappa\gamma)},&x<0,\\
  e^{-\kappa x/\gamma},&x\geq0.
 \end{cases}
\f]

For \f$\alpha<1,\beta=1,\mu\geq0\f$ the support is \f$[0,\infty)\f$; for
\f$\alpha<1,\beta=-1,\mu\leq0\f$ it is \f$(-\infty,0]\f$; otherwise it is all of
\f$\mathbb R\f$. The stable and geometric-stable laws have the same leading power tails, so
quantile searches use the same asymptotic seed.

The mean follows

\f[
 \mathbb{E}[X]=
 \begin{cases}
  \mu,&\alpha>1,\\
  +\infty,&\alpha\leq1,\ \beta=1,\\
  -\infty,&\alpha\leq1,\ \beta=-1,\\
  \text{undefined},&\alpha\leq1,\ |\beta|<1.
 \end{cases}
\f]

Only \f$\alpha=2\f$ has finite variance:

\f[
 \operatorname{Var}[X]=\mu^2+2\gamma^2.
\f]

In that case

\f[
 \operatorname{Skew}[X]=
 \frac{2(1-\kappa^6)}{(1+\kappa^4)^{3/2}},
 \qquad
 \operatorname{ExKurt}[X]=
 6-\frac{12}{(\kappa^2+\kappa^{-2})^2}.
\f]

The \f$\alpha=2\f$ median uses the exact asymmetric-Laplace quantile. The mode is zero for
\f$\alpha=1\f$ and \f$\alpha=2\f$; other medians and modes use numerical search.

## Sampling

Sampling draws the exponential mixing variable and the corresponding standard stable variate
directly. For \f$\alpha\ne1\f$ this is

\f[
 X=\mu W+\gamma W^{1/\alpha}Z,
 \qquad Z\sim S_1(\alpha,\beta,1,0),
\f]

## Fitting

with specialized but equivalent forms at \f$\alpha=1/2\f$ and \f$\alpha=1\f$. There is currently
no `fit<GeometricStableRand>` specialization or dedicated entropy member.
