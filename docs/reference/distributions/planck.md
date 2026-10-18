# Planck distribution {#planck_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Planck distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Spectral and goodness-of-fit laws</span>
  <span role="listitem"><strong>C++ type</strong><code>PlanckRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/planck_pdf.png" loading="lazy" alt="Representative Planck distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/planck_cdf.png" loading="lazy" alt="Representative Planck distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`PlanckRand(shape, scale)` uses \f$a,b>0\f$. With \f$s=a+1\f$, its density is

\f[
 f(x)=\frac{b^s}{\Gamma(s)\zeta(s)}
 \frac{x^{s-1}}{e^{bx}-1},\qquad x>0.
\f]

Expanding the Bose--Einstein denominator gives the distribution functions

\f[
 F(x)=\frac1{\Gamma(s)\zeta(s)}
 \sum_{n=1}^{\infty}\frac{\gamma(s,nbx)}{n^s},
 \qquad
 S(x)=\frac1{\Gamma(s)\zeta(s)}
 \sum_{n=1}^{\infty}\frac{\Gamma(s,nbx)}{n^s}.
\f]

The current class evaluates \f$F\f$ by adaptive quadrature and integrates the upper tail directly
for \f$S\f$ when that tail holds less than half of the mass, otherwise returning \f$1-F\f$. When
\f$a<1\f$, the density is singular at zero, so the implementation subtracts
\f$1/(bx)\f$ from \f$1/(e^{bx}-1)\f$, integrates the regular remainder, and adds the singular
power integral analytically. `expm1` and a dedicated `log(expm1)` kernel protect the density near
zero and in the far tail.

The exact mixture sampler is

\f[
 N\sim\operatorname{Zeta}(s),\qquad
 G\sim\operatorname{Gamma}(s,\text{rate }b),\qquad
 X=G/N,
\f]

## Moments and shape

with independent \f$N\f$ and \f$G\f$. Its raw moments are

\f[
 \mathbb{E}[X^r]=\frac{\Gamma(s+r)}{\Gamma(s)b^r}
 \frac{\zeta(s+r)}{\zeta(s)},\qquad r>-a.
\f]

RandLib uses this identity for the first four raw moments and forms variance, skewness, and excess
kurtosis from the corresponding central moments. The mode is

\f[
 \operatorname{Mode}[X]=
 \begin{cases}
  0,&a\leq1,\\
  \dfrac{a+W_0(-ae^{-a})}{b},&a>1,
 \end{cases}
\f]

where \f$W_0\f$ is the principal Lambert-W branch.

Mathematically, the characteristic function is the Hurwitz-zeta ratio

\f[
 \phi(t)=\frac{\zeta(s,1-it/b)}{\zeta(s)}.
\f]

## Fitting

The implementation instead uses expectation quadrature. For \f$a<1\f$, its real part applies the
same analytic subtraction used by the CDF before integrating the remainder; the imaginary part is
regular and is integrated directly. Quantiles use TOMS 748. There is currently no
`fit<PlanckRand>` specialization or dedicated entropy member.
