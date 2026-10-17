# Noncentral t distribution {#noncentral_t_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Noncentral t distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Noncentral Gaussian norms</span>
  <span role="listitem"><strong>C++ type</strong><code>NoncentralTRand</code></span>
</div>

## Definition and properties

`NoncentralTRand(degree, noncentrality)` uses \f$\nu>0\f$ and
\f$\delta\in\mathbb R\f$. Its defining representation is

\f[
 T=\frac{Z+\delta}{W},\qquad
 Z\sim N(0,1),\quad W=\sqrt{V/\nu},\quad V\sim\chi^2_\nu,
\f]

with \f$Z\f$ and \f$V\f$ independent. The denominator density is

\f[
 g_\nu(w)=\frac{2(\nu/2)^{\nu/2}}{\Gamma(\nu/2)}
 w^{\nu-1}e^{-\nu w^2/2},\qquad w>0.
\f]

Consequently the PDF, CDF, and survival function have the conditional-normal integrals

\f[
 f(x)=\int_0^\infty w\,\varphi(xw-\delta)g_\nu(w)\,dw,
\f]

\f[
 F(x)=\int_0^\infty \Phi(xw-\delta)g_\nu(w)\,dw,
 \qquad
 S(x)=\int_0^\infty \Phi(\delta-xw)g_\nu(w)\,dw.
\f]

These definitions also show that \f$F(0)=\Phi(-\delta)\f$ and
\f$S(0)=\Phi(\delta)\f$. At \f$\delta=0\f$ the class delegates to `StudentTRand`.

For ordinary noncentral cases, RandLib evaluates both probability tails with a centered
Benton--Krishnamoorthy expansion. Its Poisson-like terms are centered at
\f$\lfloor\delta^2/2\rfloor\f$, the two incomplete-beta sequences are advanced by recurrence in
both directions, and compensated summation controls cancellation. Reflection,

\f[
 F_{\nu,\delta}(x)=S_{\nu,-\delta}(-x),
\f]

reduces negative arguments to the same positive-tail kernel. Direct conditional-normal quadrature
is used when the series center is impractically large, a result leaves \f$[0,1]\f$, or a small tail
is ill-conditioned. The integrand is evaluated in the log domain.

For
\f$|x|>\max(100,10|\delta|+10)\f$, an eight-term inverse-power expansion supplies the relevant
tail. Its coefficients are positive-part normal moments
\f$\operatorname{E}[(Z+\delta)_+^r]\f$, which are themselves evaluated by log-scaled quadrature
split at the integrand saddle. This avoids the very long incomplete-beta series required in the far
tail.

The density has three complementary paths. A confluent-hypergeometric representation is summed
directly while its terms are moderate and non-cancelling. Otherwise a recurrence expresses the
density as a scaled difference of probabilities at degrees \f$\nu\f$ and \f$\nu+2\f$. When that
difference loses significance—or when \f$x\delta<0\f$—the class evaluates the defining density
integral, split at its analytic saddle. The same eight-term asymptotic expansion covers the far-tail
density.

Conditioning on \f$W\f$ also gives the characteristic function used by the implementation:

\f[
 \phi(t)=\int_0^\infty
 \exp\!\left(-\frac{t^2}{2w^2}+i\frac{t\delta}{w}\right)
 g_\nu(w)\,dw.
\f]

Its real and imaginary parts are integrated separately, with the positive envelope split at its
mode.

## Sampling

Sampling uses the defining ratio directly through a standard normal and the Nakagami
denominator.

For \f$r<\nu\f$, define

\f[
 m_r=\operatorname{E}[W^{-r}]
 =\left(\frac\nu2\right)^{r/2}
 \frac{\Gamma((\nu-r)/2)}{\Gamma(\nu/2)}.
\f]

## Moments and shape

The first two moments are

\f[
 \operatorname{E}T=\delta m_1\quad(\nu>1),
\qquad
 \operatorname{Var}T=(1+\delta^2)m_2-\delta^2m_1^2
 \quad(\nu>2).
\f]

For compact higher-moment formulas, let

\f[
 r_1=\delta m_1,\quad
 r_2=(1+\delta^2)m_2,\quad
 r_3=(\delta^3+3\delta)m_3,\quad
 r_4=(\delta^4+6\delta^2+3)m_4.
\f]

Then, when \f$\nu>3\f$ and \f$\nu>4\f$, respectively,

\f[
 \operatorname{Skew}T=
 \frac{r_3-3r_1r_2+2r_1^3}{(r_2-r_1^2)^{3/2}},
\f]

\f[
 \operatorname{ExKurt}T=
 \frac{r_4-4r_1r_3+6r_1^2r_2-3r_1^4}{(r_2-r_1^2)^2}-3.
\f]

The implementation groups these central-moment differences with fused multiply-add operations to
reduce cancellation. As with the central t law, it distinguishes divergent moments from summaries
made undefined by a missing lower moment.

## Numerical evaluation and quantiles

## Fitting

Quantiles use TOMS 748 over the direct tails, starting at the mean when it exists and at
\f$\delta\f$ otherwise. Exact probability levels at zero return zero without root finding. The mode
uses the common log-density minimizer. `fit<NoncentralTRand>` performs constrained MLE for
\f$(\nu,\delta)\f$; it starts at \f$\nu=10\f$ unless the degree is fixed, and scales the sample
median by \f$m_1\f$ to initialize \f$\delta\f$. There is no dedicated entropy member.
