# Student t distribution {#student_t_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Student t distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Classical heavy-tailed location families</span>
  <span role="listitem"><strong>C++ type</strong><code>StudentTRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/student_t_pdf.png" loading="lazy" alt="Representative Student t distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/student_t_cdf.png" loading="lazy" alt="Representative Student t distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`StudentTRand(degree, location, scale)` takes \f$\nu,\sigma>0\f$ and
\f$\mu\in\mathbb R\f$. With \f$z=(x-\mu)/\sigma\f$,

\f[
 f(x)=\frac{\Gamma((\nu+1)/2)}
 {\sigma\sqrt{\nu\pi}\,\Gamma(\nu/2)}
 \left(1+\frac{z^2}{\nu}\right)^{-(\nu+1)/2}.
\f]

For \f$z<0\f$, define \f$q=\nu/(\nu+z^2)\f$. Then

\f[
 F(x)=\frac12 I_q(\nu/2,1/2),
\f]

and the other half of the CDF and the survival function follow by symmetry. RandLib evaluates the
smaller tail directly. Near the median it uses the complementary incomplete-beta argument
\f$z^2/(\nu+z^2)\f$, while in the far tail it uses the leading log-domain term if \f$q\f$ would
underflow. Degrees one, two, and three also have simplified density branches.

The characteristic function is

\f[
 \phi(t)=e^{i\mu t}
 \frac{2^{1-\nu/2}}{\Gamma(\nu/2)}
 \left(\sqrt\nu\,\sigma|t|\right)^{\nu/2}
 K_{\nu/2}\!\left(\sqrt\nu\,\sigma|t|\right),
\f]

and is evaluated in the log domain with the modified Bessel \f$K\f$ kernel.

## Sampling

Sampling uses

\f[
 X=\mu+\sigma\frac{Z}{\sqrt{V/\nu}},\qquad
 Z\sim N(0,1),\quad V\sim\chi^2_\nu,
\f]

with the denominator generated through the Nakagami implementation. At \f$\nu=1\f$ the dedicated
Cauchy sampler is used.

## Moments and shape

The summaries are

\f[
 \mathbb{E}[X]=\mu\quad(\nu>1),\qquad
 \operatorname{Var}[X]=\frac{\nu\sigma^2}{\nu-2}\quad(\nu>2),
\f]

\f[
 \operatorname{Median}[X]=\operatorname{Mode}[X]=\mu,
 \qquad \operatorname{Skew}[X]=0\quad(\nu>3),
\f]

\f[
 \operatorname{ExKurt}[X]=\frac6{\nu-4}\quad(\nu>4).
\f]

The API distinguishes an infinite variance or fourth moment from summaries that are undefined
because a lower moment is absent. The differential entropy, although not exposed by a dedicated
member, is

\f[
 h[X]=\log\!\left[\sigma\sqrt\nu\,B(\nu/2,1/2)\right]
 +\frac{\nu+1}{2}\left[
 \psi((\nu+1)/2)-\psi(\nu/2)
 \right].
\f]

## Numerical evaluation and quantiles

Quantiles are analytic for \f$\nu=1,2,4\f$; all other degrees use the shared TOMS 748 inversion of
the direct tails.

## Fitting

`fit<StudentTRand>` performs constrained three-parameter MLE, initialized with
\f$\nu=5\f$, the sample median, and an interquartile scale estimate.
