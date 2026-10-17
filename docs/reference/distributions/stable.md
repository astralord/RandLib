# Stable distribution {#stable_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Stable distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Stable families</span>
  <span role="listitem"><strong>C++ type</strong><code>StableRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/stable_pdf.png" loading="lazy" alt="Representative Stable distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/stable_cdf.png" loading="lazy" alt="Representative Stable distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`StableRand(exponent, skewness, scale, location, parameterization)` uses
\f$\alpha\in[0.1,2]\f$, \f$\beta\in[-1,1]\f$, \f$\gamma>0\f$, and a finite location. The default
is Nolan's S1 parameterization. For \f$\alpha\ne1\f$ its characteristic function is

\f[
 \phi_1(t)=
 \exp\!\left\{
  i\mu_1t-|\gamma t|^\alpha
  \left[1-i\beta\,\operatorname{sgn}(t)\tan\frac{\pi\alpha}{2}\right]
 \right\},
\f]

and for \f$\alpha=1\f$ it is

\f[
 \phi_1(t)=
 \exp\!\left\{
  i\mu_1t-\gamma|t|
  \left[1+i\beta\frac{2}{\pi}\operatorname{sgn}(t)\log|t|\right]
 \right\}.
\f]

The explicit S0 convention is continuous in all four parameters. Writing \f$u=|\gamma t|\f$,

\f[
 \phi_0(t)=
 \exp\!\left\{
  -u^\alpha+i\mu_0t
  +i\beta\,\operatorname{sgn}(t)\tan\frac{\pi\alpha}{2}(u^\alpha-u)
 \right\},
 \qquad \alpha\ne1,
\f]

\f[
 \phi_0(t)=
 \exp\!\left\{
  -u+i\mu_0t
  -i\beta\frac{2}{\pi}\operatorname{sgn}(t)u\log u
 \right\},
 \qquad \alpha=1.
\f]

The same law has locations

\f[
 \mu_1=
 \begin{cases}
  \mu_0-\beta\gamma\tan(\pi\alpha/2),&\alpha\ne1,\\
  \mu_0-\dfrac{2}{\pi}\beta\gamma\log\gamma,&\alpha=1.
 \end{cases}
\f]

`SetParameterization` changes only which location convention is reported. It does not change the
law. `SetScale` preserves the location in the currently selected convention.

Except for special cases, the density and CDF have no elementary form. They can always be defined
from the characteristic function by Fourier inversion,

\f[
 f(x)=\frac1\pi\int_0^\infty
 \operatorname{Re}\!\left(e^{-itx}\phi(t)\right)\,dt,
\f]

\f[
 F(x)=\frac12-\frac1\pi\int_0^\infty
 \frac{\operatorname{Im}\!\left(e^{-itx}\phi(t)\right)}{t}\,dt.
\f]

RandLib does not use these improper integrals indiscriminately. The Gaussian
\f$(\alpha=2)\f$, Cauchy \f$(\alpha=1,\beta=0)\f$, and one-sided Levy
\f$(\alpha=1/2,|\beta|=1)\f$ cases dispatch to their closed forms. General S1 probabilities use
Zolotarev finite-angle integral representations. The transformed integrand is sharply localized in
the tails, so the implementation finds the unique zero of its exponent, estimates the peak width,
and places logarithmically spaced quadrature breakpoints around it.

Near the origin, analytic values and convergent power series avoid an ill-conditioned integral.
Far into a heavy tail, ten-term inverse-power series are used; totally skewed short tails use their
saddle-point asymptote. Skewed \f$\alpha=1\f$ laws use a separate finite-angle representation and
an eighth-order logarithmic far-tail expansion. For S0 parameters close to \f$\alpha=1\f$, where
the S1 location conversion becomes singular, RandLib instead evaluates the continuous S0 Fourier
kernel by tanh--sinh quadrature and switches to the exact \f$\alpha=1\f$ limit at the closest
representable exponents. `S(x)` follows direct upper-tail paths rather than subtracting \f$F(x)\f$
where cancellation would destroy relative accuracy.

For \f$\alpha<1\f$ and \f$\beta=1\f$ the support is \f$[\mu_1,\infty)\f$; for
\f$\beta=-1\f$ it is \f$(-\infty,\mu_1]\f$. All other stable laws have support \f$\mathbb R\f$.
For \f$\alpha<2\f$ the power-law tails are

\f[
 P(X>x)\sim C_\alpha(1+\beta)\gamma^\alpha x^{-\alpha},
 \qquad
 P(X<-x)\sim C_\alpha(1-\beta)\gamma^\alpha x^{-\alpha},
\f]

\f[
 C_\alpha=\frac{\Gamma(\alpha)\sin(\pi\alpha/2)}{\pi}.
\f]

Translations are asymptotically negligible in these formulas. RandLib uses the corresponding
inverse tail as the initial quantile estimate, then brackets and refines against the smaller of the
two direct tails. Gaussian, Cauchy, and Levy quantiles use their exact inverse-error or tangent
forms.

## Moments and shape

The principal moments and summaries are

\f[
 \operatorname{E}X=
 \begin{cases}
  \mu_1,&\alpha>1,\\
  +\infty,&\alpha\leq1,\ \beta=1,\\
  -\infty,&\alpha\leq1,\ \beta=-1,\\
  \text{undefined},&\alpha\leq1,\ |\beta|<1,
 \end{cases}
\f]

\f[
 \operatorname{Var}X=
 \begin{cases}
  2\gamma^2,&\alpha=2,\\
  \infty,&\alpha<2.
 \end{cases}
\f]

Only the Gaussian case reports finite skewness and excess kurtosis, both zero. Symmetric laws have
median and mode \f$\mu_1\f$. For the Levy special cases the mode is
\f$\mu_1+\beta\gamma/3\f$; other asymmetric medians and modes are found numerically from the direct
probability and log-density functions.

## Sampling

Sampling uses the Chambers--Mallows--Stuck construction. For \f$\alpha\ne1\f$, let
\f$U\sim U(-\pi/2,\pi/2)\f$, \f$W\sim\operatorname{Exp}(1)\f$,

\f[
 B=\frac1\alpha\arctan\!\left(\beta\tan\frac{\pi\alpha}{2}\right),
 \qquad
 A=\left(1+\beta^2\tan^2\frac{\pi\alpha}{2}\right)^{1/(2\alpha)}.
\f]

A standard S1 variate is

\f[
 Z=A\,
 \frac{\sin(\alpha(U+B))}{(\cos U)^{1/\alpha}}
 \left[
  \frac{\cos(U-\alpha(U+B))}{W}
 \right]^{(1-\alpha)/\alpha}.
\f]

The class returns \f$\mu_1+\gamma Z\f$, with dedicated formulas for
\f$\alpha=1\f$, \f$\alpha=1/2\f$, Gaussian, Cauchy, and Levy cases. A continuous S0
\f$\alpha=1\f$ sampling limit prevents the S1 location shift from losing precision near one.

## Fitting

`fit<StableRand>(sample, Method::Quantile)` implements McCulloch's estimator from the 5th, 25th,
50th, 75th, and 95th sample percentiles. Lookup tables estimate \f$(\alpha,\beta)\f$, after which
the theoretical quartiles determine scale and location. `Method::MLE` starts from that robust
estimate and performs constrained likelihood optimization; fixed parameters are supported.
Both methods return S0 parameters. There is no dedicated entropy member.
