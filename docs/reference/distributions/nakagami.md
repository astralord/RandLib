# Nakagami distribution {#nakagami_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Nakagami distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Nakagami family</span>
  <span role="listitem"><strong>C++ type</strong><code>NakagamiRand</code></span>
</div>

<div class="distribution-plots"><div class="distribution-plot"><img src="images/nakagami_pdf.png" loading="lazy" alt="Representative Nakagami distribution density curves"><span class="distribution-plot__caption">Representative density curves.</span></div><div class="distribution-plot"><img src="images/nakagami_cdf.png" loading="lazy" alt="Representative Nakagami distribution cumulative distribution curves"><span class="distribution-plot__caption">Corresponding cumulative distribution functions.</span></div></div>

<div class="plot-note"><span>The legends give the parameter sets used in the archived reference plots; formulas and current constructor conventions below are authoritative.</span></div>

## Definition and properties

`NakagamiRand(shape, spread)` uses \f$\mu>0\f$ and \f$\omega>0\f$. Although the conventional
wireless-fading parameter range is often restricted to \f$\mu\geq1/2\f$, RandLib implements the
normalized law for every positive shape. Its defining transformation is

\f[
 X^2\sim\operatorname{Gamma}\!\left(\mu,\frac{\mu}{\omega}\right),
\f]

where the second Gamma parameter is a rate. Therefore, for \f$x\geq0\f$,

\f[
 f(x)=\frac{2\mu^\mu}{\Gamma(\mu)\omega^\mu}
 x^{2\mu-1}e^{-\mu x^2/\omega},
\f]

\f[
 F(x)=P\!\left(\mu,\frac{\mu x^2}{\omega}\right),
 \qquad
 S(x)=Q\!\left(\mu,\frac{\mu x^2}{\omega}\right).
\f]

The direct Gamma upper tail is used for \f$S\f$. Quantiles and samples are respectively the square
roots of Gamma quantiles and Gamma variates, so they inherit the Gamma family's tail inversion and
shape-dependent sampling dispatch.

## Moments and shape

All raw moments are summarized by

\f[
 \operatorname{E}X^r=
 \frac{\Gamma(\mu+r/2)}{\Gamma(\mu)}
 \left(\frac{\omega}{\mu}\right)^{r/2},
 \qquad r>-2\mu.
\f]

In particular,

\f[
 m=\operatorname{E}X=
 \frac{\Gamma(\mu+1/2)}{\Gamma(\mu)}
 \sqrt{\frac\omega\mu},\qquad
 \operatorname{Var}X=\omega-m^2,
\f]

\f[
 \operatorname{Mode}X=
 \begin{cases}
 \sqrt{\omega(1-1/(2\mu))},&\mu>1/2,\\
 0,&0<\mu\leq1/2.
 \end{cases}
\f]

`Skewness()` and `ExcessKurtosis()` combine the same closed raw moments through order four. The
median is the square root of the corresponding Gamma median.

The characteristic function has the hypergeometric representation

\f[
 \begin{aligned}
 \phi(t)={}&{}_1F_1\!\left(\mu;\frac12;-\frac{\omega t^2}{4\mu}\right)\\
 &+it\sqrt{\frac\omega\mu}
 \frac{\Gamma(\mu+1/2)}{\Gamma(\mu)}
 {}_1F_1\!\left(\mu+\frac12;\frac32;-\frac{\omega t^2}{4\mu}\right).
 \end{aligned}
\f]

## Numerical evaluation and quantiles

## Fitting

RandLib evaluates the defining real and imaginary expectations by adaptive infinite-interval
quadrature. When \f$\mu<1/2\f$, the density is singular at zero; the real part is then regularized as
\f$1+\int_0^\infty(\cos(tx)-1)f(x)\,dx\f$. This corrects the archived document's broader
\f$\mu<1\f$ threshold. There is currently no `fit<NakagamiRand>` specialization.
