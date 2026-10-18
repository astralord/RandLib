# Normal-inverse-gamma distribution {#normal_inverse_gamma_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="bivariate_distribution_reference.html">Bivariate distributions</a>
  <span aria-hidden="true">/</span>
  <span>Normal-inverse-gamma distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Bivariate</span>
  <span role="listitem"><strong>Family</strong>Joint laws</span>
  <span role="listitem"><strong>C++ type</strong><code>NormalInverseGammaRand</code></span>
</div>

## Definition and properties

`NormalInverseGammaRand(mu, lambda, alpha, beta)` uses a finite location \f$\mu\f$, precision
\f$\lambda>0\f$, inverse-gamma shape \f$\alpha>0\f$, and inverse-gamma rate \f$\beta>0\f$. It is
defined hierarchically by

\f[
 V\sim\operatorname{InvGamma}(\alpha,\beta),
 \qquad
 X\mid V=v\sim\mathcal N\!\left(\mu,\frac{v}{\lambda}\right).
\f]

Thus the support is \f$x\in\mathbb R\f$, \f$v>0\f$, and the joint density is

\f[
 f(x,v)=
 \frac{\beta^\alpha\sqrt{\lambda}}
      {\Gamma(\alpha)\sqrt{2\pi}}
 v^{-(\alpha+3/2)}
 \exp\!\left[-\frac{2\beta+\lambda(x-\mu)^2}{2v}\right].
\f]

The marginals are

\f[
 V\sim\operatorname{InvGamma}(\alpha,\beta),
 \qquad
 X\sim t_{2\alpha}\!\left(\mu,\sqrt{\frac{\beta}{\alpha\lambda}}\right),
\f]

where the second argument of the Student law is its scale. In particular,

\f[
 \mathbb{E}[V]=\frac{\beta}{\alpha-1}\quad(\alpha>1),
 \qquad
 \operatorname{Var}[V]=
 \frac{\beta^2}{(\alpha-1)^2(\alpha-2)}\quad(\alpha>2),
\f]

and

\f[
 \mathbb{E}[X]=\mu\quad(\alpha>1/2),
 \qquad
 \operatorname{Var}[X]=\frac{\beta}{\lambda(\alpha-1)}\quad(\alpha>1).
\f]

When both marginal variances exist, \f$\operatorname{Cov}[X,V]=0\f$ and the reported correlation
is zero. The variables are not independent. `Correlation()` returns NaN for \f$\alpha\leq2\f$,
where correlation is not defined. The joint mode returned by the API is

\f[
 \operatorname{Mode}[X,V]=
 \left(\mu,\frac{\beta}{\alpha+3/2}\right).
\f]

For the joint CDF, the change of variable \f$U=\beta/V\f$ gives

\f[
 F(x,v)=\frac1{\Gamma(\alpha)}
 \int_{\beta/v}^{\infty}
 u^{\alpha-1}e^{-u}
 \Phi\!\left((x-\mu)\sqrt{\frac{\lambda u}{\beta}}\right)du.
\f]

RandLib maps the semi-infinite interval to \f$0<t<1\f$ with
\f$u=\beta/v+t/(1-t)\f$, evaluates the integral by adaptive Gauss--Kronrod quadrature, and clamps
roundoff to the two marginal probabilities. It uses exact boundary reductions for \f$x=\mu\f$,
\f$x=\pm\infty\f$, and \f$v=\infty\f$.

For completeness, the mathematical joint characteristic function can be written using the modified
Bessel function \f$K_\alpha\f$. With

\f[
 q=\frac{s^2}{2\lambda}-it,
\f]

it is

\f[
 \phi(s,t)=e^{is\mu}
 \frac{2(\beta q)^{\alpha/2}}{\Gamma(\alpha)}
 K_\alpha\!\left(2\sqrt{\beta q}\right),
\f]

interpreted by continuity at \f$q=0\f$. RandLib does not currently expose this joint CF.

## Sampling

Sampling first draws \f$V\f$ from the inverse-gamma marginal and then returns

\f[
 X=\mu+\sqrt{\frac{V}{\lambda}}\,Z,
 \qquad Z\sim\mathcal N(0,1).
\f]

## Fitting

The distribution is the conjugate prior used by `fit_bayes<NormalRand>` for an unknown normal mean
and variance. If a prior has parameters \f$(\mu_0,\lambda_0,\alpha_0,\beta_0)\f$, a sample of size
\f$m\f$ has mean \f$\bar x\f$ and population-normalized variance
\f$s^2=m^{-1}\sum_j(x_j-\bar x)^2\f$, then the returned posterior parameters are

\f[
 \lambda_m=\lambda_0+m,
 \qquad
 \mu_m=\frac{\lambda_0\mu_0+m\bar x}{\lambda_m},
 \qquad
 \alpha_m=\alpha_0+\frac m2,
\f]

\f[
 \beta_m=\beta_0+\frac m2
 \left[s^2+\frac{\lambda_0}{\lambda_m}(\mu_0-\bar x)^2\right].
\f]

The same update is used for log-normal data after taking logarithms. There is no ordinary
maximum-likelihood fitting specialization for the four-parameter joint prior itself.
