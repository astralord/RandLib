# Inverse-gamma distribution {#inverse_gamma_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Inverse-gamma distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Transformed beta and gamma laws</span>
  <span role="listitem"><strong>C++ type</strong><code>InverseGammaRand</code></span>
</div>

## Definition and properties

`InverseGammaRand(shape, rate)` uses \f$\alpha,\beta>0\f$ and the reciprocal relationship

\f[
 X=\frac1Y,\qquad Y\sim\operatorname{Gamma}(\alpha,\beta).
\f]

For \f$x>0\f$,

\f[
 f(x)=\frac{\beta^\alpha}{\Gamma(\alpha)}x^{-\alpha-1}e^{-\beta/x},
 \qquad
 F(x)=Q\!\left(\alpha,\frac\beta x\right),
 \qquad
 S(x)=P\!\left(\alpha,\frac\beta x\right).
\f]

## Numerical evaluation and quantiles

RandLib evaluates those swapped Gamma tails directly, uses reciprocal Gamma quantiles, and samples
by taking the reciprocal of a Gamma variate. The mathematical characteristic function is

\f[
 \phi(t)=\frac{2}{\Gamma(\alpha)}
 \left(\sqrt{-i\beta t}\right)^\alpha
 K_\alpha\!\left(2\sqrt{-i\beta t}\right),
\f]

with the continuous value one at \f$t=0\f$. The implementation instead uses the shared oscillatory
Fourier integrator.

## Moments and shape

The moments and shape summaries are

\f[
 \mathbb{E}[X]=\frac\beta{\alpha-1}\quad(\alpha>1),
 \qquad
 \operatorname{Var}[X]=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}
 \quad(\alpha>2),
\f]

\f[
 \operatorname{Mode}[X]=\frac\beta{\alpha+1},\qquad
 \operatorname{Median}[X]=\frac{1}{\operatorname{Median}[Y]},
\f]

\f[
 \operatorname{Skew}[X]=\frac{4\sqrt{\alpha-2}}{\alpha-3}
 \quad(\alpha>3),
 \qquad
 \operatorname{ExKurt}[X]=\frac{30\alpha-66}{(\alpha-3)(\alpha-4)}
 \quad(\alpha>4).
\f]

## Fitting

`fit<InverseGammaRand>` initializes \f$\alpha=2+\bar x^2/s^2\f$ and
\f$\beta=\bar x(\alpha-1)\f$, then runs the common constrained maximum-likelihood optimizer. The law
also serves as the conditional conjugate prior for a normal variance when the normal location is
fixed; that model is described with the normal family.
