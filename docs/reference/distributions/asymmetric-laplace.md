# Asymmetric Laplace distributions {#asymmetric_laplace_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Asymmetric Laplace distributions</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Exponential-tail location families</span>
  <span role="listitem"><strong>C++ type</strong><code>AsymmetricLaplaceRand</code>, <code>CenteredAsymmetricLaplaceRand</code></span>
</div>

## Definition and properties

`AsymmetricLaplaceRand(shift, scale, asymmetry)` uses
\f$m\in\mathbb R\f$, \f$\gamma>0\f$, and \f$\kappa>0\f$. With \f$z=x-m\f$,

\f[
 f(x)=\frac{\kappa}{\gamma(1+\kappa^2)}
 \begin{cases}
  e^{z/(\gamma\kappa)},&z<0,\\
  e^{-\kappa z/\gamma},&z\geq0,
 \end{cases}
\f]

\f[
 F(x)=
 \begin{cases}
  \dfrac{\kappa^2}{1+\kappa^2}e^{z/(\gamma\kappa)},&z<0,\\
  1-\dfrac{1}{1+\kappa^2}e^{-\kappa z/\gamma},&z\geq0,
 \end{cases}
\f]

and the direct survival implementation uses one minus the first branch below \f$m\f$ and the
exponential term in the second branch above it. The characteristic function is

\f[
 \phi(t)=\frac{e^{imt}}
 {(1-i\gamma t/\kappa)(1+i\gamma\kappa t)}.
\f]

This follows from the sampling representation

\f[
 X=m+\gamma\left(\frac{E_1}{\kappa}-\kappa E_2\right),
 \qquad E_1,E_2\stackrel{\mathrm{iid}}{\sim}\operatorname{Exponential}(1),
\f]

which is also RandLib's sampler. The summaries are

\f[
 \operatorname{E}X=m+\gamma(\kappa^{-1}-\kappa),
 \qquad
 \operatorname{Var}X=\gamma^2(\kappa^{-2}+\kappa^2),
 \qquad
 \operatorname{Mode}X=m,
\f]

\f[
 \operatorname{Skew}X=
 \frac{2(1-\kappa^6)}{(1+\kappa^4)^{3/2}},
 \qquad
 \operatorname{ExKurt}X=
 6-\frac{12}{(\kappa^2+\kappa^{-2})^2},
\f]

\f[
 h(X)=1+\log\!\left[\gamma(\kappa+\kappa^{-1})\right].
\f]

The median is obtained without root finding:

\f[
 \operatorname{Median}X=
 \begin{cases}
  m+\gamma\kappa\log\!\left[(1+\kappa^{-2})/2\right],&\kappa>1,\\
  m-(\gamma/\kappa)\log\!\left[(1+\kappa^2)/2\right],&0<\kappa\leq1.
 \end{cases}
\f]

All other quantiles use the corresponding logarithm of the appropriate exponential branch.

`CenteredAsymmetricLaplaceRand(scale, asymmetry)` fixes \f$m=0\f$. “Centered” therefore refers to
the density's join and mode, not its mean unless \f$\kappa=1\f$. It exposes the exponential-family
sufficient statistic
\f$(\max(X,0),\max(-X,0))\f$ with natural parameters
\f$(-\kappa/\gamma,-1/(\gamma\kappa))\f$.

There is no free `fit` specialization for either asymmetric class. The deprecated member API can
fit the scale from the weighted absolute deviations, fit the shift as the sample quantile whose
level is \f$\kappa^2/(1+\kappa^2)\f$, and fit asymmetry from the resulting one-sided sums.
