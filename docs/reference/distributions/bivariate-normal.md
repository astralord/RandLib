# Bivariate normal distribution {#bivariate_normal_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="bivariate_distribution_reference.html">Bivariate distributions</a>
  <span aria-hidden="true">/</span>
  <span>Bivariate normal distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Bivariate</span>
  <span role="listitem"><strong>Family</strong>Joint laws</span>
  <span role="listitem"><strong>C++ type</strong><code>BivariateNormalRand</code></span>
</div>

## Definition and properties

`BivariateNormalRand(mu1, mu2, sigma1, sigma2, rho)` uses locations
\f$\mu_1,\mu_2\in\mathbb R\f$, positive marginal scales \f$\sigma_1,\sigma_2\f$, and correlation
\f$-1\leq\rho\leq1\f$. For \f$|\rho|<1\f$, put

\f[
 a=\frac{x-\mu_1}{\sigma_1},
 \qquad
 b=\frac{y-\mu_2}{\sigma_2}.
\f]

The joint density on \f$\mathbb R^2\f$ is

\f[
 f(x,y)=
 \frac{1}{2\pi\sigma_1\sigma_2\sqrt{1-\rho^2}}
 \exp\!\left[
  -\frac{a^2-2\rho ab+b^2}{2(1-\rho^2)}
 \right].
\f]

The marginals are

\f[
 X\sim\mathcal N(\mu_1,\sigma_1^2),
 \qquad
 Y\sim\mathcal N(\mu_2,\sigma_2^2),
\f]

and

\f[
 \operatorname{E}\!\begin{bmatrix}X\\Y\end{bmatrix}
 =\begin{bmatrix}\mu_1\\\mu_2\end{bmatrix},
 \qquad
 \operatorname{Cov}\!\begin{bmatrix}X\\Y\end{bmatrix}
 =\begin{bmatrix}
   \sigma_1^2&\rho\sigma_1\sigma_2\\
   \rho\sigma_1\sigma_2&\sigma_2^2
  \end{bmatrix}.
\f]

Although no joint-CF member is exposed, the mathematical characteristic function is

\f[
 \phi(s,t)=
 \exp\!\left(
  i(s\mu_1+t\mu_2)
  -\frac12\{s^2\sigma_1^2+2\rho st\sigma_1\sigma_2+t^2\sigma_2^2\}
 \right).
\f]

For \f$\rho=0\f$, the CDF is evaluated as the product of the two normal CDFs. Otherwise RandLib
uses an angular representation due to Genz. With \f$\Phi\f$ the standard-normal CDF, define

\f[
 A=\begin{cases}
    \min\{\Phi(a),\Phi(b)\},&\rho>0,\\
    \max\{\Phi(a)-\Phi(-b),0\},&\rho<0,
   \end{cases}
\f]

and, for \f$|\rho|<1\f$,

\f[
 I=\int_{\arcsin\rho}^{\operatorname{sgn}(\rho)\pi/2}
 \exp\!\left[-\frac12
  \left\{a^2+\left(a\tan\theta-b\sec\theta\right)^2\right\}
 \right]d\theta.
\f]

Then the implementation returns

\f[
 F(x,y)=A-\frac{I}{2\pi}.
\f]

The integral is evaluated by the shared adaptive Gauss--Kronrod routine. At \f$\rho=1\f$ the law
is supported on \f$b=a\f$ and \f$F=\min\{\Phi(a),\Phi(b)\}\f$; at \f$\rho=-1\f$ it is supported
on \f$b=-a\f$ and \f$F=\max\{\Phi(a)-\Phi(-b),0\}\f$. For these singular endpoints, `f` uses
infinity on the support line and zero elsewhere as a computational representation of the Dirac
mass. `Mode()` returns \f$(\mu_1,\mu_2)\f$.

## Sampling

Sampling draws two independent standard normals \f$Z_1,Z_2\f$ and applies

\f[
 X=\mu_1+\sigma_1Z_1,
 \qquad
 Y=\mu_2+\sigma_2\left(\rho Z_1+\sqrt{1-\rho^2}\,Z_2\right).
\f]

This construction also gives the exact singular limits. There is currently no fitting
specialization for `BivariateNormalRand`.
