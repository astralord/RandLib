# Bivariate distributions {#bivariate_distribution_reference}

RandLib represents a bivariate observation as `Pair<T>`. `F({x, y})` is the lower-orthant
probability

\f[
 F(x,y)=\Pr[X\leq x,\,Y\leq y].
\f]

The inherited `S({x, y})` is \f$1-F(x,y)\f$, the probability of the complement of that lower
orthant. It is not the upper-right joint survival probability \f$\Pr[X>x,Y>y]\f$. `Mean()` returns
the two marginal means, while `Covariance()` returns the tuple
\f$(\operatorname{Var}X,\operatorname{Cov}(X,Y),\operatorname{Var}Y)\f$. The public bivariate API
does not define vector quantiles or expose a joint characteristic-function member.

## Bivariate normal distribution (`BivariateNormalRand`)

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

Sampling draws two independent standard normals \f$Z_1,Z_2\f$ and applies

\f[
 X=\mu_1+\sigma_1Z_1,
 \qquad
 Y=\mu_2+\sigma_2\left(\rho Z_1+\sqrt{1-\rho^2}\,Z_2\right).
\f]

This construction also gives the exact singular limits. There is currently no fitting
specialization for `BivariateNormalRand`.

## Normal-inverse-gamma distribution (`NormalInverseGammaRand`)

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
 \operatorname{E}V=\frac{\beta}{\alpha-1}\quad(\alpha>1),
 \qquad
 \operatorname{Var}V=
 \frac{\beta^2}{(\alpha-1)^2(\alpha-2)}\quad(\alpha>2),
\f]

and

\f[
 \operatorname{E}X=\mu\quad(\alpha>1/2),
 \qquad
 \operatorname{Var}X=\frac{\beta}{\lambda(\alpha-1)}\quad(\alpha>1).
\f]

When both marginal variances exist, \f$\operatorname{Cov}(X,V)=0\f$ and the reported correlation
is zero. The variables are not independent. `Correlation()` returns NaN for \f$\alpha\leq2\f$,
where correlation is not defined. The joint mode returned by the API is

\f[
 \operatorname{mode}(X,V)=
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

Sampling first draws \f$V\f$ from the inverse-gamma marginal and then returns

\f[
 X=\mu+\sqrt{\frac{V}{\lambda}}\,Z,
 \qquad Z\sim\mathcal N(0,1).
\f]

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

## Trinomial distribution (`TrinomialRand`)

`TrinomialRand(n, p1, p2)` represents two counts from a three-category multinomial experiment. It
requires \f$n>0\f$, \f$p_1,p_2\geq0\f$, and \f$p_0=1-p_1-p_2\geq0\f$. With

\f[
 x\geq0,
 \qquad y\geq0,
 \qquad x+y\leq n,
\f]

the joint mass is

\f[
 P(X=x,Y=y)=
 \frac{n!}{x!\,y!\,(n-x-y)!}
 p_1^x p_2^y p_0^{n-x-y}.
\f]

Zero-count factors use the convention \f$0\log0=0\f$ in `logP`; a positive count in a category of
zero probability gives log mass \f$-\infty\f$. The joint probability-generating and characteristic
functions are

\f[
 G(z_1,z_2)=(p_0+p_1z_1+p_2z_2)^n,
 \qquad
 \phi(s,t)=(p_0+p_1e^{is}+p_2e^{it})^n.
\f]

The marginals are \f$X\sim\operatorname{Binomial}(n,p_1)\f$ and
\f$Y\sim\operatorname{Binomial}(n,p_2)\f$. Their joint moments are

\f[
 \operatorname{E}X=np_1,
 \qquad
 \operatorname{E}Y=np_2,
\f]

\f[
 \operatorname{Var}X=np_1(1-p_1),
 \qquad
 \operatorname{Var}Y=np_2(1-p_2),
 \qquad
 \operatorname{Cov}(X,Y)=-np_1p_2,
\f]

and, when both marginal variances are positive,

\f[
 \operatorname{Corr}(X,Y)=
 -\sqrt{\frac{p_1p_2}{(1-p_1)(1-p_2)}}.
\f]

`Correlation()` returns NaN when either marginal is deterministic; the covariance returned by the
shared bivariate base is still exactly zero when a marginal variance is zero.

The CDF conditions on whichever marginal requires the shorter outer sum. For example,

\f[
 F(x,y)=
 \sum_{j=0}^{y}
 \Pr(Y=j)\,
 \Pr\!\left[
  \operatorname{Binomial}\!\left(n-j,\frac{p_1}{1-p_2}\right)\leq x
 \right],
\f]

with the symmetric expression used when it is shorter. The conditional binomial tails are evaluated
through regularized incomplete-beta functions. Deterministic \f$p_1=1\f$ and \f$p_2=1\f$ cases
take exact branches rather than forming \f$0/0\f$.

Sampling uses the same conditional factorization. If \f$p_1>p_2\f$, it draws

\f[
 X\sim\operatorname{Binomial}(n,p_1),
 \qquad
 Y\mid X=x\sim
 \operatorname{Binomial}\!\left(n-x,\frac{p_2}{1-p_1}\right);
\f]

otherwise it draws \f$Y\f$ first and then uses

\f[
 X\mid Y=y\sim
 \operatorname{Binomial}\!\left(n-y,\frac{p_1}{1-p_2}\right).
\f]

The deterministic category endpoints return \f$(n,0)\f$ or \f$(0,n)\f$ without consuming random
bits. `Mode()` starts from the floored mean counts and performs a log-PMF hill climb over the six
one-count transfers on the probability simplex. Multinomial log-concavity makes any resulting local
maximum a global mode; strict comparisons provide a deterministic choice when several modes tie.
There is currently no fitting or entropy specialization for `TrinomialRand`.
