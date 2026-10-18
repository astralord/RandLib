# Trinomial distribution {#trinomial_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="bivariate_distribution_reference.html">Bivariate distributions</a>
  <span aria-hidden="true">/</span>
  <span>Trinomial distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Bivariate</span>
  <span role="listitem"><strong>Family</strong>Joint laws</span>
  <span role="listitem"><strong>C++ type</strong><code>TrinomialRand</code></span>
</div>

## Definition and properties

`TrinomialRand(n, p1, p2)` represents two counts from a three-category multinomial experiment. It
requires \f$n>0\f$, \f$p_1,p_2\geq0\f$, and \f$p_0=1-p_1-p_2\geq0\f$. With

\f[
 x\geq0,
 \qquad y\geq0,
 \qquad x+y\leq n,
\f]

the joint mass is

\f[
 \Pr[X=x,Y=y]=
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

## Moments and shape

\f$Y\sim\operatorname{Binomial}(n,p_2)\f$. Their joint moments are

\f[
 \mathbb{E}[X]=np_1,
 \qquad
 \mathbb{E}[Y]=np_2,
\f]

\f[
 \operatorname{Var}[X]=np_1(1-p_1),
 \qquad
 \operatorname{Var}[Y]=np_2(1-p_2),
 \qquad
 \operatorname{Cov}[X,Y]=-np_1p_2,
\f]

and, when both marginal variances are positive,

\f[
 \operatorname{Corr}[X,Y]=
 -\sqrt{\frac{p_1p_2}{(1-p_1)(1-p_2)}}.
\f]

`Correlation()` returns NaN when either marginal is deterministic; the covariance returned by the
shared bivariate base is still exactly zero when a marginal variance is zero.

The CDF conditions on whichever marginal requires the shorter outer sum. For example,

\f[
 F(x,y)=
 \sum_{j=0}^{y}
 \Pr[Y=j]\,
 \Pr\!\left[
  \operatorname{Binomial}\!\left(n-j,\frac{p_1}{1-p_2}\right)\leq x
 \right],
\f]

with the symmetric expression used when it is shorter. The conditional binomial tails are evaluated
through regularized incomplete-beta functions. Deterministic \f$p_1=1\f$ and \f$p_2=1\f$ cases
take exact branches rather than forming \f$0/0\f$.

## Sampling

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
