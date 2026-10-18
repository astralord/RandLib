# Log-stable distribution {#log_stable_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Log-stable distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Stable families</span>
  <span role="listitem"><strong>C++ type</strong><code>LogStableRand</code></span>
</div>

## Definition and properties

If \f$X\f$ is the selected S0- or S1-parameterized stable law, `LogStableRand` defines

\f[
 Y=e^X.
\f]

For \f$y>0\f$,

\f[
 f_Y(y)=\frac{f_X(\log y)}{y},
 \qquad
 F_Y(y)=F_X(\log y),
 \qquad
 S_Y(y)=S_X(\log y),
\f]

and its quantiles and median are the exponentials of their stable counterparts. A one-sided
stable support transforms endpoint by endpoint: a finite lower endpoint \f$\mu_1\f$ becomes
\f$e^{\mu_1}\f$, a finite upper endpoint becomes \f$e^{\mu_1}\f$, and an infinite lower endpoint
becomes zero.

## Moments and shape

Positive integer moments exist for the Gaussian case and for non-Gaussian laws with
\f$\beta=-1\f$. Let \f$\mu_1\f$ be the S1 location. For \f$r>0\f$ their logarithms are

\f[
 \log\mathbb{E}[Y^r]=
 \begin{cases}
 r\mu_1-\dfrac{(\gamma r)^\alpha}{\cos(\pi\alpha/2)},&\alpha\ne1,\\[7pt]
 r\mu_1+\dfrac{2\gamma}{\pi}r\log r,&\alpha=1.
 \end{cases}
\f]

## Numerical evaluation and quantiles

The implementation evaluates the first four raw moments in the log domain, then forms variance,
skewness, and excess kurtosis with `expm1`-based differences. For \f$\alpha=2\f$ this reduces to
the log-normal law with underlying variance \f$2\gamma^2\f$. For other skewnesses, `Mean()` and
`Variance()` return infinity and standardized higher moments return NaN.

The Jacobian means that the mode is not generally \f$e^{\operatorname{Mode}[X]}\f$. RandLib returns
zero when \f$\alpha<2\f$ and \f$\beta<1\f$, where the transformed density diverges at the boundary.
For the Gaussian case,

\f[
 \operatorname{Mode}[Y]=e^{\mu_1-2\gamma^2}.
\f]

The remaining totally right-skewed cases minimize
\f$-\log f_X(x)+x\f$ in log space.

## Sampling

Sampling simply exponentiates a stable draw; finite-precision
underflow and overflow are therefore possible for sufficiently heavy tails.

## Fitting

`fit<LogStableRand>` requires positive observations, applies the stable McCulloch or MLE estimator
to their logarithms, and returns S0 parameters. `CF(0)` is one, but nonzero characteristic-
function arguments deliberately throw `randlib::not_implemented`: no production-quality
oscillatory transform is included. There is no dedicated entropy member.
