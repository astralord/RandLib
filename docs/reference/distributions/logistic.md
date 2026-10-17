# Logistic distribution {#logistic_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Logistic distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Exponential-tail location families</span>
  <span role="listitem"><strong>C++ type</strong><code>LogisticRand</code></span>
</div>

## Definition and properties

`LogisticRand(location, scale)` uses \f$\mu\in\mathbb R\f$ and \f$s>0\f$. With
\f$z=(x-\mu)/s\f$,

\f[
 f(x)=\frac{e^{-z}}{s(1+e^{-z})^2}
 =\frac{e^{-|z|}}{s(1+e^{-|z|})^2},
\f]

\f[
 F(x)=\frac1{1+e^{-z}},\qquad
 S(x)=\frac1{1+e^z},\qquad
 Q(p)=\mu+s\log\frac{p}{1-p}.
\f]

The density uses the symmetric second form to avoid overflow, `logf` uses a softplus kernel, and
both quantile tails use `log1p`.

## Sampling

Sampling applies the logit transform to an open-interval uniform
variate. The characteristic function and summaries are

\f[
 \phi(t)=e^{i\mu t}\frac{\pi st}{\sinh(\pi st)},
\f]

\f[
 \operatorname{E}X=\operatorname{Median}X=\operatorname{Mode}X=\mu,
 \qquad
 \operatorname{Var}X=\frac{\pi^2s^2}{3},
\f]

\f[
 \operatorname{Skew}X=0,\qquad
 \operatorname{ExKurt}X=\frac65,\qquad
 h(X)=2+\log s.
\f]

## Fitting

`fit<LogisticRand>` initializes location with the sample mean and scale with
\f$\sqrt{3s_n^2}/\pi\f$, then runs the common constrained MLE optimizer. With scale fixed, the
deprecated `FitLocation` member instead solves the location score by safeguarded Newton iteration.
