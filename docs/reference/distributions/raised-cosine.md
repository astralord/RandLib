# Raised-cosine distributions {#raised_cosine_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Raised-cosine distributions</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Bounded algebraic distributions</span>
  <span role="listitem"><strong>C++ type</strong><code>RaisedCosineRand</code>, <code>RaabGreenRand</code></span>
</div>

## Definition and properties

`RaisedCosineRand(location, scale)` takes \f$\mu\in\mathbb R\f$ and \f$s>0\f$, with support
\f$[\mu-s,\mu+s]\f$. Writing \f$y=(x-\mu)/s\f$,

\f[
 f(x)=\frac{1+\cos(\pi y)}{2s},\qquad
 F(x)=\frac12\left(1+y+\frac{\sin(\pi y)}{\pi}\right),
 \qquad |y|<1,
\f]

and the direct survival formula changes the two plus signs after one to minus signs. The
characteristic function is

\f[
 \phi(t)=e^{i\mu t}\frac{\pi^2\sin(st)}{st\,[\pi^2-(st)^2]}.
\f]

At \f$st=\pi\f$ the displayed expression is indeterminate. The implementation rewrites it with a
cardinal-sine factor, preserving the analytic limit instead of evaluating \f$0/0\f$ (B33).

The mean, median, and mode all equal \f$\mu\f$, and

\f[
 \operatorname{Var}X=s^2\left(\frac13-\frac{2}{\pi^2}\right),\qquad
 \operatorname{Skew}X=0,
\f]

\f[
 \operatorname{ExKurt}X=
 \frac{6(90-\pi^4)}{5(\pi^2-6)^2}.
\f]

The CDF has no elementary inverse, so quantiles use the shared bracketed TOMS 748 solver. Sampling
uses Devroye's alternating-series accept/reject method rather than numerically inverting the CDF.
`RaabGreenRand` is the fixed \f$\mu=0,s=\pi\f$ specialization. Neither class currently has a
free `fit` specialization.
