# Rademacher distribution {#rademacher_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Rademacher distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Finite categorical laws</span>
  <span role="listitem"><strong>C++ type</strong><code>RademacherRand</code></span>
</div>

## Definition and properties

The parameter-free Rademacher law has

\f[
 \Pr[X=-1]=\Pr[X=1]=\frac12,
\f]

\f[
 F(k)=
 \begin{cases}
  0,&k<-1,\\
  1/2,&-1\leq k<1,\\
  1,&k\geq1.
 \end{cases}
\f]

Its generating functions and summaries are

\f[
 G(z)=\frac{z+z^{-1}}2,
 \qquad
 \phi(t)=\cos t,
\f]

\f[
 \mathbb{E}[X]=0,\qquad
 \operatorname{Var}[X]=1,\qquad
 \operatorname{Skew}[X]=0,\qquad
 \operatorname{ExKurt}[X]=-2.
\f]

Both support points are modes and medians under common set-valued definitions; RandLib consistently
returns the lower one, \f$-1\f$. The entropy is \f$\log2\f$.

## Sampling

Sampling maps one unbiased engine bit
to \f$-1\f$ or \f$1\f$. There are no parameters to fit.
