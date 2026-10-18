# Categorical distribution {#categorical_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Categorical distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Finite categorical laws</span>
  <span role="listitem"><strong>C++ type</strong><code>CategoricalRand</code></span>
</div>

## Definition and properties

`CategoricalRand(probabilities)` stores a nonempty vector
\f$(p_0,\ldots,p_{K-1})\f$ with nonnegative entries summing to one. Its support is
\f$\{0,\ldots,K-1\}\f$ and

\f[
 \Pr[X=k]=p_k,\qquad
 F(k)=\sum_{j=0}^{k}p_j.
\f]

The implementation sums the shorter index side of the CDF: forward from zero when \f$k\f$ is in
the lower half, or backward as one minus the upper probabilities otherwise. Its generating
functions are finite polynomials,

\f[
 G(z)=\sum_{k=0}^{K-1}p_kz^k,
 \qquad
 \phi(t)=\sum_{k=0}^{K-1}p_ke^{itk}.
\f]

## Moments and shape

The first two moments are

\f[
 \mathbb{E}[X]=\sum_{k=0}^{K-1}kp_k,
 \qquad
 \operatorname{Var}[X]=
 \sum_{k=0}^{K-1}k^2p_k-\left(\sum_{k=0}^{K-1}kp_k\right)^2.
\f]

Higher standardized moments use the shared exact finite-support expectation sums. `Mode()`
returns the first index attaining \f$\max_kp_k\f$, and the median is the first index whose
cumulative probability reaches one half. The mathematical entropy is

\f[
 H(X)=-\sum_{k=0}^{K-1}p_k\log p_k,
 \qquad 0\log0=0,
\f]

although the class does not expose a dedicated entropy member.

## Sampling

Sampling is inverse transform with a guide table (Chen and Asau, 1974). Setting the probabilities
stores their partial sums and, for each \f$j<K\f$, the first category whose partial sum reaches
\f$j/K\f$. One uniform draw \f$U\f$ starts the search at the entry for \f$\lfloor KU\rfloor\f$, so
the expected number of comparisons is below two for every probability vector, after
\f$O(K)\f$ setup. The result equals plain sequential inversion of the same uniform, so the engine
advancement and the sampled values are unchanged from linear search.

## Fitting

There is currently no `fit<CategoricalRand>` specialization.
