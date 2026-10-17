# Truncated Zipf distribution {#truncated_zipf_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="discrete_distribution_reference.html">Discrete distributions</a>
  <span aria-hidden="true">/</span>
  <span>Truncated Zipf distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Discrete</span>
  <span role="listitem"><strong>Family</strong>Discrete power-law families</span>
  <span role="listitem"><strong>C++ type</strong><code>ZipfRand</code></span>
</div>

## Definition and properties

`ZipfRand(s, n)` is the zeta power law truncated to \f$\{1,\ldots,n\}\f$. RandLib requires
\f$s>1\f$ and \f$n>0\f$. With \f$H_{n,s}=\sum_{j=1}^n j^{-s}\f$,

\f[
 P(X=k)=\frac{k^{-s}}{H_{n,s}},
 \qquad
 F(k)=\frac{H_{k,s}}{H_{n,s}}\quad(1\leq k<n).
\f]

Its generating functions are the finite sums

\f[
 G(z)=\frac1{H_{n,s}}\sum_{k=1}^n\frac{z^k}{k^s},
 \qquad
 \phi(t)=\frac1{H_{n,s}}\sum_{k=1}^n\frac{e^{ikt}}{k^s}.
\f]

The CF implementation evaluates the latter sum directly. All raw moments exist and satisfy

\f[
 \operatorname{E}(X^m)=\frac{H_{n,s-m}}{H_{n,s}}.
\f]

In particular,

\f[
 \operatorname{E}X=\frac{H_{n,s-1}}{H_{n,s}},
 \qquad
 \operatorname{Var}X=
 \frac{H_{n,s-2}}{H_{n,s}}-
 \left(\frac{H_{n,s-1}}{H_{n,s}}\right)^2.
\f]

Skewness and excess kurtosis are computed from these raw moments through order four. They are
undefined for the degenerate case \f$n=1\f$. The mode is one and the median uses the generic
discrete quantile search.

## Sampling

Sampling is exact CDF inversion. Construction caches the first \f$\min(n,16)\f$ cumulative masses;
a uniform in that prefix is resolved from the table, while a uniform in the remaining tail extends
the cumulative sum one integer at a time. There is currently no fitting specialization or dedicated
entropy member.
