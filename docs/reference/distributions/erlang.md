# Erlang distribution {#erlang_distribution}

<div class="distribution-breadcrumbs" aria-label="Breadcrumb">
  <a href="continuous_distribution_reference.html">Continuous distributions</a>
  <span aria-hidden="true">/</span>
  <span>Erlang distribution</span>
</div>

<div class="distribution-meta" role="list" aria-label="Distribution summary">
  <span role="listitem"><strong>Kind</strong>Continuous</span>
  <span role="listitem"><strong>Family</strong>Gamma family</span>
  <span role="listitem"><strong>C++ type</strong><code>ErlangRand</code></span>
</div>

## Definition and properties

`ErlangRand(shape, rate)` restricts the Gamma shape to a positive integer \f$k\f$:

\f[
 X\sim\operatorname{Erlang}(k,\beta)
 =\operatorname{Gamma}(k,\beta)
 \stackrel{d}{=}\sum_{j=1}^{k}E_j,
 \qquad E_j\sim\operatorname{Exponential}(\beta).
\f]

## Fitting

It otherwise uses all Gamma formulas and algorithms. For integer shape the sampler explicitly sums
\f$k\f$ exponential variates when \f$k\leq3\f$ and uses the general Gamma dispatch for larger
values. `fit<ErlangRand>` profiles the rate as \f$k/\bar x\f$ for every candidate integer shape and
selects the largest likelihood. The candidate ceiling is derived from the moment estimate
\f$\bar x^2/s^2\f$ and bounded by the implementation; fixed shape and rate constraints are honored.
