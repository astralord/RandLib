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

<div class="reference-callout">
  <strong>3 focused entries</strong>
  <span>Choose a family below; each distribution opens on its own page with definitions, formulas, numerical methods, and fitting notes.</span>
</div>

## Joint laws {#bivariate_joint_laws}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Bivariate</span>
    <a href="bivariate_normal_distribution.html"><strong>Bivariate normal distribution</strong></a>
    <code>BivariateNormalRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Bivariate</span>
    <a href="normal_inverse_gamma_distribution.html"><strong>Normal-inverse-gamma distribution</strong></a>
    <code>NormalInverseGammaRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Bivariate</span>
    <a href="trinomial_distribution.html"><strong>Trinomial distribution</strong></a>
    <code>TrinomialRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Related guides

- @ref mathematical_reference "Mathematical-reference conventions"
- @ref numerical_methods_reference "Shared numerical methods"
- @ref md_docs_2fitting "Fitting API"
