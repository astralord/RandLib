# Continuous distributions {#continuous_distribution_reference}

This page documents every concrete continuous distribution in RandLib. Related special cases share
one section so that their parameter transformations and common numerical implementation remain
visible. Every detailed entry is checked against the public headers and test evidence.

<div class="reference-callout">
  <strong>48 focused entries</strong>
  <span>Choose a family below; each distribution opens on its own page with definitions, formulas, numerical methods, and fitting notes.</span>
</div>

## Beta family {#continuous_beta_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="beta_distribution.html"><strong>Beta distribution</strong></a>
    <code>BetaRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="arcsine_subfamily_distribution.html"><strong>Arcsine subfamily</strong></a>
    <code>ArcsineRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="balding_nichols_distribution.html"><strong>Balding--Nichols distribution</strong></a>
    <code>BaldingNicholsRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="uniform_distribution.html"><strong>Uniform distribution</strong></a>
    <code>UniformRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Gamma family {#continuous_gamma_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="gamma_distribution.html"><strong>Gamma distribution</strong></a>
    <code>GammaRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="chi_squared_distribution.html"><strong>Chi-squared distribution</strong></a>
    <code>ChiSquaredRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="erlang_distribution.html"><strong>Erlang distribution</strong></a>
    <code>ErlangRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="exponential_distribution.html"><strong>Exponential distribution</strong></a>
    <code>ExponentialRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Transformed beta and gamma laws {#continuous_transformed_beta_gamma_laws}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="beta_prime_distribution.html"><strong>Beta-prime distribution</strong></a>
    <code>BetaPrimeRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="fisher_f_distribution.html"><strong>Fisher F distribution</strong></a>
    <code>FisherFRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="inverse_gamma_distribution.html"><strong>Inverse-gamma distribution</strong></a>
    <code>InverseGammaRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Normal family {#continuous_normal_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="normal_distribution.html"><strong>Normal distribution</strong></a>
    <code>NormalRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="log_normal_distribution.html"><strong>Log-normal distribution</strong></a>
    <code>LogNormalRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="exponentially_modified_gaussian_distribution.html"><strong>Exponentially modified Gaussian distribution</strong></a>
    <code>ExponentiallyModifiedGaussianRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Nakagami family {#continuous_nakagami_family}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="nakagami_distribution.html"><strong>Nakagami distribution</strong></a>
    <code>NakagamiRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="chi_distribution.html"><strong>Chi distribution</strong></a>
    <code>ChiRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="maxwell_boltzmann_distribution.html"><strong>Maxwell--Boltzmann distribution</strong></a>
    <code>MaxwellBoltzmannRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="rayleigh_distribution.html"><strong>Rayleigh distribution</strong></a>
    <code>RayleighRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Exponential-tail location families {#continuous_exponential_tail_location_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="asymmetric_laplace_distribution.html"><strong>Asymmetric Laplace distributions</strong></a>
    <code>AsymmetricLaplaceRand, CenteredAsymmetricLaplaceRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="laplace_distribution.html"><strong>Laplace distribution</strong></a>
    <code>LaplaceRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="logistic_distribution.html"><strong>Logistic distribution</strong></a>
    <code>LogisticRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="hyperbolic_secant_distribution.html"><strong>Hyperbolic secant distribution</strong></a>
    <code>SechRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Bounded algebraic distributions {#continuous_bounded_algebraic}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="triangular_distribution.html"><strong>Triangular distribution</strong></a>
    <code>TriangularRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="raised_cosine_distribution.html"><strong>Raised-cosine distributions</strong></a>
    <code>RaisedCosineRand, RaabGreenRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="wigner_semicircle_distribution.html"><strong>Wigner semicircle distribution</strong></a>
    <code>WignerSemicircleRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="irwin_hall_distribution.html"><strong>Irwin--Hall distribution</strong></a>
    <code>IrwinHallRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="degenerate_distribution.html"><strong>Degenerate distribution</strong></a>
    <code>DegenerateRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Extreme-value and power-tail families {#continuous_extreme_value_power_tail_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="pareto_distribution.html"><strong>Pareto distribution</strong></a>
    <code>ParetoRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="weibull_distribution.html"><strong>Weibull distribution</strong></a>
    <code>WeibullRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="gumbel_distribution.html"><strong>Gumbel distribution</strong></a>
    <code>GumbelRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="frechet_distribution.html"><strong>Frechet distribution</strong></a>
    <code>FrechetRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="generalized_extreme_value_distribution.html"><strong>Generalized extreme-value distribution</strong></a>
    <code>GeneralizedExtremeValueRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="generalized_pareto_distribution.html"><strong>Generalized Pareto distribution</strong></a>
    <code>GeneralizedParetoRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Classical heavy-tailed location families {#continuous_classical_heavy_tailed_location_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="cauchy_distribution.html"><strong>Cauchy distribution</strong></a>
    <code>CauchyRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="student_t_distribution.html"><strong>Student t distribution</strong></a>
    <code>StudentTRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="levy_distribution.html"><strong>Levy distribution</strong></a>
    <code>LevyRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="inverse_gaussian_distribution.html"><strong>Inverse-Gaussian distribution</strong></a>
    <code>InverseGaussianRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Noncentral Gaussian norms {#continuous_noncentral_gaussian_norms}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="noncentral_chi_squared_distribution.html"><strong>Noncentral chi-squared distribution</strong></a>
    <code>NoncentralChiSquaredRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="rice_distribution.html"><strong>Rice distribution</strong></a>
    <code>RiceRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="noncentral_t_distribution.html"><strong>Noncentral t distribution</strong></a>
    <code>NoncentralTRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Spectral and goodness-of-fit laws {#continuous_spectral_goodness_of_fit_laws}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="kolmogorov_distribution.html"><strong>Kolmogorov distribution</strong></a>
    <code>KolmogorovSmirnovRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="marchenko_pastur_distribution.html"><strong>Marchenko--Pastur distribution</strong></a>
    <code>MarchenkoPasturRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="planck_distribution.html"><strong>Planck distribution</strong></a>
    <code>PlanckRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
</div>

## Stable families {#continuous_stable_families}

<div class="distribution-grid">
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="stable_distribution.html"><strong>Stable distribution</strong></a>
    <code>StableRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="stable_count_distribution.html"><strong>Stable Count distribution</strong></a>
    <code>StableCountRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="holtsmark_landau_distribution.html"><strong>Holtsmark and Landau distributions</strong></a>
    <code>HoltsmarkRand, LandauRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="geometric_stable_distribution.html"><strong>Geometric-stable distribution</strong></a>
    <code>GeometricStableRand</code>
    <span class="distribution-card__meta">Density/mass and CDF plots</span>
  </div>
  <div class="distribution-card">
    <span class="distribution-card__kind">Continuous</span>
    <a href="log_stable_distribution.html"><strong>Log-stable distribution</strong></a>
    <code>LogStableRand</code>
    <span class="distribution-card__meta">Formula reference</span>
  </div>
</div>

## Related guides

- @ref mathematical_reference "Mathematical-reference conventions"
- @ref numerical_methods_reference "Shared numerical methods"
- @ref md_docs_2fitting "Fitting API"
