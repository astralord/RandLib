#ifndef PROBABILITY_DISTRIBUTIONS_H
#define PROBABILITY_DISTRIBUTIONS_H

#include "ProbabilityDistribution.h"
#include "univariate/BasicRandGenerator.h"

/// DISCRETE
#include "univariate/discrete/DiscreteDistribution.h"
#include "univariate/discrete/BernoulliRand.h"
#include "univariate/discrete/BetaBinomialRand.h"
#include "univariate/discrete/BinomialRand.h"
#include "univariate/discrete/GeometricRand.h"
#include "univariate/discrete/HyperGeometricRand.h"
#include "univariate/discrete/NegativeBinomialRand.h"
#include "univariate/discrete/LogarithmicRand.h"
#include "univariate/discrete/PoissonRand.h"
#include "univariate/discrete/RademacherRand.h"
#include "univariate/discrete/SkellamRand.h"
#include "univariate/discrete/UniformDiscreteRand.h"
#include "univariate/discrete/YuleRand.h"
#include "univariate/discrete/ZetaRand.h"
#include "univariate/discrete/ZipfRand.h"

/// CONTINUOUS
#include "univariate/continuous/ContinuousDistribution.h"
#include "univariate/continuous/BetaPrimeRand.h"
#include "univariate/continuous/BetaRand.h"
#include "univariate/continuous/CauchyRand.h"
#include "univariate/continuous/ExponentialRand.h"
#include "univariate/continuous/ExponentialNormalRand.h"
#include "univariate/continuous/FisherSnedecorRand.h"
#include "univariate/continuous/FrechetRand.h"
#include "univariate/continuous/GammaRand.h"
#include "univariate/continuous/GeometricStableRand.h"
#include "univariate/continuous/GumbelRand.h"
#include "univariate/continuous/InverseGammaRand.h"
#include "univariate/continuous/IrwinHallRand.h"
#include "univariate/continuous/LaplaceRand.h"
#include "univariate/continuous/LevyRand.h"
#include "univariate/continuous/LogCauchyRand.h"
#include "univariate/continuous/LogisticRand.h"
#include "univariate/continuous/LogNormalRand.h"
#include "univariate/continuous/NakagamiRand.h"
#include "univariate/continuous/NormalRand.h"
#include "univariate/continuous/MaxwellBoltzmannRand.h"
#include "univariate/continuous/ParetoRand.h"
#include "univariate/continuous/PlanckRand.h"
#include "univariate/continuous/RaisedCosineRand.h"
#include "univariate/continuous/RayleighRand.h"
#include "univariate/continuous/SechRand.h"
#include "univariate/continuous/StableRand.h"
#include "univariate/continuous/StudentTRand.h"
#include "univariate/continuous/UniformRand.h"
#include "univariate/continuous/TriangularRand.h"
#include "univariate/continuous/VonMisesRand.h"
#include "univariate/continuous/WaldRand.h"
#include "univariate/continuous/WeibullRand.h"
#include "univariate/continuous/WignerSemicircleRand.h"

/// SINGULAR
#include "univariate/singular/SingularDistribution.h"
#include "univariate/singular/CantorRand.h"

/// BIVARIATE
#include "multivariate/BivariateProbabilityDistribution.h"
#include "multivariate/NormalInverseGammaRand.h"

enum DISTRIBUTION {
    /// DISCRETE
    BERNOULLI_RAND,
    BETA_BINOMIAL_RAND,
    BINOMIAL_RAND,
    GEOMETRIC_RAND,
    HYPERGEOMETRIC_RAND,
    NEGATIVEBINOMIAL_RAND,
    LOGARITHMIC_RAND,
    POISSON_RAND,
    RADEMACHER_RAND,
    SKELLAM_RAND,
    UNIFORM_DISCRETE_RAND,
    YULE_RAND,
    ZETA_RAND,
    ZIPF_RAND,
  
    ///CONTINUOUS
    ARCSINE_RAND,
    BETA_PRIME_RAND,
    BETA_RAND,
    CAUCHY_RAND,
    CHI_SQUARED_RAND,
    ERLANG_RAND,
    EXPONENTIAL_RAND,
    EXPONENTIAL_NORMAL_RAND,
    F_RAND,
    FRECHET_RAND,
    GAMMA_RAND,
    GEOMETRIC_STABLE_RAND,
    GUMBEL_RAND,
    INVERSE_GAMMA_RAND,
    IRWIN_HALL_RAND,
    LAPLACE_RAND,
    LEVY_RAND,
    LOGISTIC_RAND,
    LOGNORMAL_RAND,
    NAKAGAMI_RAND,
    NORMAL_RAND,
    MAXWELL_BOLTZMANN_RAND,
    PARETO_RAND,
    PLANCK_RAND,
    RAISED_COSINE_RAND,
    RAYLEIGH_RAND,
    SECH_RAND,
    STABLE_RAND,
    STUDENT_T_RAND,
    UNIFORM_RAND,
    TRIANGULAR_RAND,
    VON_MISES_RAND,
    WALD_RAND,
    WEIBULL_RAND,
    WIGNER_SEMICIRCLE_RAND,

    /// SINGULAR
    CANTOR_RAND,

    ///BIVARIATE
    NORMAL_INVERSE_GAMMA_RAND
};

#endif // PROBABILITY_DISTRIBUTIONS_H
