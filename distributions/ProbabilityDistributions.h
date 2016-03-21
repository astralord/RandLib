#ifndef PROBABILITY_DISTRIBUTIONS_H
#define PROBABILITY_DISTRIBUTIONS_H

#include "ProbabilityDistribution.h"
#include "BasicRandGenerator.h"

/// DISCRETE
#include "discrete/DiscreteDistribution.h"
#include "discrete/BernoulliRand.h"
#include "discrete/BinomialRand.h"
#include "discrete/GeometricRand.h"
#include "discrete/HyperGeometricRand.h"
#include "discrete/NegativeBinomialRand.h"
#include "discrete/LogarithmicRand.h"
#include "discrete/PoissonRand.h"
#include "discrete/RademacherRand.h"
#include "discrete/SkellamRand.h"
#include "discrete/UniformDiscreteRand.h"
#include "discrete/YuleRand.h"
#include "discrete/ZetaRand.h"
#include "discrete/ZipfRand.h"

/// CONTINUOUS
#include "continuous/ContinuousDistribution.h"
#include "continuous/BetaPrimeRand.h"
#include "continuous/BetaRand.h"
#include "continuous/CauchyRand.h"
#include "continuous/ExponentialRand.h"
#include "continuous/FisherSnedecorRand.h"
#include "continuous/FrechetRand.h"
#include "continuous/GammaRand.h"
#include "continuous/GeometricStableRand.h"
#include "continuous/GumbelRand.h"
#include "continuous/InverseGammaRand.h"
#include "continuous/IrwinHallRand.h"
#include "continuous/LaplaceRand.h"
#include "continuous/LevyRand.h"
#include "continuous/LogisticRand.h"
#include "continuous/LogNormalRand.h"
#include "continuous/NakagamiRand.h"
#include "continuous/NormalRand.h"
#include "continuous/MaxwellBoltzmannRand.h"
#include "continuous/ParetoRand.h"
#include "continuous/PlanckRand.h"
#include "continuous/RaisedCosineRand.h"
#include "continuous/RayleighRand.h"
#include "continuous/SechRand.h"
#include "continuous/StableRand.h"
#include "continuous/StudentTRand.h"
#include "continuous/UniformRand.h"
#include "continuous/TriangularRand.h"
#include "continuous/VonMisesRand.h"
#include "continuous/WaldRand.h"
#include "continuous/WeibullRand.h"
#include "continuous/WignerSemicircleRand.h"

/// SINGULAR
#include "singular/SingularDistribution.h"
#include "singular/CantorRand.h"

enum DISTRIBUTION {
    /// DISCRETE
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
    CANTOR_RAND
};

#endif // PROBABILITY_DISTRIBUTIONS_H
