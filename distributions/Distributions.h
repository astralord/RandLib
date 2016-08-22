#ifndef RANDLIB_H
#define RANDLIB_H

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
#include "univariate/continuous/DegenerateRand.h"
#include "univariate/continuous/ExponentialRand.h"
#include "univariate/continuous/ExponentiallyModifiedGaussianRand.h"
#include "univariate/continuous/FisherSnedecorRand.h"
#include "univariate/continuous/FrechetRand.h"
#include "univariate/continuous/GammaRand.h"
#include "univariate/continuous/GeometricStableRand.h"
#include "univariate/continuous/GumbelRand.h"
#include "univariate/continuous/InverseGammaRand.h"
#include "univariate/continuous/IrwinHallRand.h"
#include "univariate/continuous/LaplaceRand.h"
#include "univariate/continuous/LevyRand.h"
#include "univariate/continuous/LimitingDistribution.h"
#include "univariate/continuous/LogCauchyRand.h"
#include "univariate/continuous/LogisticRand.h"
#include "univariate/continuous/LogNormalRand.h"
#include "univariate/continuous/NakagamiRand.h"
#include "univariate/continuous/NoncentralChiSquared.h"
#include "univariate/continuous/NormalRand.h"
#include "univariate/continuous/ParetoRand.h"
#include "univariate/continuous/PlanckRand.h"
#include "univariate/continuous/RaisedCosineRand.h"
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
#include "multivariate/BivariateNormalRand.h"

/// SUM
#include "univariate/SumRand.h"

#endif // RANDLIB_H
