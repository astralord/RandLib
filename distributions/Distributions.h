#ifndef RANDLIB_H
#define RANDLIB_H

#include "ProbabilityDistribution.h"
#include "univariate/BasicRandGenerator.h"

/// UNIVARIATE
#include "univariate/UnivariateDistribution.h"

/// CONTINUOUS
#include "univariate/continuous/ContinuousDistribution.h"
#include "univariate/continuous/BetaPrimeRand.h"
#include "univariate/continuous/BetaRand.h"
#include "univariate/continuous/CauchyRand.h"
#include "univariate/continuous/DegenerateRand.h"
#include "univariate/continuous/ExponentialRand.h"
#include "univariate/continuous/ExponentiallyModifiedGaussianRand.h"
#include "univariate/continuous/FisherFRand.h"
#include "univariate/continuous/FrechetRand.h"
#include "univariate/continuous/GammaRand.h"
#include "univariate/continuous/GeometricStableRand.h"
#include "univariate/continuous/GumbelRand.h"
#include "univariate/continuous/InverseGammaRand.h"
#include "univariate/continuous/InverseGaussianRand.h"
#include "univariate/continuous/IrwinHallRand.h"
#include "univariate/continuous/KolmogorovSmirnovRand.h"
#include "univariate/continuous/LaplaceRand.h"
#include "univariate/continuous/LevyRand.h"
#include "univariate/continuous/LogisticRand.h"
#include "univariate/continuous/LogNormalRand.h"
#include "univariate/continuous/MarchenkoPasturRand.h"
#include "univariate/continuous/NakagamiRand.h"
#include "univariate/continuous/NoncentralChiSquaredRand.h"
#include "univariate/continuous/NormalRand.h"
#include "univariate/continuous/ParetoRand.h"
#include "univariate/continuous/PlanckRand.h"
#include "univariate/continuous/RaisedCosineRand.h"
#include "univariate/continuous/SechRand.h"
#include "univariate/continuous/StableRand.h"
#include "univariate/continuous/StudentTRand.h"
#include "univariate/continuous/UniformRand.h"
#include "univariate/continuous/TriangularRand.h"
#include "univariate/continuous/WeibullRand.h"
#include "univariate/continuous/WignerSemicircleRand.h"

/// CIRCULAR
#include "univariate/continuous/circular/VonMisesRand.h"
#include "univariate/continuous/circular/WrappedExponentialRand.h"

/// DISCRETE
#include "univariate/discrete/DiscreteDistribution.h"
#include "univariate/discrete/BernoulliRand.h"
#include "univariate/discrete/BetaBinomialRand.h"
#include "univariate/discrete/BinomialRand.h"
#include "univariate/discrete/CategoricalRand.h"
#include "univariate/discrete/GeometricRand.h"
#include "univariate/discrete/HyperGeometricRand.h"
#include "univariate/discrete/NegativeBinomialRand.h"
#include "univariate/discrete/NegativeHyperGeometricRand.h"
#include "univariate/discrete/LogarithmicRand.h"
#include "univariate/discrete/PoissonRand.h"
#include "univariate/discrete/RademacherRand.h"
#include "univariate/discrete/SkellamRand.h"
#include "univariate/discrete/UniformDiscreteRand.h"
#include "univariate/discrete/YuleRand.h"
#include "univariate/discrete/ZetaRand.h"
#include "univariate/discrete/ZipfRand.h"

/// SINGULAR
#include "univariate/singular/SingularDistribution.h"
#include "univariate/singular/CantorRand.h"

/// BIVARIATE
#include "bivariate/ContinuousBivariateDistribution.h"
#include "bivariate/NormalInverseGammaRand.h"
#include "bivariate/BivariateNormalRand.h"
#include "bivariate/TrinomialRand.h"

#endif // RANDLIB_H
