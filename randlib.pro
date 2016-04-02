QT       -= gui

TARGET = randlib
TEMPLATE = lib

CONFIG += c++14

DEFINES += RANDLIB_LIBRARY

SOURCES += \
    distributions/ProbabilityDistribution.cpp \
    distributions/univariate/BasicRandGenerator.cpp \
    distributions/univariate/continuous/BetaRand.cpp \
    distributions/univariate/continuous/CauchyRand.cpp \
    distributions/univariate/continuous/ExponentialRand.cpp \
    distributions/univariate/continuous/FisherSnedecorRand.cpp \
    distributions/univariate/continuous/GammaRand.cpp \
    distributions/univariate/continuous/LaplaceRand.cpp \
    distributions/univariate/continuous/LevyRand.cpp \
    distributions/univariate/continuous/LogNormalRand.cpp \
    distributions/univariate/continuous/NormalRand.cpp \
    distributions/univariate/continuous/StableRand.cpp \
    distributions/univariate/continuous/StudentTRand.cpp \
    distributions/univariate/continuous/UniformRand.cpp \
    distributions/univariate/continuous/ParetoRand.cpp \
    distributions/univariate/continuous/ContinuousDistribution.cpp \
    distributions/univariate/continuous/RayleighRand.cpp \
    distributions/univariate/continuous/TriangularRand.cpp \
    distributions/univariate/continuous/LogisticRand.cpp \
    distributions/univariate/discrete/DiscreteDistribution.cpp \
    distributions/univariate/discrete/UniformDiscreteRand.cpp \
    distributions/univariate/discrete/PoissonRand.cpp \
    distributions/univariate/continuous/NakagamiRand.cpp \
    distributions/univariate/continuous/WaldRand.cpp \
    distributions/univariate/continuous/WeibullRand.cpp \
    distributions/univariate/discrete/RademacherRand.cpp \
    math/RandMath.cpp \
    distributions/univariate/discrete/BernoulliRand.cpp \
    distributions/univariate/discrete/GeometricRand.cpp \
    distributions/univariate/discrete/BinomialRand.cpp \
    distributions/univariate/continuous/GeometricStableRand.cpp \
    distributions/univariate/continuous/MaxwellBoltzmannRand.cpp \
    distributions/univariate/continuous/BetaPrimeRand.cpp \
    distributions/univariate/discrete/NegativeBinomialRand.cpp \
    distributions/univariate/discrete/HyperGeometricRand.cpp \
    distributions/univariate/discrete/ZipfRand.cpp \
    distributions/univariate/discrete/YuleRand.cpp \
    distributions/univariate/continuous/VonMisesRand.cpp \
    distributions/univariate/continuous/SechRand.cpp \
    distributions/univariate/continuous/WignerSemicircleRand.cpp \
    distributions/univariate/continuous/GumbelRand.cpp \
    distributions/univariate/discrete/LogarithmicRand.cpp \
    distributions/univariate/discrete/ZetaRand.cpp \
    distributions/univariate/singular/SingularDistribution.cpp \
    distributions/univariate/singular/CantorRand.cpp \
    distributions/univariate/continuous/RaisedCosineRand.cpp \
    distributions/univariate/continuous/FrechetRand.cpp \
    distributions/univariate/discrete/SkellamRand.cpp \
    distributions/univariate/continuous/PlanckRand.cpp \
    distributions/univariate/continuous/IrwinHallRand.cpp \
    distributions/univariate/continuous/InverseGammaRand.cpp \
    distributions/univariate/UnivariateProbabilityDistribution.cpp \
    distributions/multivariate/BivariateProbabilityDistribution.cpp \
    distributions/multivariate/NormalInverseGammaRand.cpp \
    distributions/univariate/discrete/BetaBinomialRand.cpp \
    distributions/univariate/continuous/ExponentialNormalRand.cpp \
    math/Matrix.cpp \
    distributions/multivariate/MultivariateProbabilityDistribution.cpp \
    distributions/univariate/continuous/LogCauchyRand.cpp

HEADERS +=\
    randlib_global.h \
    distributions/ProbabilityDistributions.h \
    distributions/ProbabilityDistribution.h \
    distributions/univariate/BasicRandGenerator.h \
    distributions/univariate/continuous/BetaRand.h \
    distributions/univariate/continuous/CauchyRand.h \
    distributions/univariate/continuous/ExponentialRand.h \
    distributions/univariate/continuous/FisherSnedecorRand.h \
    distributions/univariate/continuous/GammaRand.h \
    distributions/univariate/continuous/LaplaceRand.h \
    distributions/univariate/continuous/LevyRand.h \
    distributions/univariate/continuous/LogNormalRand.h \
    distributions/univariate/continuous/NormalRand.h \
    distributions/univariate/continuous/StableRand.h \
    distributions/univariate/continuous/StudentTRand.h \
    distributions/univariate/continuous/UniformRand.h \
    distributions/univariate/continuous/ParetoRand.h \
    distributions/univariate/continuous/ContinuousDistribution.h \
    distributions/univariate/continuous/RayleighRand.h \
    distributions/univariate/continuous/TriangularRand.h \
    distributions/univariate/continuous/LogisticRand.h \
    distributions/univariate/discrete/DiscreteDistribution.h \
    distributions/univariate/discrete/UniformDiscreteRand.h \
    distributions/univariate/discrete/PoissonRand.h \
    distributions/univariate/continuous/NakagamiRand.h \
    distributions/univariate/continuous/WaldRand.h \
    distributions/univariate/continuous/WeibullRand.h \
    distributions/univariate/discrete/RademacherRand.h \
    math/RandMath.h \
    distributions/univariate/discrete/BernoulliRand.h \
    distributions/univariate/discrete/GeometricRand.h \
    distributions/univariate/discrete/BinomialRand.h \
    distributions/univariate/continuous/GeometricStableRand.h \
    distributions/univariate/continuous/MaxwellBoltzmannRand.h \
    distributions/univariate/continuous/BetaPrimeRand.h \
    distributions/univariate/discrete/NegativeBinomialRand.h \
    distributions/univariate/discrete/HyperGeometricRand.h \
    distributions/univariate/discrete/ZipfRand.h \
    distributions/univariate/discrete/YuleRand.h \
    distributions/univariate/continuous/VonMisesRand.h \
    math/Constants.h \
    distributions/univariate/continuous/SechRand.h \
    distributions/univariate/continuous/WignerSemicircleRand.h \
    distributions/univariate/continuous/GumbelRand.h \
    distributions/univariate/discrete/LogarithmicRand.h \
    distributions/univariate/discrete/ZetaRand.h \
    distributions/univariate/singular/SingularDistribution.h \
    distributions/univariate/singular/CantorRand.h \
    distributions/univariate/continuous/RaisedCosineRand.h \
    distributions/univariate/continuous/FrechetRand.h \
    distributions/univariate/discrete/SkellamRand.h \
    distributions/univariate/continuous/PlanckRand.h \
    distributions/univariate/continuous/IrwinHallRand.h \
    distributions/univariate/continuous/InverseGammaRand.h \
    distributions/univariate/UnivariateProbabilityDistribution.h \
    distributions/multivariate/BivariateProbabilityDistribution.h \
    distributions/multivariate/NormalInverseGammaRand.h \
    distributions/univariate/discrete/BetaBinomialRand.h \
    distributions/univariate/continuous/ExponentialNormalRand.h \
    math/Matrix.h \
    distributions/multivariate/MultivariateProbabilityDistribution.h \
    distributions/univariate/continuous/LogCauchyRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

