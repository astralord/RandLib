QT       -= gui

TARGET = randlib
TEMPLATE = lib

CONFIG += c++14

DEFINES += RANDLIB_LIBRARY

SOURCES += \
    distributions/ProbabilityDistribution.cpp \
    distributions/BasicRandGenerator.cpp \
    distributions/continuous/BetaRand.cpp \
    distributions/continuous/CauchyRand.cpp \
    distributions/continuous/ExponentialRand.cpp \
    distributions/continuous/FisherSnedecorRand.cpp \
    distributions/continuous/GammaRand.cpp \
    distributions/continuous/LaplaceRand.cpp \
    distributions/continuous/LevyRand.cpp \
    distributions/continuous/LogNormalRand.cpp \
    distributions/continuous/NormalRand.cpp \
    distributions/continuous/StableRand.cpp \
    distributions/continuous/StudentTRand.cpp \
    distributions/continuous/UniformRand.cpp \
    distributions/continuous/ParetoRand.cpp \
    distributions/continuous/ContinuousDistribution.cpp \
    distributions/continuous/RayleighRand.cpp \
    distributions/continuous/TriangularRand.cpp \
    distributions/continuous/LogisticRand.cpp \
    distributions/discrete/DiscreteDistribution.cpp \
    distributions/discrete/UniformDiscreteRand.cpp \
    distributions/discrete/PoissonRand.cpp \
    distributions/continuous/NakagamiRand.cpp \
    distributions/continuous/WaldRand.cpp \
    distributions/continuous/WeibullRand.cpp \
    distributions/discrete/RademacherRand.cpp \
    math/RandMath.cpp \
    distributions/discrete/BernoulliRand.cpp \
    distributions/discrete/GeometricRand.cpp \
    distributions/discrete/BinomialRand.cpp \
    distributions/continuous/GeometricStableRand.cpp \
    distributions/continuous/MaxwellBoltzmannRand.cpp \
    distributions/continuous/BetaPrimeRand.cpp \
    distributions/discrete/NegativeBinomialRand.cpp \
    distributions/discrete/HyperGeometricRand.cpp \
    distributions/discrete/ZipfRand.cpp \
    distributions/discrete/YuleRand.cpp \
    distributions/continuous/VonMisesRand.cpp \
    distributions/continuous/ArcsineRand.cpp \
    distributions/continuous/SechRand.cpp \
    distributions/continuous/WignerSemicircleRand.cpp \
    distributions/continuous/GumbelRand.cpp \
    distributions/discrete/LogarithmicRand.cpp \
    distributions/discrete/ZetaRand.cpp \
    distributions/singular/SingularDistribution.cpp \
    distributions/singular/CantorRand.cpp \
    distributions/continuous/RaisedCosineRand.cpp \
    distributions/continuous/FrechetRand.cpp \
    distributions/discrete/SkellamRand.cpp \
    distributions/continuous/PlanckRand.cpp \
    distributions/continuous/IrwinHallRand.cpp \
    distributions/continuous/InverseGammaRand.cpp

HEADERS +=\
    randlib_global.h \
    distributions/ProbabilityDistributions.h \
    distributions/ProbabilityDistribution.h \
    distributions/BasicRandGenerator.h \
    distributions/continuous/BetaRand.h \
    distributions/continuous/CauchyRand.h \
    distributions/continuous/ExponentialRand.h \
    distributions/continuous/FisherSnedecorRand.h \
    distributions/continuous/GammaRand.h \
    distributions/continuous/LaplaceRand.h \
    distributions/continuous/LevyRand.h \
    distributions/continuous/LogNormalRand.h \
    distributions/continuous/NormalRand.h \
    distributions/continuous/StableRand.h \
    distributions/continuous/StudentTRand.h \
    distributions/continuous/UniformRand.h \
    distributions/continuous/ParetoRand.h \
    distributions/continuous/ContinuousDistribution.h \
    distributions/continuous/RayleighRand.h \
    distributions/continuous/TriangularRand.h \
    distributions/continuous/LogisticRand.h \
    distributions/discrete/DiscreteDistribution.h \
    distributions/discrete/UniformDiscreteRand.h \
    distributions/discrete/PoissonRand.h \
    distributions/continuous/NakagamiRand.h \
    distributions/continuous/WaldRand.h \
    distributions/continuous/WeibullRand.h \
    distributions/discrete/RademacherRand.h \
    math/RandMath.h \
    distributions/discrete/BernoulliRand.h \
    distributions/discrete/GeometricRand.h \
    distributions/discrete/BinomialRand.h \
    distributions/continuous/GeometricStableRand.h \
    distributions/continuous/MaxwellBoltzmannRand.h \
    distributions/continuous/BetaPrimeRand.h \
    distributions/discrete/NegativeBinomialRand.h \
    distributions/discrete/HyperGeometricRand.h \
    distributions/discrete/ZipfRand.h \
    distributions/discrete/YuleRand.h \
    distributions/continuous/VonMisesRand.h \
    distributions/continuous/ArcsineRand.h \
    math/Constants.h \
    distributions/continuous/SechRand.h \
    distributions/continuous/WignerSemicircleRand.h \
    distributions/continuous/GumbelRand.h \
    distributions/discrete/LogarithmicRand.h \
    distributions/discrete/ZetaRand.h \
    distributions/singular/SingularDistribution.h \
    distributions/singular/CantorRand.h \
    distributions/continuous/RaisedCosineRand.h \
    distributions/continuous/FrechetRand.h \
    distributions/discrete/SkellamRand.h \
    distributions/continuous/PlanckRand.h \
    distributions/continuous/IrwinHallRand.h \
    distributions/continuous/InverseGammaRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

