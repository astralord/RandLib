QT       -= gui

TARGET = randlib
TEMPLATE = lib

CONFIG += c++11

DEFINES += RANDLIB_LIBRARY

SOURCES += \
    RandMath.cpp \
    variables/RandomVariable.cpp \
    variables/BasicRandGenerator.cpp \
    variables/continuous/BetaRand.cpp \
    variables/continuous/CauchyRand.cpp \
    variables/continuous/ChiSquaredRand.cpp \
    variables/continuous/ExponentialRand.cpp \
    variables/continuous/FisherSnedecorRand.cpp \
    variables/continuous/GammaRand.cpp \
    variables/continuous/LaplaceRand.cpp \
    variables/continuous/LevyRand.cpp \
    variables/continuous/LogNormalRand.cpp \
    variables/continuous/NormalRand.cpp \
    variables/continuous/StableRand.cpp \
    variables/continuous/StudentTRand.cpp \
    variables/continuous/UniformRand.cpp \
    variables/continuous/ParetoRand.cpp \
    variables/continuous/ErlangRand.cpp \
    variables/continuous/ContinuousRand.cpp \
    variables/continuous/RayleighRand.cpp \
    variables/continuous/TriangularRand.cpp \
    variables/continuous/LogisticRand.cpp \
    variables/discrete/DiscreteRand.cpp \
    variables/discrete/UniformDiscreteRand.cpp \
    variables/discrete/PoissonRand.cpp \
    processes/StochasticProcess.cpp \
    processes/WienerProcess.cpp \
    variables/continuous/NakagamiRand.cpp \
    variables/continuous/WaldRand.cpp

HEADERS +=\
    randlib_global.h \
    RandMath.h \
    variables/RandomGenerators.h \
    variables/RandomVariable.h \
    variables/BasicRandGenerator.h \
    variables/continuous/BetaRand.h \
    variables/continuous/CauchyRand.h \
    variables/continuous/ChiSquaredRand.h \
    variables/continuous/ExponentialRand.h \
    variables/continuous/FisherSnedecorRand.h \
    variables/continuous/GammaRand.h \
    variables/continuous/LaplaceRand.h \
    variables/continuous/LevyRand.h \
    variables/continuous/LogNormalRand.h \
    variables/continuous/NormalRand.h \
    variables/continuous/StableRand.h \
    variables/continuous/StudentTRand.h \
    variables/continuous/UniformRand.h \
    variables/continuous/ParetoRand.h \
    variables/continuous/ErlangRand.h \
    variables/continuous/ContinuousRand.h \
    variables/continuous/RayleighRand.h \
    variables/continuous/TriangularRand.h \
    variables/continuous/LogisticRand.h \
    variables/discrete/DiscreteRand.h \
    variables/discrete/UniformDiscreteRand.h \
    variables/discrete/PoissonRand.h \
    processes/StochasticProcess.h \
    processes/WienerProcess.h \
    variables/continuous/NakagamiRand.h \
    variables/continuous/WaldRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

