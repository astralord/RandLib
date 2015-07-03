QT       -= gui

TARGET = randlib
TEMPLATE = lib

CONFIG += c++11

DEFINES += RANDLIB_LIBRARY

SOURCES += \
    RandomVariable.cpp \
    continuous/BetaRand.cpp \
    continuous/CauchyRand.cpp \
    continuous/ChiSquaredRand.cpp \
    continuous/ExponentialRand.cpp \
    continuous/FisherSnedecorRand.cpp \
    continuous/GammaRand.cpp \
    continuous/LaplaceRand.cpp \
    continuous/LevyRand.cpp \
    continuous/LogNormalRand.cpp \
    continuous/NormalRand.cpp \
    continuous/StableRand.cpp \
    continuous/StudentTRand.cpp \
    continuous/UniformRand.cpp \
    continuous/ParetoRand.cpp \
    discrete/PoissonRand.cpp \
    RandMath.cpp \
    continuous/ErlangRand.cpp \
    continuous/ContinuousRand.cpp \
    discrete/DiscreteRand.cpp \
    continuous/RayleighRand.cpp \
    discrete/UniformDiscreteRand.cpp \
    BasicRandGenerator.cpp \
    continuous/TriangularRand.cpp

HEADERS +=\
    randlib_global.h \
    RandomGenerators.h \
    continuous/BetaRand.h \
    continuous/CauchyRand.h \
    continuous/ChiSquaredRand.h \
    continuous/ExponentialRand.h \
    continuous/FisherSnedecorRand.h \
    continuous/GammaRand.h \
    continuous/LaplaceRand.h \
    continuous/LevyRand.h \
    continuous/LogNormalRand.h \
    continuous/NormalRand.h \
    continuous/StableRand.h \
    continuous/StudentTRand.h \
    continuous/UniformRand.h \
    continuous/ParetoRand.h \
    discrete/PoissonRand.h \
    RandMath.h \
    RandomVariable.h \
    continuous/ErlangRand.h \
    continuous/ContinuousRand.h \
    discrete/DiscreteRand.h \
    continuous/RayleighRand.h \
    discrete/UniformDiscreteRand.h \
    BasicRandGenerator.h \
    continuous/TriangularRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

