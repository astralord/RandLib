QT       -= gui

TARGET = randlib
TEMPLATE = lib

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
    continuous/ParetoRand.cpp

HEADERS +=\
    randlib_global.h \
    RandomVariable.h \
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
    continuous/ParetoRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

CONFIG += c++11
