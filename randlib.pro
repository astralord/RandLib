QT       -= gui

TARGET = randlib
TEMPLATE = lib

DEFINES += RANDLIB_LIBRARY

SOURCES += \
    RandomVariable.cpp \
    UniformRand.cpp \
    NormalRand.cpp \
    CauchyRand.cpp \
    LevyRand.cpp \
    ExponentialRand.cpp \
    LaplaceRand.cpp \
    LogNormalRand.cpp \
    ChiSquaredRand.cpp \
    StudentTRand.cpp \
    BetaRand.cpp \
    GammaRand.cpp \
    StableRand.cpp \
    FisherSnedecorRand.cpp

HEADERS +=\
    randlib_global.h \
    RandomVariable.h \
    UniformRand.h \
    NormalRand.h \
    CauchyRand.h \
    LevyRand.h \
    ExponentialRand.h \
    LaplaceRand.h \
    LogNormalRand.h \
    ChiSquaredRand.h \
    RandomGenerators.h \
    StudentTRand.h \
    BetaRand.h \
    GammaRand.h \
    StableRand.h \
    FisherSnedecorRand.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

CONFIG += c++11
