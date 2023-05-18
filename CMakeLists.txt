cmake_minimum_required(VERSION 3.19)

project(RandLib
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(RandLib SHARED)

target_compile_features(RandLib
  PUBLIC
    cxx_std_17
)

target_include_directories(RandLib PUBLIC ${RandLib_SOURCE_DIR})

if(MSVC)
  target_compile_options(RandLib
    PRIVATE
      /MP
  )
endif()

if(UNIX)
  target_link_libraries(RandLib
    PRIVATE
      ${CMAKE_DL_LIBS}
  )
endif()

set(RANDLIB_HDRS
    # Overarching
    ${RandLib_SOURCE_DIR}/RandLib.h
    ${RandLib_SOURCE_DIR}/RandLib_global.h

    # Math
    ${RandLib_SOURCE_DIR}/math/BetaMath.h
    ${RandLib_SOURCE_DIR}/math/Constants.h
    ${RandLib_SOURCE_DIR}/math/GammaMath.h
    ${RandLib_SOURCE_DIR}/math/NumericMath.h
    ${RandLib_SOURCE_DIR}/math/RandMath.h

    # Abstract Distribution Classes
    ${RandLib_SOURCE_DIR}/distributions/Distributions.h
    ${RandLib_SOURCE_DIR}/distributions/ProbabilityDistribution.h

    # Bivariate Distributions
    ${RandLib_SOURCE_DIR}/distributions/bivariate/BivariateDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/bivariate/BivariateNormalRand.h
    ${RandLib_SOURCE_DIR}/distributions/bivariate/ContinuousBivariateDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/bivariate/DiscreteBivariateDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/bivariate/NormalInverseGammaRand.h
    ${RandLib_SOURCE_DIR}/distributions/bivariate/TrinomialRand.h

    # Univariate Distributions
    ## Helper
    ${RandLib_SOURCE_DIR}/distributions/univariate/BasicRandGenerator.h

    ## Abstract Parents
    ${RandLib_SOURCE_DIR}/distributions/univariate/ExponentialFamily.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/UnivariateDistribution.h

    ## Continuous
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/BetaPrimeRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/BetaRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/CauchyRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/CircularDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/VonMisesRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/WrappedExponentialRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ContinuousDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/DegenerateRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ExponentiallyModifiedGaussianRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ExponentialRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/FisherFRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/FrechetRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GammaRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GeometricStableRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GumbelRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/InverseGammaRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/InverseGaussianRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/IrwinHallRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/KolmogorovSmirnovRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LaplaceRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LevyRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LogisticRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LogNormalRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/MarchenkoPasturRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NakagamiRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NoncentralChiSquaredRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NormalRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ParetoRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/PlanckRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/RaisedCosineRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/SechRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/StableRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/StudentTRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/TriangularRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/UniformRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/WeibullRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/WignerSemicircleRand.h

    ## Discrete
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BernoulliRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BetaBinomialRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BinomialRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/CategoricalRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/DiscreteDistribution.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/GeometricRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/HyperGeometricRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/LogarithmicRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/NegativeBinomialRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/NegativeHyperGeometricRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/PoissonRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/RademacherRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/SkellamRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/UniformDiscreteRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/YuleRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/ZetaRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/ZipfRand.h

    ## Singular
    ${RandLib_SOURCE_DIR}/distributions/univariate/singular/CantorRand.h
    ${RandLib_SOURCE_DIR}/distributions/univariate/singular/SingularDistribution.h
)

set(RANDLIB_SRCS
    # Math
    ${RandLib_SOURCE_DIR}/math/BetaMath.cpp
    ${RandLib_SOURCE_DIR}/math/GammaMath.cpp
    ${RandLib_SOURCE_DIR}/math/NumericMath.cpp
    ${RandLib_SOURCE_DIR}/math/RandMath.cpp

    # Abstract Distribution
    ${RandLib_SOURCE_DIR}/distributions/ProbabilityDistribution.cpp

    # Bivariate Distributions
    ${RandLib_SOURCE_DIR}/distributions/bivariate/BivariateDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/bivariate/BivariateNormalRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/bivariate/ContinuousBivariateDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/bivariate/DiscreteBivariateDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/bivariate/NormalInverseGammaRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/bivariate/TrinomialRand.cpp

    # Univariate Distributions
    ## Helper
    ${RandLib_SOURCE_DIR}/distributions/univariate/BasicRandGenerator.cpp

    ## Abstract Parents
    ${RandLib_SOURCE_DIR}/distributions/univariate/ExponentialFamily.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/UnivariateDistribution.cpp

    ## Continuous
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/BetaPrimeRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/BetaRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/CauchyRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/CircularDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/VonMisesRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/circular/WrappedExponentialRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ContinuousDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/DegenerateRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ExponentiallyModifiedGaussianRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ExponentialRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/FisherFRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/FrechetRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GammaRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GeometricStableRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/GumbelRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/InverseGammaRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/InverseGaussianRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/IrwinHallRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/KolmogorovSmirnovRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LaplaceRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LevyRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LogisticRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/LogNormalRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/MarchenkoPasturRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NakagamiRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NoncentralChiSquaredRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/NormalRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/ParetoRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/PlanckRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/RaisedCosineRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/SechRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/StableRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/StudentTRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/TriangularRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/UniformRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/WeibullRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/continuous/WignerSemicircleRand.cpp

    ## Discrete
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BernoulliRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BetaBinomialRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/BinomialRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/CategoricalRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/DiscreteDistribution.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/GeometricRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/HyperGeometricRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/LogarithmicRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/NegativeBinomialRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/NegativeHyperGeometricRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/PoissonRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/RademacherRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/SkellamRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/UniformDiscreteRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/YuleRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/ZetaRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/discrete/ZipfRand.cpp

    ## Singular
    ${RandLib_SOURCE_DIR}/distributions/univariate/singular/CantorRand.cpp
    ${RandLib_SOURCE_DIR}/distributions/univariate/singular/SingularDistribution.cpp
)

target_sources(RandLib
  PRIVATE
    ${RANDLIB_HDRS}
    ${RANDLIB_SRCS}
)
