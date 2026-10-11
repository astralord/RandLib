cmake_minimum_required(VERSION 3.21)

project(RandLib
    VERSION 1.2.0
    DESCRIPTION "Probability distributions: sampling, properties and fitting"
    HOMEPAGE_URL "https://github.com/astralord/RandLib"
    LANGUAGES CXX)

include(GNUInstallDirs)

option(RANDLIB_BUILD_TESTS "Build the RandLib test suite" ${PROJECT_IS_TOP_LEVEL})
option(RANDLIB_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(RANDLIB_INSTALL "Generate install and package-config rules" ${PROJECT_IS_TOP_LEVEL})
set(RANDLIB_GENERATOR "JKISS" CACHE STRING "Engine behind RandGenerator: JKISS, JLKISS64 or PCG32")
set_property(CACHE RANDLIB_GENERATOR PROPERTY STRINGS JKISS JLKISS64 PCG32)
option(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS
    "Use RandLib's own Bessel, beta and zeta functions even when the standard library provides them" OFF)

add_library(RandLib
    distributions/ProbabilityDistribution.cpp
    distributions/bivariate/BivariateDistribution.cpp
    distributions/bivariate/BivariateNormalRand.cpp
    distributions/bivariate/ContinuousBivariateDistribution.cpp
    distributions/bivariate/DiscreteBivariateDistribution.cpp
    distributions/bivariate/NormalInverseGammaRand.cpp
    distributions/bivariate/TrinomialRand.cpp
    distributions/univariate/BasicRandGenerator.cpp
    distributions/univariate/ExponentialFamily.cpp
    distributions/univariate/UnivariateDistribution.cpp
    distributions/univariate/continuous/BetaPrimeRand.cpp
    distributions/univariate/continuous/BetaRand.cpp
    distributions/univariate/continuous/CauchyRand.cpp
    distributions/univariate/continuous/ContinuousDistribution.cpp
    distributions/univariate/continuous/DegenerateRand.cpp
    distributions/univariate/continuous/ExponentialRand.cpp
    distributions/univariate/continuous/ExponentiallyModifiedGaussianRand.cpp
    distributions/univariate/continuous/FisherFRand.cpp
    distributions/univariate/continuous/FrechetRand.cpp
    distributions/univariate/continuous/GammaRand.cpp
    distributions/univariate/continuous/GeometricStableRand.cpp
    distributions/univariate/continuous/GumbelRand.cpp
    distributions/univariate/continuous/InverseGammaRand.cpp
    distributions/univariate/continuous/InverseGaussianRand.cpp
    distributions/univariate/continuous/IrwinHallRand.cpp
    distributions/univariate/continuous/KolmogorovSmirnovRand.cpp
    distributions/univariate/continuous/LaplaceRand.cpp
    distributions/univariate/continuous/LevyRand.cpp
    distributions/univariate/continuous/LogNormalRand.cpp
    distributions/univariate/continuous/LogisticRand.cpp
    distributions/univariate/continuous/MarchenkoPasturRand.cpp
    distributions/univariate/continuous/NakagamiRand.cpp
    distributions/univariate/continuous/NoncentralChiSquaredRand.cpp
    distributions/univariate/continuous/NormalRand.cpp
    distributions/univariate/continuous/ParetoRand.cpp
    distributions/univariate/continuous/PlanckRand.cpp
    distributions/univariate/continuous/RaisedCosineRand.cpp
    distributions/univariate/continuous/SechRand.cpp
    distributions/univariate/continuous/StableRand.cpp
    distributions/univariate/continuous/StudentTRand.cpp
    distributions/univariate/continuous/TriangularRand.cpp
    distributions/univariate/continuous/UniformRand.cpp
    distributions/univariate/continuous/WeibullRand.cpp
    distributions/univariate/continuous/WignerSemicircleRand.cpp
    distributions/univariate/continuous/circular/CircularDistribution.cpp
    distributions/univariate/continuous/circular/VonMisesRand.cpp
    distributions/univariate/continuous/circular/WrappedExponentialRand.cpp
    distributions/univariate/discrete/BernoulliRand.cpp
    distributions/univariate/discrete/BetaBinomialRand.cpp
    distributions/univariate/discrete/BinomialRand.cpp
    distributions/univariate/discrete/CategoricalRand.cpp
    distributions/univariate/discrete/DiscreteDistribution.cpp
    distributions/univariate/discrete/GeometricRand.cpp
    distributions/univariate/discrete/HyperGeometricRand.cpp
    distributions/univariate/discrete/LogarithmicRand.cpp
    distributions/univariate/discrete/NegativeBinomialRand.cpp
    distributions/univariate/discrete/NegativeHyperGeometricRand.cpp
    distributions/univariate/discrete/PoissonRand.cpp
    distributions/univariate/discrete/RademacherRand.cpp
    distributions/univariate/discrete/SkellamRand.cpp
    distributions/univariate/discrete/UniformDiscreteRand.cpp
    distributions/univariate/discrete/YuleRand.cpp
    distributions/univariate/discrete/ZetaRand.cpp
    distributions/univariate/discrete/ZipfRand.cpp
    distributions/univariate/singular/CantorRand.cpp
    distributions/univariate/singular/SingularDistribution.cpp
    math/BetaMath.cpp
    math/GammaMath.cpp
    math/NumericMath.cpp
    math/RandMath.cpp
    math/SpecialFunctions.cpp
)
add_library(RandLib::RandLib ALIAS RandLib)

target_compile_features(RandLib PUBLIC cxx_std_17)
set_target_properties(RandLib PROPERTIES
    CXX_EXTENSIONS OFF
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
target_include_directories(RandLib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/RandLib>)
target_compile_definitions(RandLib PRIVATE RANDLIB_LIBRARY)
if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(RandLib PUBLIC RANDLIB_STATIC)
endif()
if(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
    target_compile_definitions(RandLib PUBLIC RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
endif()
if(RANDLIB_GENERATOR STREQUAL "JLKISS64")
    target_compile_definitions(RandLib PUBLIC RANDLIB_JLKISS64)
elseif(RANDLIB_GENERATOR STREQUAL "PCG32")
    target_compile_definitions(RandLib PUBLIC RANDLIB_PCG32)
elseif(NOT RANDLIB_GENERATOR STREQUAL "JKISS")
    message(FATAL_ERROR "RANDLIB_GENERATOR must be JKISS, JLKISS64 or PCG32, not ${RANDLIB_GENERATOR}")
endif()

if(MSVC)
    # implicit narrowing conversions (C4244, C4267, C4305) are not warned about by GCC and Clang either
    target_compile_options(RandLib PRIVATE /W4 /permissive- /bigobj /wd4244 /wd4267 /wd4305)
    if(RANDLIB_WARNINGS_AS_ERRORS)
        target_compile_options(RandLib PRIVATE /WX)
    endif()
else()
    target_compile_options(RandLib PRIVATE -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic)
    if(RANDLIB_WARNINGS_AS_ERRORS)
        target_compile_options(RandLib PRIVATE -Werror)
    endif()
endif()

if(RANDLIB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(RANDLIB_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS RandLib
        EXPORT RandLibTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES RandLib.h RandLib_global.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/RandLib)
    install(DIRECTORY distributions math
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/RandLib
        FILES_MATCHING PATTERN "*.h")

    set(RANDLIB_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/RandLib)
    install(EXPORT RandLibTargets NAMESPACE RandLib:: DESTINATION ${RANDLIB_CMAKE_DIR})
    configure_package_config_file(cmake/RandLibConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfig.cmake
        INSTALL_DESTINATION ${RANDLIB_CMAKE_DIR})
    write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RandLibConfigVersion.cmake
        COMPATIBILITY SameMajorVersion)
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfigVersion.cmake
        DESTINATION ${RANDLIB_CMAKE_DIR})
    export(EXPORT RandLibTargets NAMESPACE RandLib:: FILE ${CMAKE_CURRENT_BINARY_DIR}/RandLibTargets.cmake)
endif()
