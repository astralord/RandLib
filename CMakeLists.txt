cmake_minimum_required(VERSION 3.21)

project(RandLib
    VERSION 2.0.0
    DESCRIPTION "Probability distributions: sampling, properties and fitting"
    HOMEPAGE_URL "https://github.com/astralord/RandLib"
    LANGUAGES CXX)

include(GNUInstallDirs)

option(RANDLIB_BUILD_TESTS "Build the RandLib test suite" ${PROJECT_IS_TOP_LEVEL})
option(RANDLIB_BUILD_EXAMPLES "Build the RandLib examples" ${PROJECT_IS_TOP_LEVEL})
option(RANDLIB_BUILD_BENCHMARKS "Build the RandLib sampling benchmarks" OFF)
option(RANDLIB_BUILD_DOCS "Build the Doxygen API reference" OFF)
option(RANDLIB_WARNINGS_AS_ERRORS "Treat compiler warnings as errors in RandLib's own targets" OFF)
option(RANDLIB_INSTALL "Generate install and package-config rules" ${PROJECT_IS_TOP_LEVEL})
option(RANDLIB_PRECOMPILED "Compile the double and int instantiations once into a static library" OFF)
option(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS
    "Use RandLib's own Bessel, beta and zeta functions even when the standard library provides them" OFF)

add_library(RandLib INTERFACE)
add_library(RandLib::RandLib ALIAS RandLib)
target_compile_features(RandLib INTERFACE cxx_std_20)
target_include_directories(RandLib INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
    target_compile_definitions(RandLib INTERFACE RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
endif()

if(RANDLIB_PRECOMPILED)
    add_library(RandLibPrecompiled STATIC src/precompiled.cpp)
    target_include_directories(RandLibPrecompiled PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
    target_compile_features(RandLibPrecompiled PUBLIC cxx_std_20)
    target_compile_definitions(RandLibPrecompiled PUBLIC RANDLIB_PRECOMPILED)
    if(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
        target_compile_definitions(RandLibPrecompiled PUBLIC RANDLIB_PORTABLE_SPECIAL_FUNCTIONS)
    endif()
    if(MSVC)
        target_compile_options(RandLibPrecompiled PRIVATE /bigobj)
    endif()
    target_link_libraries(RandLib INTERFACE RandLibPrecompiled)
endif()

# Warning flags for the targets built in this project (tests, examples, benchmarks)
add_library(randlib_warnings INTERFACE)
if(MSVC)
    # implicit narrowing conversions (C4244, C4267, C4305) are not warned about by GCC and Clang either
    target_compile_options(randlib_warnings INTERFACE /W4 /permissive- /bigobj /wd4244 /wd4267 /wd4305)
    if(RANDLIB_WARNINGS_AS_ERRORS)
        target_compile_options(randlib_warnings INTERFACE /WX)
    endif()
else()
    target_compile_options(randlib_warnings INTERFACE -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic)
    if(RANDLIB_WARNINGS_AS_ERRORS)
        target_compile_options(randlib_warnings INTERFACE -Werror)
    endif()
endif()

if(RANDLIB_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(RANDLIB_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(RANDLIB_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(RANDLIB_BUILD_DOCS)
    find_package(Doxygen REQUIRED)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    include(FetchContent)
    FetchContent_Declare(doxygen_awesome
        URL https://github.com/jothepro/doxygen-awesome-css/archive/d52eafe3e9303399fda15661f3d7bb8fe3d7eabc.tar.gz
        URL_HASH SHA256=6d58bfc1546cab9ce606494a3df857271b9cc712db9b9b77a9abf0567d61c17d
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(doxygen_awesome)
    FetchContent_GetProperties(doxygen_awesome SOURCE_DIR RANDLIB_DOXYGEN_AWESOME_DIR)
    configure_file(docs/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(randlib_docs ALL
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tools/docs/check_reference_coverage.py
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        COMMENT "Generating RandLib API documentation"
        VERBATIM)
endif()

if(RANDLIB_INSTALL)
    include(CMakePackageConfigHelpers)

    install(TARGETS RandLib EXPORT RandLibTargets)
    if(RANDLIB_PRECOMPILED)
        install(TARGETS RandLibPrecompiled EXPORT RandLibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
    endif()
    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    if(RANDLIB_PRECOMPILED)
        set(RANDLIB_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/RandLib)
    else()
        set(RANDLIB_CMAKE_DIR ${CMAKE_INSTALL_DATADIR}/cmake/RandLib)
    endif()
    install(EXPORT RandLibTargets NAMESPACE RandLib:: DESTINATION ${RANDLIB_CMAKE_DIR})
    configure_package_config_file(cmake/RandLibConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfig.cmake
        INSTALL_DESTINATION ${RANDLIB_CMAKE_DIR})
    if(RANDLIB_PRECOMPILED)
        write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RandLibConfigVersion.cmake
            COMPATIBILITY SameMajorVersion)
    else()
        write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/RandLibConfigVersion.cmake
            COMPATIBILITY SameMajorVersion
            ARCH_INDEPENDENT)
    endif()
    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/RandLibConfigVersion.cmake
        DESTINATION ${RANDLIB_CMAKE_DIR})
    export(EXPORT RandLibTargets NAMESPACE RandLib:: FILE ${CMAKE_CURRENT_BINARY_DIR}/RandLibTargets.cmake)
endif()
