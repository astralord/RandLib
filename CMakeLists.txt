cmake_minimum_required(VERSION 3.21)

project(RandLib
    VERSION 2.0.0
    DESCRIPTION "Probability distributions: sampling, properties and fitting"
    HOMEPAGE_URL "https://github.com/astralord/RandLib"
    LANGUAGES CXX)

include(GNUInstallDirs)

option(RANDLIB_BUILD_TESTS "Build the RandLib test suite" ${PROJECT_IS_TOP_LEVEL})
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
    target_compile_options(randlib_warnings INTERFACE /W4 /permissive- /bigobj)
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
