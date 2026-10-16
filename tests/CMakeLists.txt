find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    include(FetchContent)
    FetchContent_Declare(Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.7.1.tar.gz
        URL_HASH SHA256=c991b247a1a0d7bb9c39aa35faf0fe9e19764213f28ffba3109388e62ee0269c
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
endif()
include(Catch)

add_library(randlib_test_support STATIC
    support/data.cpp
    support/registry.cpp
    support/stats.cpp)
target_include_directories(randlib_test_support PUBLIC support)
target_link_libraries(randlib_test_support PUBLIC RandLib::RandLib Catch2::Catch2)
target_compile_definitions(randlib_test_support PUBLIC RANDLIB_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_link_libraries(randlib_test_support PUBLIC randlib_warnings)
if(MSVC)
    target_compile_options(randlib_test_support PUBLIC /bigobj)
    target_compile_definitions(randlib_test_support PUBLIC _CRT_SECURE_NO_WARNINGS)
endif()

add_executable(randlib_unit_tests
    unit/test_coherence.cpp
    unit/test_engines.cpp
    unit/test_fit_api.cpp
    unit/test_golden.cpp
    unit/test_mode_median.cpp
    unit/test_numerics.cpp
    unit/test_parameters.cpp
    unit/test_reference.cpp
    unit/test_regressions.cpp
    unit/test_stable.cpp
    unit/test_reproducibility.cpp)
target_link_libraries(randlib_unit_tests PRIVATE randlib_test_support Catch2::Catch2WithMain)
catch_discover_tests(randlib_unit_tests PROPERTIES LABELS unit)

add_executable(randlib_statistical_tests
    statistical/test_fitting.cpp
    statistical/test_sampling.cpp)
target_link_libraries(randlib_statistical_tests PRIVATE randlib_test_support Catch2::Catch2WithMain)
catch_discover_tests(randlib_statistical_tests
    TEST_SPEC "~[discrete]"
    PROPERTIES LABELS statistical TIMEOUT 3600)
catch_discover_tests(randlib_statistical_tests
    TEST_SPEC "[discrete]"
    PROPERTIES LABELS statistical TIMEOUT 60)
