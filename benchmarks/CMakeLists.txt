find_package(benchmark 1.9 QUIET)
if(NOT benchmark_FOUND)
    include(FetchContent)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(google_benchmark
        URL https://github.com/google/benchmark/archive/refs/tags/v1.9.5.tar.gz
        URL_HASH SHA256=9631341c82bac4a288bef951f8b26b41f69021794184ece969f8473977eaa340
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(google_benchmark)
endif()

add_executable(randlib_benchmarks sampling.cpp)
target_link_libraries(randlib_benchmarks PRIVATE RandLib::RandLib benchmark::benchmark_main randlib_warnings)
