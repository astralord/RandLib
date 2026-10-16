cmake_minimum_required(VERSION 3.21)
project(RandLibPackageTest LANGUAGES CXX)

find_package(RandLib CONFIG REQUIRED)

add_executable(randlib_package_test main.cpp)
target_compile_features(randlib_package_test PRIVATE cxx_std_20)
target_link_libraries(randlib_package_test PRIVATE RandLib::RandLib)
