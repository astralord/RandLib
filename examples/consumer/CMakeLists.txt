cmake_minimum_required(VERSION 3.21)
project(RandLibConsumer LANGUAGES CXX)

find_package(RandLib 1.2 REQUIRED)

add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE RandLib::RandLib)
