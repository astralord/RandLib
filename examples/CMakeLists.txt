find_package(Threads REQUIRED)

set(RANDLIB_EXAMPLES
    sampling_std_engine
    parallel_streams
    fit_model_selection
    bayesian_update)

foreach(example IN LISTS RANDLIB_EXAMPLES)
    add_executable(randlib_example_${example} ${example}.cpp)
    target_link_libraries(randlib_example_${example} PRIVATE RandLib::RandLib randlib_warnings Threads::Threads)
    if(RANDLIB_BUILD_TESTS)
        add_test(NAME example_${example} COMMAND randlib_example_${example})
        set_tests_properties(example_${example} PROPERTIES LABELS "example;unit")
    endif()
endforeach()
