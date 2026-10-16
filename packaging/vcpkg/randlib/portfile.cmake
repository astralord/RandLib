get_filename_component(SOURCE_PATH "${CURRENT_PORT_DIR}/../../.." ABSOLUTE)
if(NOT EXISTS "${SOURCE_PATH}/include/randlib/randlib.hpp")
    message(FATAL_ERROR "This overlay port must remain inside a RandLib source checkout")
endif()

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DRANDLIB_BUILD_TESTS=OFF
        -DRANDLIB_BUILD_EXAMPLES=OFF
        -DRANDLIB_BUILD_BENCHMARKS=OFF
        -DRANDLIB_BUILD_DOCS=OFF
)
vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME RandLib CONFIG_PATH share/cmake/RandLib)

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
