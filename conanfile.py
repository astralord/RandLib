import os

from conan import ConanFile
from conan.tools.files import copy


class RandLibRecipe(ConanFile):
    name = "randlib"
    version = "2.0.0"
    package_type = "header-library"
    license = "MIT"
    author = "Aleksandr Samarin"
    url = "https://github.com/astralord/RandLib"
    homepage = "https://github.com/astralord/RandLib"
    description = "Header-only C++ probability distributions, fitting, and heavy-tail tools"
    topics = ("probability", "statistics", "random", "stable-distributions", "fitting")

    exports_sources = "include/*", "LICENSE"
    no_copy_source = True

    def package(self):
        copy(self, "*.hpp", self.source_folder, self.package_folder)
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))

    def package_id(self):
        self.info.clear()

    def package_info(self):
        self.cpp_info.bindirs = []
        self.cpp_info.libdirs = []
        self.cpp_info.set_property("cmake_file_name", "RandLib")
        self.cpp_info.set_property("cmake_target_name", "RandLib::RandLib")
