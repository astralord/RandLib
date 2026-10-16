# Packaging and the single header

The CMake install exports `RandLib::RandLib`. In addition to direct CMake consumption, the source
tree contains locally testable vcpkg and Conan packages.

## vcpkg overlay

The overlay intentionally packages the current checkout, which keeps pre-release validation
independent of a published archive:

```sh
vcpkg install randlib --overlay-ports=packaging/vcpkg
```

The upstream vcpkg submission will replace the local source path with the final release archive and
its SHA-512 checksum.

## Conan 2

```sh
conan profile detect --force
conan create . --build=missing
```

The recipe is a platform-independent header package and the `test_package/` consumer verifies the
generated `RandLib::RandLib` CMake target.

## Single header

Generate the release artifact deterministically from the public umbrella header:

```sh
python3 tools/amalgamate.py --output randlib_single.hpp
```

The generator recursively expands internal `randlib/...` includes, preserves standard-library
includes and preprocessor conditions, and emits each project header once. Tagged v2 releases compile
the result, attach it to the GitHub Release, and publish a SHA-256 checksum.
