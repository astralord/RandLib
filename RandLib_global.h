#ifndef RANDLIB_GLOBAL_H
#define RANDLIB_GLOBAL_H

#if defined(RANDLIB_STATIC)
    #define RANDLIBSHARED_EXPORT
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define RANDLIB_DECL_EXPORT __declspec(dllexport)
        #define RANDLIB_DECL_IMPORT __declspec(dllimport)
    #else
        #define RANDLIB_DECL_EXPORT __attribute__((visibility("default")))
        #define RANDLIB_DECL_IMPORT __attribute__((visibility("default")))
    #endif

    #ifdef RANDLIB_LIBRARY
        #define RANDLIBSHARED_EXPORT RANDLIB_DECL_EXPORT
    #else
        #define RANDLIBSHARED_EXPORT RANDLIB_DECL_IMPORT
    #endif
#endif

#endif // RANDLIB_GLOBAL_H
