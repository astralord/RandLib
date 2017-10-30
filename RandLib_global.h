#ifndef RANDLIB_GLOBAL_H
#define RANDLIB_GLOBAL_H

#if defined(_WIN32) || defined(_WIN64)
    #define Q_DECL_EXPORT __declspec(dllexport)
    #define Q_DECL_IMPORT __declspec(dllimport)
#else
    #define Q_DECL_EXPORT __attribute__((visibility("default")))
    #define Q_DECL_IMPORT __attribute__((visibility("default")))
#endif

#ifdef RANDLIB_LIBRARY
#define RANDLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#define RANDLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // RANDLIB_GLOBAL_H
