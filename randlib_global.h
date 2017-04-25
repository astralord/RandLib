#ifndef RANDLIB_GLOBAL_H
#define RANDLIB_GLOBAL_H

#include <algorithm>

#ifdef _WIN32
    #define Q_DECL_EXPORT __declspec(dllexport)
    #define Q_DECL_IMPORT __declspec(dllimport)
#else
    #define Q_DECL_EXPORT
    #define Q_DECL_IMPORT
#endif

#ifdef RANDLIB_LIBRARY
#define RANDLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#define RANDLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // RANDLIB_GLOBAL_H
