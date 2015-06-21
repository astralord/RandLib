#ifndef RANDLIB_GLOBAL_H
#define RANDLIB_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(RANDLIB_LIBRARY)
#  define RANDLIBSHARED_EXPORT Q_DECL_EXPORT
#else
#  define RANDLIBSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // RANDLIB_GLOBAL_H
