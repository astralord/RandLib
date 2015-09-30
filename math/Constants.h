#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>

#ifndef INFINITY
#include <limits>
static long double INFINITY = std::numeric_limits<long double>::infinity();
#endif

#ifndef NAN
#include <limits>
static long double NAN = std::numeric_limits<long double>::quiet_NaN();
#endif

#ifndef M_E
constexpr long double M_E         = 2.71828182845904523536;
#endif

#ifndef M_LOG2E
constexpr long double M_LOG2E     = 1.44269504088896340760;
#endif

#ifndef M_LOG10E
constexpr long double M_LOG10E    = 0.43429448190325182765;
#endif

#ifndef M_LN2
constexpr long double M_LN2       = 0.69314718055994530942;
#endif

#ifndef M_LN10
constexpr long double M_LN10      = 2.30258509299404568402;
#endif

#ifndef M_PI
constexpr long double M_PI        = 3.14159265358979323846;
#endif

#ifndef M_PI_2
constexpr long double M_PI_2      = 1.57079632679489661923;
#endif

#ifndef M_PI_4
constexpr long double M_PI_4      = 0.78539816339744830962;
#endif

#ifndef M_1_PI
constexpr long double M_1_PI      = 0.31830988618379067154;
#endif

#ifndef M_2_PI
constexpr long double M_2_PI      = 0.63661977236758134308;
#endif

#ifndef M_2_SQRTPI
constexpr long double M_2_SQRTPI  = 1.12837916709551257390;
#endif

#ifndef M_SQRT2
constexpr long double M_SQRT2     = 1.41421356237309504880;
#endif

#ifndef M_SQRT1_2
constexpr long double M_SQRT1_2   = 0.70710678118654752440;
#endif

#ifndef M_SQRT3
constexpr long double M_SQRT3     = 1.73205080756887729353;
#endif

#ifndef M_SQRT5
constexpr long double M_SQRT5     = 2.23606797749978969640;
#endif

#ifndef M_SQRTPI
constexpr long double M_SQRTPI    = 1.77245385090551602730;
#endif

#ifndef M_SQRT2PI
constexpr long double M_SQRT2PI   = 2.50662827463100050242;
#endif

#ifndef M_1_SQRTPI
constexpr long double M_1_SQRTPI  = 0.56418958354775628695;
#endif

#ifndef M_1_SQRT2PI
constexpr long double M_1_SQRT2PI = 0.39894228040143267794;
#endif

#ifndef M_1_E
constexpr long double M_1_E       = 0.36787944117144232160;
#endif

#ifndef M_LNPI
constexpr long double M_LNPI      = 1.14472988585494001741;
#endif

#ifndef M_EULER
constexpr long double M_EULER     = 0.57721566490153286061;
#endif

#ifndef M_APERY
constexpr long double M_APERY     = 1.20205690315959428539;
#endif

#ifndef M_CATALAN
constexpr long double M_CATALAN   = 0.91596559417721901505;
#endif

#ifndef M_PI_SQ
constexpr long double M_PI_SQ     = 9.86960440108935861883;
#endif


#endif // CONSTANTS_H
