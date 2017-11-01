#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <climits>
#include <string>

constexpr long double MIN_POSITIVE = 1e-21l;

typedef std::pair <double, double> DoublePair;
typedef std::tuple <double, double, double> DoubleTriplet;
typedef std::pair <int, int> IntPair;
typedef std::string String;

#ifndef INFINITY
#include <limits>
long double INFINITY = std::numeric_limits<long double>::infinity()l;
#endif

#ifndef NAN
#include <limits>
long double NAN = std::numeric_limits<long double>::quiet_NaN()l;
#endif

#ifndef M_E
constexpr long double M_E         = 2.71828182845904523536l;
#endif

#ifndef M_LOG2E
constexpr long double M_LOG2E     = 1.44269504088896340760l;
#endif

#ifndef M_LOG10E
constexpr long double M_LOG10E    = 0.43429448190325182765l;
#endif

#ifndef M_LN2
constexpr long double M_LN2       = 0.69314718055994530942l;
#endif

#ifndef M_LN3
constexpr long double M_LN3       = 1.09861228866810969139l;
#endif

#ifndef M_LN10
constexpr long double M_LN10      = 2.30258509299404568402l;
#endif

#ifndef M_PI
constexpr long double M_PI        = 3.14159265358979323846l;
#endif

#ifndef M_PI_2
constexpr long double M_PI_2      = 1.57079632679489661923l;
#endif

#ifndef M_1_PI
constexpr long double M_1_PI      = 0.31830988618379067154l;
#endif

#ifndef M_2_PI
constexpr long double M_2_PI      = 0.63661977236758134308l;
#endif

#ifndef M_2_SQRTPI
constexpr long double M_2_SQRTPI  = 1.12837916709551257390l;
#endif

#ifndef M_SQRT2
constexpr long double M_SQRT2     = 1.41421356237309504880l;
#endif

#ifndef M_SQRT1_2
constexpr long double M_SQRT1_2   = 0.70710678118654752440l;
#endif

#ifndef M_SQRT3
constexpr long double M_SQRT3     = 1.73205080756887729353l;
#endif

#ifndef M_SQRT5
constexpr long double M_SQRT5     = 2.23606797749978969640l;
#endif

#ifndef M_SQRTPI
constexpr long double M_SQRTPI    = 1.77245385090551602730l;
#endif

#ifndef M_SQRT2PI
constexpr long double M_SQRT2PI   = 2.50662827463100050242l;
#endif

#ifndef M_1_SQRTPI
constexpr long double M_1_SQRTPI  = 0.56418958354775628695l;
#endif

#ifndef M_1_SQRT2PI
constexpr long double M_1_SQRT2PI = 0.39894228040143267794l;
#endif

#ifndef M_1_E
constexpr long double M_1_E       = 0.36787944117144232160l;
#endif

#ifndef M_LNPI
constexpr long double M_LNPI      = 1.14472988585494001741l;
#endif

#ifndef M_EULER
constexpr long double M_EULER     = 0.57721566490153286061l;
#endif

#ifndef M_APERY
constexpr long double M_APERY     = 1.20205690315959428539l;
#endif

#ifndef M_CATALAN
constexpr long double M_CATALAN   = 0.91596559417721901505l;
#endif

#ifndef M_PI_SQ
constexpr long double M_PI_SQ     = 9.86960440108935861883l;
#endif

#endif // CONSTANTS_H
