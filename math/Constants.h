#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cmath>
#include <climits>
#include <string>
#include <tuple>
#include <utility>

namespace randlib::consts
{
constexpr long double e              = 2.71828182845904523536028747135266250l;
constexpr long double log2e          = 1.44269504088896340735992468100189214l;
constexpr long double log10e         = 0.434294481903251827651128918916605082l;
constexpr long double ln2            = 0.693147180559945309417232121458176568l;
constexpr long double ln3            = 1.09861228866810969139524523692252570l;
constexpr long double ln10           = 2.30258509299404568401799145468436421l;
constexpr long double lnpi           = 1.14472988584940017414342735135305871l;
constexpr long double pi             = 3.14159265358979323846264338327950288l;
constexpr long double half_pi        = 1.57079632679489661923132169163975144l;
constexpr long double inv_pi         = 0.318309886183790671537767526745028724l;
constexpr long double two_inv_pi     = 0.636619772367581343075535053490057448l;
constexpr long double two_inv_sqrtpi = 1.12837916709551257389615890312154517l;
constexpr long double sqrt2          = 1.41421356237309504880168872420969808l;
constexpr long double inv_sqrt2      = 0.707106781186547524400844362104849039l;
constexpr long double sqrt3          = 1.73205080756887729352744634150587237l;
constexpr long double sqrt5          = 2.23606797749978969640917366873127624l;
constexpr long double sqrtpi         = 1.77245385090551602729816748334114518l;
constexpr long double sqrt2pi        = 2.50662827463100050241576528481104525l;
constexpr long double inv_sqrtpi     = 0.564189583547756286948079451560772586l;
constexpr long double inv_sqrt2pi    = 0.398942280401432677939946059934381868l;
constexpr long double inv_e          = 0.367879441171442321595523770161460867l;
constexpr long double egamma         = 0.577215664901532860606512090082402431l;
constexpr long double apery          = 1.20205690315959428539973816151144999l;
constexpr long double catalan        = 0.915965594177219015054603514932384111l;
constexpr long double pi_sq          = 9.86960440108935861883449099987615114l;
}

constexpr long double MIN_POSITIVE = 1e-21l;

template <typename T>
using Pair = std::pair<T, T>;

template < typename T >
Pair<T> operator +(const Pair<T> & l, const Pair<T> & r) {
    return {l.first + r.first, l.second + r.second};
}

template < typename T >
Pair<T> operator -(const Pair<T> & l, const Pair<T> & r) {
    return {l.first - r.first, l.second - r.second};
}

template < typename T >
T operator *(const Pair<T> & l, const Pair<T> & r) {
    return l.first * r.first + l.second * r.second;
}

template <typename T>
using Triplet = std::tuple<T, T, T>;

using IntPair = Pair<int>;
using DoublePair = Pair<double>;
using LongDoublePair = Pair<long double>;

using IntTriplet = Triplet<int>;
using DoubleTriplet = Triplet<double>;
using LongDoubleTriplet = Triplet<long double>;

using String = std::string;

/// Legacy names, kept for source compatibility in the 1.x series.
/// They are only defined when the C library does not already provide them.
/// New code should use randlib::consts.
#ifndef M_E
constexpr long double M_E = randlib::consts::e;
#endif
#ifndef M_LOG2E
constexpr long double M_LOG2E = randlib::consts::log2e;
#endif
#ifndef M_LOG10E
constexpr long double M_LOG10E = randlib::consts::log10e;
#endif
#ifndef M_LN2
constexpr long double M_LN2 = randlib::consts::ln2;
#endif
#ifndef M_LN3
constexpr long double M_LN3 = randlib::consts::ln3;
#endif
#ifndef M_LN10
constexpr long double M_LN10 = randlib::consts::ln10;
#endif
#ifndef M_LNPI
constexpr long double M_LNPI = randlib::consts::lnpi;
#endif
#ifndef M_PI
constexpr long double M_PI = randlib::consts::pi;
#endif
#ifndef M_PI_2
constexpr long double M_PI_2 = randlib::consts::half_pi;
#endif
#ifndef M_1_PI
constexpr long double M_1_PI = randlib::consts::inv_pi;
#endif
#ifndef M_2_PI
constexpr long double M_2_PI = randlib::consts::two_inv_pi;
#endif
#ifndef M_2_SQRTPI
constexpr long double M_2_SQRTPI = randlib::consts::two_inv_sqrtpi;
#endif
#ifndef M_SQRT2
constexpr long double M_SQRT2 = randlib::consts::sqrt2;
#endif
#ifndef M_SQRT1_2
constexpr long double M_SQRT1_2 = randlib::consts::inv_sqrt2;
#endif
#ifndef M_SQRT3
constexpr long double M_SQRT3 = randlib::consts::sqrt3;
#endif
#ifndef M_SQRT5
constexpr long double M_SQRT5 = randlib::consts::sqrt5;
#endif
#ifndef M_SQRTPI
constexpr long double M_SQRTPI = randlib::consts::sqrtpi;
#endif
#ifndef M_SQRT2PI
constexpr long double M_SQRT2PI = randlib::consts::sqrt2pi;
#endif
#ifndef M_1_SQRTPI
constexpr long double M_1_SQRTPI = randlib::consts::inv_sqrtpi;
#endif
#ifndef M_1_SQRT2PI
constexpr long double M_1_SQRT2PI = randlib::consts::inv_sqrt2pi;
#endif
#ifndef M_1_E
constexpr long double M_1_E = randlib::consts::inv_e;
#endif
#ifndef M_EULER
constexpr long double M_EULER = randlib::consts::egamma;
#endif
#ifndef M_APERY
constexpr long double M_APERY = randlib::consts::apery;
#endif
#ifndef M_CATALAN
constexpr long double M_CATALAN = randlib::consts::catalan;
#endif
#ifndef M_PI_SQ
constexpr long double M_PI_SQ = randlib::consts::pi_sq;
#endif

#endif // CONSTANTS_H
