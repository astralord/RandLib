#pragma once

#include <climits>
#include <cmath>
#include <numbers>
#include <string>

#include <randlib/core/fwd.hpp>

namespace randlib
{

namespace consts
{
inline constexpr long double e = std::numbers::e_v<long double>;
inline constexpr long double log2e = std::numbers::log2e_v<long double>;
inline constexpr long double log10e = std::numbers::log10e_v<long double>;
inline constexpr long double ln2 = std::numbers::ln2_v<long double>;
inline constexpr long double ln3 = 1.09861228866810969139524523692252570l;
inline constexpr long double ln10 = std::numbers::ln10_v<long double>;
inline constexpr long double lnpi = 1.14472988584940017414342735135305871l;
inline constexpr long double pi = std::numbers::pi_v<long double>;
inline constexpr long double half_pi = 1.57079632679489661923132169163975144l;
inline constexpr long double inv_pi = std::numbers::inv_pi_v<long double>;
inline constexpr long double two_inv_pi = 0.636619772367581343075535053490057448l;
inline constexpr long double two_inv_sqrtpi = 1.12837916709551257389615890312154517l;
inline constexpr long double sqrt2 = std::numbers::sqrt2_v<long double>;
inline constexpr long double inv_sqrt2 = 0.707106781186547524400844362104849039l;
inline constexpr long double sqrt3 = std::numbers::sqrt3_v<long double>;
inline constexpr long double sqrt5 = 2.23606797749978969640917366873127624l;
inline constexpr long double sqrtpi = 1.77245385090551602729816748334114518l;
inline constexpr long double sqrt2pi = 2.50662827463100050241576528481104525l;
inline constexpr long double inv_sqrtpi = std::numbers::inv_sqrtpi_v<long double>;
inline constexpr long double inv_sqrt2pi = 0.398942280401432677939946059934381868l;
inline constexpr long double inv_e = 0.367879441171442321595523770161460867l;
inline constexpr long double egamma = std::numbers::egamma_v<long double>;
inline constexpr long double apery = 1.20205690315959428539973816151144999l;
inline constexpr long double catalan = 0.915965594177219015054603514932384111l;
inline constexpr long double pi_sq = 9.86960440108935861883449099987615114l;
} // namespace consts

namespace math::detail
{
inline constexpr long double minimumTolerance = 1e-21L;
} // namespace math::detail

} // namespace randlib
