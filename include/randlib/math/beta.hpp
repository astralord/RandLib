#pragma once

#include <concepts>
#include <type_traits>

#include <randlib/math/constants.hpp>

namespace randlib
{

/// Beta-related functions

namespace math
{

/**
 * Calculate logarithm of beta function
 * @param a positive parameter
 * @param b positive parameter
 * @return log(B(a, b)) = log(Γ(a)) + log(Γ(b)) - log(Γ(a + b))
 */
template <std::floating_point Real>
Real logBeta(Real a, std::type_identity_t<Real> b);

template <typename Left, typename Right>
    requires(std::is_arithmetic_v<Left> && std::is_arithmetic_v<Right> &&
             (std::integral<Left> || std::integral<Right> || !std::same_as<Left, Right>))
auto logBeta(Left a, Right b) -> std::common_type_t<double, Left, Right>;

/**
 * Fast calculation of regularized beta function, using precalculated values
 * @param x
 * @param a
 * @param b
 * @param logBetaFun log(B(a, b))
 * @param logX log(x)
 * @param log1mX log(1-x)
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
template <std::floating_point Real>
Real ibeta(Real x, std::type_identity_t<Real> a, std::type_identity_t<Real> b, std::type_identity_t<Real> logBetaFun,
           std::type_identity_t<Real> logX, std::type_identity_t<Real> log1mX);

/**
 * Calculate regularized beta function
 * @param x
 * @param a
 * @param b
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
template <std::floating_point Real>
Real ibeta(Real x, std::type_identity_t<Real> a, std::type_identity_t<Real> b);

} // namespace math

} // namespace randlib

#include <randlib/math/functions.hpp>
#include <randlib/math/gamma.hpp>
#include <randlib/math/quadrature.hpp>
#include <randlib/math/roots.hpp>
#include <randlib/math/optimize.hpp>
#include <randlib/math/special_functions.hpp>

namespace randlib
{

namespace math
{

template <std::floating_point Real>
inline Real logBeta(Real a, std::type_identity_t<Real> b)
{
    if (a <= 0 || b <= 0)
        throw std::invalid_argument("Parameters (a, b) of beta function should be both positive");
    /// Stirling corrections avoid both overflow of Γ and cancellation between large log Γ values
    using Work = std::conditional_t<(sizeof(Real) < sizeof(long double)), long double, Real>;
    Work p = std::min(static_cast<Work>(a), static_cast<Work>(b));
    Work q = std::max(static_cast<Work>(a), static_cast<Work>(b));
    if (p >= 10) {
        Work corr = lgammaCorrection(p) + lgammaCorrection(q) - lgammaCorrection(p + q);
        Work result = -Work{0.5} * std::log(q) +
                      Work{0.5} * static_cast<Work>(randlib::consts::ln2 + randlib::consts::lnpi) + corr +
                      (p - Work{0.5}) * std::log(p / (p + q)) + q * std::log1p(-p / (p + q));
        return static_cast<Real>(result);
    }
    if (q >= 10) {
        Work corr = lgammaCorrection(q) - lgammaCorrection(p + q);
        Work result = std::lgamma(p) + corr + p - p * std::log(p + q) + (q - Work{0.5}) * std::log1p(-p / (p + q));
        return static_cast<Real>(result);
    }
    return static_cast<Real>(std::log(std::tgamma(p) * (std::tgamma(q) / std::tgamma(p + q))));
}

template <typename Left, typename Right>
    requires(std::is_arithmetic_v<Left> && std::is_arithmetic_v<Right> &&
             (std::integral<Left> || std::integral<Right> || !std::same_as<Left, Right>))
inline auto logBeta(Left a, Right b) -> std::common_type_t<double, Left, Right>
{
    using Real = std::common_type_t<double, Left, Right>;
    return logBeta(static_cast<Real>(a), static_cast<Real>(b));
}

template <std::floating_point Real>
inline Real ibetaPowerSeries1(Real x, Real a, Real b, Real logBetaFun, Real logX, Real log1mX)
{
    Real addon = 1.0, sum = 0.0;
    /// avoid underflow
    Real u = (x > 1e-5) ? -x / (1.0 - x) : -std::exp(logX - log1mX);
    int n = 1;
    do {
        addon *= (n - b) * u / (a + n);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(sum));
    Real y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y += std::log1p(sum);
    return std::exp(y) / a;
}

template <std::floating_point Real>
inline Real ibetacontinuedFraction1(Real x, Real a, Real b, int number, Real logBetaFun, Real logX, Real log1mX)
{
    int n = number;
    /// avoid underflow
    Real u = (x > 1e-5) ? x / (1.0 - x) : std::exp(logX - log1mX);
    Real frac = 0.0;
    while (n > 1) {
        Real cn = (b - n + 1) / (a + n - 1) * u;
        frac = cn / (1 + cn - frac);
        --n;
    }
    Real c1 = b / a;
    Real F = c1 / (1 + c1 - frac);
    Real y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / b;
    y *= F / (1.0 - F);
    return y;
}

template <std::floating_point Real>
inline Real ibetaPowerSeries2(Real x, Real a, Real b, Real logBetaFun, Real logX, Real log1mX)
{
    Real addon = 1.0, sum = 0.0;
    int n = 0;
    do {
        addon *= (n + a + b) * x / (a + n + 1);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(sum));
    Real y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y += std::log1p(sum);
    return std::exp(y) / a;
}

template <std::floating_point Real>
inline Real ibetaContinuedFraction2(Real x, Real a, Real b, int number, Real logBetaFun, Real logX, Real log1mX)
{
    int n = number;
    Real frac = 0.0;
    while (n > 1) {
        Real cn = (a + b + n - 2) / (a + n - 1) * x;
        frac = cn / (1 + cn - frac);
        --n;
    }
    Real c1 = (a + b - 1) / a;
    Real F = c1 / (1 + c1 - frac);
    Real y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / (a + b - 1);
    y *= F / (1.0 - F);
    return y;
}

/// Continued fraction of I(x, a, b) (DiDonato and Morris, also Numerical Recipes betacf) by the modified Lentz
/// method; it converges for x below (a + 1) / (a + b + 2) in far fewer terms than b for large shapes
template <std::floating_point Real>
inline Real ibetaContinuedFractionLentz(Real x, Real a, Real b, Real logBetaFun, Real logX, Real log1mX)
{
    constexpr Real epsilon = std::numeric_limits<Real>::epsilon();
    constexpr Real tiny = std::numeric_limits<Real>::min() / epsilon;
    auto guard = [](Real value) { return (std::fabs(value) < tiny) ? tiny : value; };
    Real c = 1, d = 1 / guard(1 - (a + b) * x / (a + 1));
    Real fraction = d;
    for (int m = 1; m < 1000000; ++m) {
        const Real even = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
        d = 1 / guard(1 + even * d);
        c = guard(1 + even / c);
        fraction *= d * c;
        const Real odd = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
        d = 1 / guard(1 + odd * d);
        c = guard(1 + odd / c);
        const Real step = d * c;
        fraction *= step;
        if (std::fabs(step - 1) <= epsilon)
            break;
    }
    return std::exp(a * logX + b * log1mX - logBetaFun) * fraction / a;
}

namespace detail
{

/// Beyond this size of a shape the finite continued fractions of ibetaLowerTail accumulate rounding errors, and
/// its power series converge at the rate x close to 1
inline constexpr double IBETA_LARGEST_FINITE_SHAPE = 4096;

} // namespace detail

/// I(x, a, b) for x not above the mean a / (a + b) and shapes up to IBETA_LARGEST_FINITE_SHAPE
template <std::floating_point Real>
inline Real ibetaLowerTail(Real x, Real a, Real b, Real logBetaFun, Real logX, Real log1mX)
{
    if (b == 1.0)
        return std::exp(a * logX);
    if (a == 1.0)
        return -std::expm1(b * log1mX);

    /// We truncate after 0.45 instead of 0.5 here
    /// in order to avoid x being too close to 0.5
    if (x < 0.45) {
        if (b < 1)
            return ibetaPowerSeries1(x, a, b, logBetaFun, logX, log1mX);
        int number = static_cast<int>(std::floor(b));
        /// equivalent continued fraction #1
        Real ecd = ibetacontinuedFraction1(x, a, b, number, logBetaFun, logX, log1mX);
        Real apn = a + static_cast<Real>(number), bmn = b - static_cast<Real>(number);
        Real residue = (b == number) ? 0.0 : ibetaPowerSeries1(x, apn, bmn, logBeta(apn, bmn), logX, log1mX);
        return ecd + residue;
    }

    if (b < 1)
        return ibetaPowerSeries2(x, a, b, logBetaFun, logX, log1mX);
    int number = static_cast<int>(std::floor(b));
    /// equivalent continued fraction #2
    Real ecd = ibetaContinuedFraction2(x, a, b, number, logBetaFun, logX, log1mX);
    Real apn = a + static_cast<Real>(number);
    Real residue = ibetaPowerSeries2(x, apn, b, logBeta(apn, b), logX, log1mX);
    return ecd + residue;
}

template <std::floating_point Real>
inline Real ibeta(Real x, std::type_identity_t<Real> a, std::type_identity_t<Real> b,
                  std::type_identity_t<Real> logBetaFun, std::type_identity_t<Real> logX,
                  std::type_identity_t<Real> log1mX)
{
    /// Check special and incorrect values of parameters
    if (!(a > 0 && b > 0))
        throw std::invalid_argument(
            "Parameters (a, b) of incomplete beta function should be both positive, but they're equal to " +
            std::to_string(a) + " and " + std::to_string(b) + " respectively");
    if (!(x >= 0.0 && x <= 1.0))
        throw std::invalid_argument(
            "Argument x of incomplete beta function should be in interval [0, 1], but it's equal to " +
            std::to_string(x));
    /// x may be a rounded complement such as 1 - t for tiny t, so 1 - x is recovered from log(1 - x)
    if (x == 0.0 || logX == -INFINITY)
        return 0.0;
    if (log1mX == -INFINITY)
        return 1.0;
    const Real y = (x < 0.5) ? Real{1} - x : std::exp(log1mX);

    /// For large shapes the continued fraction converges below c = (a + 1) / (a + b + 2). For a tiny x the
    /// complementary fraction at 1 - x would lose its digits, and the direct one still converges up to about
    /// four standard deviations above c, so the side whose argument is tiny is preferred there
    if (std::max(a, b) > detail::IBETA_LARGEST_FINITE_SHAPE) {
        constexpr Real TINY_ARGUMENT = 1e-4, DEVIATIONS = 4;
        const Real c = (a + 1) / (a + b + 2);
        const bool directExtends = x < TINY_ARGUMENT && x < c * (1 + DEVIATIONS / std::sqrt(a + 1));
        const bool complementExtends = y < TINY_ARGUMENT && y < (1 - c) * (1 + DEVIATIONS / std::sqrt(b + 1));
        const bool direct = (x < c) ? !complementExtends : directExtends;
        if (direct)
            return ibetaContinuedFractionLentz<Real>(x, a, b, logBetaFun, logX, log1mX);
        return Real{1} - ibetaContinuedFractionLentz<Real>(y, b, a, logBetaFun, log1mX, logX);
    }

    /// If x is greater than mean of Beta distribution,
    /// convergence of complementary distribution function is faster.
    /// Close to 1 the series below converge only at the rate x, so the argument is never moved there: a
    /// tiny argument beyond a tinier mean is evaluated directly, and one within this distance of 1 under
    /// the mean (for a tiny second shape) through the complement, losing at most ε / I(x, a, b)
    constexpr Real SLOW_SERIES_DISTANCE = 1e-3;
    const bool complement = (x > a / (a + b)) ? x >= SLOW_SERIES_DISTANCE : y < SLOW_SERIES_DISTANCE;
    if (complement)
        return Real{1} - ibetaLowerTail<Real>(y, b, a, logBetaFun, log1mX, logX);
    return ibetaLowerTail<Real>(x, a, b, logBetaFun, logX, log1mX);
}

template <std::floating_point Real>
inline Real ibeta(Real x, std::type_identity_t<Real> a, std::type_identity_t<Real> b)
{
    /// Check special and incorrect values of parameters
    if (a <= 0 || b <= 0 || x < 0.0 || x > 1.0)
        throw std::invalid_argument("Parameters (a, b) of incomplete beta function should be both positive");
    if (x < 0.0 || x > 1.0)
        throw std::invalid_argument("Argument x of incomplete beta function should be in interval [0, 1]");
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    if (b == 1.0)
        return std::pow(x, a);
    if (a == 1.0) {
        Real y = b * std::log1p(-x);
        return -std::expm1(y);
    }
    Real logBetaFun = logBeta(a, b);
    Real logX = std::log(x), log1mX = std::log1p(-x);
    return ibeta(x, a, b, logBetaFun, logX, log1mX);
}

} // namespace math

} // namespace randlib
