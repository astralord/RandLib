#pragma once

#include <randlib/math/constants.hpp>

namespace randlib
{

/// Beta-related functions

namespace math
{

/**
 * @fn logBeta
 * Calculate logarithm of beta function
 * @param a positive parameter
 * @param b positive parameter
 * @return log(B(a, b)) = log(Γ(a)) + log(Γ(b)) - log(Γ(a + b))
 */
long double logBeta(long double a, long double b);

/**
 * @fn ibeta
 * Fast calculation of regularized beta function, using precalculated values
 * @param x
 * @param a
 * @param b
 * @param logBetaFun log(B(a, b))
 * @param logX log(x)
 * @param log1mX log(1-x)
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
double ibeta(double x, double a, double b, double logBetaFun, double logX, double log1mX);

/**
 * @fn ibeta
 * Calculate regularized beta function
 * @param x
 * @param a
 * @param b
 * @return I(x, a, b) = B(x, a, b) / B(a, b)
 */
double ibeta(double x, double a, double b);

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

inline long double logBeta(long double a, long double b)
{
    if (a <= 0 || b <= 0)
        throw std::invalid_argument("Parameters (a, b) of beta function should be both positive");
    /// Stirling corrections avoid both overflow of Γ and cancellation between large log Γ values
    long double p = std::min(a, b), q = std::max(a, b);
    if (p >= 10) {
        long double corr = lgammaCorrection(p) + lgammaCorrection(q) - lgammaCorrection(p + q);
        return -0.5l * std::log(q) + 0.5l * (randlib::consts::ln2 + randlib::consts::lnpi) + corr +
               (p - 0.5l) * std::log(p / (p + q)) + q * std::log1p(-p / (p + q));
    }
    if (q >= 10) {
        long double corr = lgammaCorrection(q) - lgammaCorrection(p + q);
        return std::lgamma(p) + corr + p - p * std::log(p + q) + (q - 0.5l) * std::log1p(-p / (p + q));
    }
    return std::log(std::tgamma(p) * (std::tgamma(q) / std::tgamma(p + q)));
}

inline double ibetaPowerSeries1(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    double addon = 1.0, sum = 0.0;
    /// avoid underflow
    double u = (x > 1e-5) ? -x / (1.0 - x) : -std::exp(logX - log1mX);
    int n = 1;
    do {
        addon *= (n - b) * u / (a + n);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(sum));
    double y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y += std::log1pl(sum);
    return std::exp(y) / a;
}

inline double ibetacontinuedFraction1(double x, double a, double b, int number, double logBetaFun, double logX,
                                      double log1mX)
{
    int n = number;
    /// avoid underflow
    double u = (x > 1e-5) ? x / (1.0 - x) : std::exp(logX - log1mX);
    double frac = 0.0;
    while (n > 1) {
        double cn = (b - n + 1) / (a + n - 1) * u;
        frac = cn / (1 + cn - frac);
        --n;
    }
    double c1 = b / a;
    double F = c1 / (1 + c1 - frac);
    double y = a * logX;
    y += (b - 1) * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / b;
    y *= F / (1.0 - F);
    return y;
}

inline double ibetaPowerSeries2(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    double addon = 1.0, sum = 0.0;
    int n = 0;
    do {
        addon *= (n + a + b) * x / (a + n + 1);
        sum += addon;
        ++n;
    } while (std::fabs(addon) > math::detail::minimumTolerance * std::fabs(sum));
    double y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y += std::log1pl(sum);
    return std::exp(y) / a;
}

inline double ibetaContinuedFraction2(double x, double a, double b, int number, double logBetaFun, double logX,
                                      double log1mX)
{
    int n = number;
    double frac = 0.0;
    while (n > 1) {
        double cn = (a + b + n - 2) / (a + n - 1) * x;
        frac = cn / (1 + cn - frac);
        --n;
    }
    double c1 = (a + b - 1) / a;
    double F = c1 / (1 + c1 - frac);
    double y = a * logX;
    y += b * log1mX;
    y -= logBetaFun;
    y = std::exp(y) / (a + b - 1);
    y *= F / (1.0 - F);
    return y;
}

inline double ibeta(double x, double a, double b, double logBetaFun, double logX, double log1mX)
{
    /// Check special and incorrect values of parameters
    if (a <= 0 || b <= 0 || x < 0.0 || x > 1.0)
        throw std::invalid_argument(
            "Parameters (a, b) of incomplete beta function should be both positive, but they're equal to " +
            std::to_string(a) + " and " + std::to_string(b) + " respectively");
    if (x < 0.0 || x > 1.0)
        throw std::invalid_argument(
            "Argument x of incomplete beta function should be in interval [0, 1], but it's equal to " +
            std::to_string(x));
    if (x == 0.0)
        return 0.0;
    if (x == 1.0)
        return 1.0;
    if (b == 1.0)
        return std::exp(a * logX);
    if (a == 1.0)
        return -std::expm1l(b * log1mX);

    /// If x is greater than mean of Beta distribution,
    /// convergence of complementary distribution function is faster
    if (x > a / (a + b))
        return 1.0 - ibeta(1.0 - x, b, a);

    /// We truncate after 0.45 instead of 0.5 here
    /// in order to avoid x being too close to 0.5
    if (x < 0.45) {
        if (b < 1)
            return ibetaPowerSeries1(x, a, b, logBetaFun, logX, log1mX);
        double number = std::floor(b);
        /// equivalent continued fraction #1
        double ecd = ibetacontinuedFraction1(x, a, b, number, logBetaFun, logX, log1mX);
        double apn = a + number, bmn = b - number;
        double residue = (b == number) ? 0.0 : ibetaPowerSeries1(x, apn, bmn, logBeta(apn, bmn), logX, log1mX);
        return ecd + residue;
    }

    if (b < 1)
        return ibetaPowerSeries2(x, a, b, logBetaFun, logX, log1mX);
    double number = std::floor(b);
    /// equivalent continued fraction #2
    double ecd = ibetaContinuedFraction2(x, a, b, number, logBetaFun, logX, log1mX);
    double apn = a + number;
    double residue = ibetaPowerSeries2(x, apn, b, logBeta(apn, b), logX, log1mX);
    return ecd + residue;
}

inline double ibeta(double x, double a, double b)
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
        double y = b * std::log1pl(-x);
        return -std::expm1l(y);
    }
    double logBetaFun = logBeta(a, b);
    double logX = std::log(x), log1mX = std::log1pl(-x);
    return ibeta(x, a, b, logBetaFun, logX, log1mX);
}

} // namespace math

} // namespace randlib
