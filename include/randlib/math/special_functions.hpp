#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <randlib/math/constants.hpp>

namespace randlib
{

/// Special functions that the C++17 standard library provides only optionally.
///
/// libstdc++ and MSVC implement the C++17 mathematical special functions, libc++
/// does not. These wrappers forward to std:: when it is available and fall back to
/// RandLib's own implementations otherwise (or always, with
/// RANDLIB_PORTABLE_SPECIAL_FUNCTIONS defined).

namespace math
{

/**
 * @fn besselI
 * @param nu order, ν ≥ 0
 * @param x argument, x ≥ 0
 * @return modified Bessel function of the first kind I_ν(x)
 */
long double besselI(double nu, double x);

/**
 * @fn besselK
 * @param nu order, ν ≥ 0
 * @param x argument, x > 0
 * @return modified Bessel function of the second kind K_ν(x)
 */
long double besselK(double nu, double x);

/**
 * @fn riemannZeta
 * @param s argument, s > 1
 * @return Riemann zeta function ζ(s)
 */
long double riemannZeta(double s);

/**
 * @fn lgammaCorrection
 * @param x argument, x ≥ 10
 * @return log Γ(x) - ((x - 1/2) log x - x + log √(2π)), the remainder of Stirling's formula
 */
double lgammaCorrection(double x);

namespace portable
{
long double besselI(double nu, double x);
long double besselK(double nu, double x);
long double riemannZeta(double s);
}

}

} // namespace randlib

#include <randlib/math/beta.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/math/gamma.hpp>
#include <randlib/math/numeric.hpp>

namespace randlib
{

#if !defined(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS) && defined(__cpp_lib_math_special_functions)
#define RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
#endif

namespace math
{

constexpr double EPS = std::numeric_limits<double>::epsilon();
constexpr int SERIES_MAX_ITER = 100000;

/**
 * @fn reciprocalGammaTerms
 * Values needed by Temme's series for |μ| ≤ 1/2:
 * gam1 = (1/Γ(1-μ) - 1/Γ(1+μ)) / (2μ), gam2 = (1/Γ(1-μ) + 1/Γ(1+μ)) / 2,
 * gampl = 1/Γ(1+μ), gammi = 1/Γ(1-μ).
 */
inline void reciprocalGammaTerms(double mu, double &gam1, double &gam2, double &gampl, double &gammi)
{
    gampl = 1.0 / std::tgamma(1.0 + mu);
    gammi = 1.0 / std::tgamma(1.0 - mu);
    gam2 = 0.5 * (gammi + gampl);
    if (std::fabs(mu) < 0.01) {
        /// Odd part of the Taylor series 1/Γ(1+z) = Σ a_{k+1} z^k
        static constexpr double a2 = 0.5772156649015328606, a4 = -0.0420026350340952355,
                                a6 = -0.0421977345555443367, a8 = 0.0072189432466630995;
        double muSq = mu * mu;
        gam1 = -(a2 + muSq * (a4 + muSq * (a6 + muSq * a8)));
    }
    else {
        gam1 = (gammi - gampl) / (2 * mu);
    }
}

/**
 * @fn besselIK
 * Temme's method with Steed's continued fraction (Numerical Recipes, "bessik")
 * @param nu order ν ≥ 0
 * @param x argument x > 0
 * @param ri I_ν(x)
 * @param rk K_ν(x)
 */
inline void besselIK(double nu, double x, double &ri, double &rk)
{
    constexpr double FPMIN = std::numeric_limits<double>::min() / EPS;
    constexpr double XMIN = 2.0;
    const int nl = static_cast<int>(nu + 0.5);
    const double xmu = nu - nl, xmu2 = xmu * xmu;
    const double xi = 1.0 / x, xi2 = 2.0 * xi;

    /// CF1 for I'_ν / I_ν by the modified Lentz method
    double h = std::max(nu * xi, FPMIN);
    double b = xi2 * nu, d = 0.0, c = h;
    for (int i = 0; i < SERIES_MAX_ITER; ++i) {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;
        double del = c * d;
        h *= del;
        if (std::fabs(del - 1.0) < EPS)
            break;
    }

    /// Downward recurrence for I to order μ
    double ril = FPMIN, ripl = h * ril;
    const double ril1 = ril;
    double fact = nu * xi;
    for (int l = nl - 1; l >= 0; --l) {
        double ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    const double f = ripl / ril;

    double rkmu = 0.0, rk1 = 0.0;
    if (x < XMIN) {
        /// Temme's series for K_μ and K_{μ+1}
        double x2 = 0.5 * x, pimu = randlib::consts::pi * xmu;
        double factor = (std::fabs(pimu) < EPS) ? 1.0 : pimu / std::sin(pimu);
        double dd = -std::log(x2), e = xmu * dd;
        double factor2 = (std::fabs(e) < EPS) ? 1.0 : std::sinh(e) / e;
        double gam1, gam2, gampl, gammi;
        reciprocalGammaTerms(xmu, gam1, gam2, gampl, gammi);
        double ff = factor * (gam1 * std::cosh(e) + gam2 * factor2 * dd);
        double sum = ff;
        e = std::exp(e);
        double p = 0.5 * e / gampl, q = 0.5 / (e * gammi), cc = 1.0, x2Sq = x2 * x2;
        double sum1 = p;
        for (int i = 1; i <= SERIES_MAX_ITER; ++i) {
            ff = (i * ff + p + q) / (i * i - xmu2);
            cc *= x2Sq / i;
            p /= i - xmu;
            q /= i + xmu;
            double del = cc * ff;
            sum += del;
            sum1 += cc * (p - i * ff);
            if (std::fabs(del) < std::fabs(sum) * EPS)
                break;
        }
        rkmu = sum;
        rk1 = sum1 * xi2;
    }
    else {
        /// Steed's CF2 for K_μ and K_{μ+1}
        double bb = 2.0 * (1.0 + x), dd = 1.0 / bb, hh = dd, delh = dd;
        double q1 = 0.0, q2 = 1.0, a1 = 0.25 - xmu2, q = a1, cc = a1, a = -a1;
        double s = 1.0 + q * delh;
        for (int i = 1; i < SERIES_MAX_ITER; ++i) {
            a -= 2 * i;
            cc = -a * cc / (i + 1.0);
            double qnew = (q1 - bb * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += cc * qnew;
            bb += 2.0;
            dd = 1.0 / (bb + a * dd);
            delh = (bb * dd - 1.0) * delh;
            hh += delh;
            double dels = q * delh;
            s += dels;
            if (std::fabs(dels / s) < EPS)
                break;
        }
        hh *= a1;
        rkmu = std::sqrt(randlib::consts::pi / (2.0 * x)) * std::exp(-x) / s;
        rk1 = rkmu * (xmu + x + 0.5 - hh) * xi;
    }

    /// Wronskian gives I_μ, then scale the recurrence and recur K upwards
    double rkmup = xmu * xi * rkmu - rk1;
    double rimu = xi / (f * rkmu - rkmup);
    ri = (rimu * ril1) / ril;
    for (int i = 1; i <= nl; ++i) {
        double rktemp = (xmu + i) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    rk = rkmu;
}

inline void checkBesselArguments(double nu, double x)
{
    if (nu < 0.0)
        throw std::domain_error("Bessel function: order should be non-negative, but it's equal to " + std::to_string(nu));
    if (x < 0.0)
        throw std::domain_error("Bessel function: argument should be non-negative, but it's equal to " + std::to_string(x));
}

namespace portable
{

inline long double besselI(double nu, double x)
{
    checkBesselArguments(nu, x);
    if (x == 0.0)
        return (nu == 0.0) ? 1.0 : 0.0;
    if (x > 700.0) {
        /// I_ν(x) overflows; keep the behaviour of std::cyl_bessel_i
        return std::numeric_limits<long double>::infinity();
    }
    if (x < 5.0 || x * x < 4.0 * (nu + 1.0)) {
        /// Power series I_ν(x) = (x/2)^ν Σ (x²/4)^k / (k! Γ(ν + k + 1)); every term is positive
        double y = 0.25 * x * x, term = 1.0, sum = 1.0;
        for (int k = 1; k < SERIES_MAX_ITER; ++k) {
            term *= y / (k * (nu + k));
            sum += term;
            if (term < EPS * sum)
                break;
        }
        return std::exp(nu * std::log(0.5 * x) - std::lgamma(nu + 1.0)) * sum;
    }
    double ri = 0.0, rk = 0.0;
    besselIK(nu, x, ri, rk);
    return ri;
}

inline long double besselK(double nu, double x)
{
    checkBesselArguments(nu, x);
    if (x == 0.0)
        return std::numeric_limits<long double>::infinity();
    if (x > 745.0)
        return 0.0;
    double ri = 0.0, rk = 0.0;
    besselIK(nu, x, ri, rk);
    return rk;
}

inline long double riemannZeta(double s)
{
    if (!(s > 1.0))
        throw std::domain_error("Riemann zeta function: argument should be greater than 1, but it's equal to " + std::to_string(s));
    if (s > 60.0)
        return 1.0L + std::pow(2.0L, -s) + std::pow(3.0L, -s);

    /// Euler-Maclaurin summation with N terms and Bernoulli corrections up to B_20
    static constexpr int N = 20;
    static constexpr long double B2jOverFact[] = {
        1.0L / 12, -1.0L / 720, 1.0L / 30240, -1.0L / 1209600, 1.0L / 47900160,
        -691.0L / 1307674368000.0L, 1.0L / 74724249600.0L, -3617.0L / 10670622842880000.0L,
        43867.0L / 5109094217170944000.0L, -174611.0L / 802857662698291200000.0L
    };
    long double sum = 0.0L;
    for (int k = 1; k < N; ++k)
        sum += std::pow(static_cast<long double>(k), -static_cast<long double>(s));
    long double nPow = std::pow(static_cast<long double>(N), -static_cast<long double>(s));
    sum += N * nPow / (s - 1.0L) + 0.5L * nPow;
    long double rising = s, term = nPow / N;
    for (int j = 0; j < 10; ++j) {
        sum += B2jOverFact[j] * rising * term;
        rising *= (s + 2 * j + 1) * (s + 2 * j + 2);
        term /= static_cast<long double>(N) * N;
    }
    return sum;
}

}

inline long double besselI(double nu, double x)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::cyl_bessel_il(nu, x);
#else
    return portable::besselI(nu, x);
#endif
}

inline long double besselK(double nu, double x)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::cyl_bessel_kl(nu, x);
#else
    return portable::besselK(nu, x);
#endif
}

inline long double riemannZeta(double s)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::riemann_zetal(s);
#else
    return portable::riemannZeta(s);
#endif
}

inline double lgammaCorrection(double x)
{
    if (x < 10.0)
        throw std::domain_error("Stirling correction: argument should be at least 10, but it's equal to " + std::to_string(x));
    static constexpr double c[] = {1.0 / 12, -1.0 / 360, 1.0 / 1260, -1.0 / 1680, 1.0 / 1188, -691.0 / 360360, 1.0 / 156};
    double xInv = 1.0 / x, xInvSq = xInv * xInv;
    double sum = c[6];
    for (int i = 5; i >= 0; --i)
        sum = sum * xInvSq + c[i];
    return sum * xInv;
}

}

} // namespace randlib
