#pragma once

#include <cmath>
#include <concepts>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

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
template <std::floating_point Real>
Real besselI(std::type_identity_t<Real> nu, Real x);

/**
 * @fn besselK
 * @param nu order, ν ≥ 0
 * @param x argument, x > 0
 * @return modified Bessel function of the second kind K_ν(x)
 */
template <std::floating_point Real>
Real besselK(std::type_identity_t<Real> nu, Real x);

/**
 * @fn riemannZeta
 * @param s argument, s > 1
 * @return Riemann zeta function ζ(s)
 */
template <std::floating_point Real>
Real riemannZeta(Real s);

/**
 * @fn lgammaCorrection
 * @param x argument, x ≥ 10
 * @return log Γ(x) - ((x - 1/2) log x - x + log √(2π)), the remainder of Stirling's formula
 */
template <std::floating_point Real>
Real lgammaCorrection(Real x);

namespace portable
{
template <std::floating_point Real>
Real besselI(std::type_identity_t<Real> nu, Real x);
template <std::floating_point Real>
Real besselK(std::type_identity_t<Real> nu, Real x);
template <std::floating_point Real>
Real riemannZeta(Real s);
} // namespace portable

} // namespace math

} // namespace randlib

#include <randlib/math/beta.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/math/gamma.hpp>
#include <randlib/math/quadrature.hpp>
#include <randlib/math/roots.hpp>
#include <randlib/math/optimize.hpp>

namespace randlib
{

#if !defined(RANDLIB_PORTABLE_SPECIAL_FUNCTIONS) && defined(__cpp_lib_math_special_functions)
#define RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
#endif

namespace math
{

namespace detail
{

template <std::floating_point Real>
inline constexpr Real machineEpsilon = std::numeric_limits<Real>::epsilon();
constexpr int maxSeriesIterations = 100000;

/**
 * @fn reciprocalGammaTerms
 * Values needed by Temme's series for |μ| ≤ 1/2:
 * gam1 = (1/Γ(1-μ) - 1/Γ(1+μ)) / (2μ), gam2 = (1/Γ(1-μ) + 1/Γ(1+μ)) / 2,
 * gampl = 1/Γ(1+μ), gammi = 1/Γ(1-μ).
 */
template <std::floating_point Real>
inline void reciprocalGammaTerms(Real mu, Real &gam1, Real &gam2, Real &gampl, Real &gammi)
{
    gampl = 1.0 / std::tgamma(1.0 + mu);
    gammi = 1.0 / std::tgamma(1.0 - mu);
    gam2 = 0.5 * (gammi + gampl);
    if (std::fabs(mu) < 0.01) {
        /// Odd part of the Taylor series 1/Γ(1+z) = Σ a_{k+1} z^k
        static constexpr Real a2 = 0.5772156649015328606L, a4 = -0.0420026350340952355L, a6 = -0.0421977345555443367L,
                              a8 = 0.0072189432466630995L;
        Real muSq = mu * mu;
        gam1 = -(a2 + muSq * (a4 + muSq * (a6 + muSq * a8)));
    } else {
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
template <std::floating_point Real>
inline void besselIK(Real nu, Real x, Real &ri, Real &rk)
{
    constexpr Real FPMIN = std::numeric_limits<Real>::min() / machineEpsilon<Real>;
    constexpr Real XMIN = 2;
    const int nl = static_cast<int>(nu + 0.5);
    const Real xmu = nu - nl, xmu2 = xmu * xmu;
    const Real xi = 1 / x, xi2 = 2 * xi;

    /// CF1 for I'_ν / I_ν by the modified Lentz method
    Real h = std::max(nu * xi, FPMIN);
    Real b = xi2 * nu, d = 0, c = h;
    for (int i = 0; i < maxSeriesIterations; ++i) {
        b += xi2;
        d = 1.0 / (b + d);
        c = b + 1.0 / c;
        Real del = c * d;
        h *= del;
        if (std::fabs(del - 1) < machineEpsilon<Real>)
            break;
    }

    /// Downward recurrence for I to order μ
    Real ril = FPMIN, ripl = h * ril;
    const Real ril1 = ril;
    Real fact = nu * xi;
    for (int l = nl - 1; l >= 0; --l) {
        Real ritemp = fact * ril + ripl;
        fact -= xi;
        ripl = fact * ritemp + ril;
        ril = ritemp;
    }
    const Real f = ripl / ril;

    Real rkmu = 0, rk1 = 0;
    if (x < XMIN) {
        /// Temme's series for K_μ and K_{μ+1}
        Real x2 = Real{0.5} * x, pimu = static_cast<Real>(randlib::consts::pi) * xmu;
        Real factor = (std::fabs(pimu) < machineEpsilon<Real>) ? Real{1} : pimu / std::sin(pimu);
        Real dd = -std::log(x2), e = xmu * dd;
        Real factor2 = (std::fabs(e) < machineEpsilon<Real>) ? Real{1} : std::sinh(e) / e;
        Real gam1, gam2, gampl, gammi;
        reciprocalGammaTerms(xmu, gam1, gam2, gampl, gammi);
        Real ff = factor * (gam1 * std::cosh(e) + gam2 * factor2 * dd);
        Real sum = ff;
        e = std::exp(e);
        Real p = Real{0.5} * e / gampl, q = Real{0.5} / (e * gammi), cc = 1, x2Sq = x2 * x2;
        Real sum1 = p;
        for (int i = 1; i <= maxSeriesIterations; ++i) {
            ff = (i * ff + p + q) / (i * i - xmu2);
            cc *= x2Sq / i;
            p /= i - xmu;
            q /= i + xmu;
            Real del = cc * ff;
            sum += del;
            sum1 += cc * (p - i * ff);
            if (std::fabs(del) < std::fabs(sum) * machineEpsilon<Real>)
                break;
        }
        rkmu = sum;
        rk1 = sum1 * xi2;
    } else {
        /// Steed's CF2 for K_μ and K_{μ+1}
        Real bb = 2 * (1 + x), dd = 1 / bb, hh = dd, delh = dd;
        Real q1 = 0, q2 = 1, a1 = Real{0.25} - xmu2, q = a1, cc = a1, a = -a1;
        Real s = 1 + q * delh;
        for (int i = 1; i < maxSeriesIterations; ++i) {
            a -= 2 * i;
            cc = -a * cc / (i + 1.0);
            Real qnew = (q1 - bb * q2) / a;
            q1 = q2;
            q2 = qnew;
            q += cc * qnew;
            bb += 2.0;
            dd = 1.0 / (bb + a * dd);
            delh = (bb * dd - 1.0) * delh;
            hh += delh;
            Real dels = q * delh;
            s += dels;
            if (std::fabs(dels / s) < machineEpsilon<Real>)
                break;
        }
        hh *= a1;
        rkmu = std::sqrt(static_cast<Real>(randlib::consts::pi) / (2 * x)) * std::exp(-x) / s;
        rk1 = rkmu * (xmu + x + 0.5 - hh) * xi;
    }

    /// Wronskian gives I_μ, then scale the recurrence and recur K upwards
    Real rkmup = xmu * xi * rkmu - rk1;
    Real rimu = xi / (f * rkmu - rkmup);
    ri = (rimu * ril1) / ril;
    for (int i = 1; i <= nl; ++i) {
        Real rktemp = (xmu + i) * xi2 * rk1 + rkmu;
        rkmu = rk1;
        rk1 = rktemp;
    }
    rk = rkmu;
}

template <std::floating_point Real>
inline void checkBesselArguments(Real nu, Real x)
{
    if (nu < 0.0)
        throw std::domain_error("Bessel function: order should be non-negative, but it's equal to " +
                                std::to_string(nu));
    if (x < 0.0)
        throw std::domain_error("Bessel function: argument should be non-negative, but it's equal to " +
                                std::to_string(x));
}

} // namespace detail

namespace portable
{

template <std::floating_point Real>
inline Real besselI(std::type_identity_t<Real> nu, Real x)
{
    detail::checkBesselArguments(nu, x);
    if (x == 0.0)
        return (nu == 0.0) ? 1.0 : 0.0;
    const Real overflowLimit = std::numeric_limits<Real>::max_exponent == std::numeric_limits<double>::max_exponent
                                   ? Real{700}
                                   : std::log(std::numeric_limits<Real>::max()) - Real{10};
    if (x > overflowLimit) {
        /// I_ν(x) overflows; keep the behaviour of std::cyl_bessel_i
        return std::numeric_limits<Real>::infinity();
    }
    if (x < 5.0 || x * x < 4.0 * (nu + 1.0)) {
        /// Power series I_ν(x) = (x/2)^ν Σ (x²/4)^k / (k! Γ(ν + k + 1)); every term is positive
        Real y = Real{0.25} * x * x, term = 1, sum = 1;
        for (int k = 1; k < detail::maxSeriesIterations; ++k) {
            term *= y / (k * (nu + k));
            sum += term;
            if (term < detail::machineEpsilon<Real> * sum)
                break;
        }
        return std::exp(nu * std::log(0.5 * x) - std::lgamma(nu + 1.0)) * sum;
    }
    Real ri = 0, rk = 0;
    detail::besselIK(nu, x, ri, rk);
    return ri;
}

template <std::floating_point Real>
inline Real besselK(std::type_identity_t<Real> nu, Real x)
{
    detail::checkBesselArguments(nu, x);
    if (x == 0.0)
        return std::numeric_limits<Real>::infinity();
    const Real underflowLimit = std::numeric_limits<Real>::min_exponent == std::numeric_limits<double>::min_exponent
                                    ? Real{745}
                                    : -std::log(std::numeric_limits<Real>::denorm_min());
    if (x > underflowLimit)
        return 0.0;
    Real ri = 0, rk = 0;
    detail::besselIK(nu, x, ri, rk);
    return rk;
}

template <std::floating_point Real>
inline Real riemannZeta(Real s)
{
    if (!(s > 1.0))
        throw std::domain_error("Riemann zeta function: argument should be greater than 1, but it's equal to " +
                                std::to_string(s));
    if (s > 60.0)
        return Real{1} + std::pow(Real{2}, -s) + std::pow(Real{3}, -s);

    /// Euler-Maclaurin summation with N terms and Bernoulli corrections up to B_20
    static constexpr int N = 20;
    static constexpr Real B2jOverFact[] = {Real{1} / 12,
                                           -Real{1} / 720,
                                           Real{1} / 30240,
                                           -Real{1} / 1209600,
                                           Real{1} / 47900160,
                                           -Real{691} / Real{1307674368000},
                                           Real{1} / Real{74724249600},
                                           -Real{3617} / Real{10670622842880000},
                                           Real{43867} / Real{5109094217170944000},
                                           -Real{174611} / static_cast<Real>(802857662698291200000.0L)};
    Real sum = 0;
    for (int k = 1; k < N; ++k)
        sum += std::pow(static_cast<Real>(k), -s);
    Real nPow = std::pow(static_cast<Real>(N), -s);
    sum += N * nPow / (s - 1) + Real{0.5} * nPow;
    Real rising = s, term = nPow / N;
    for (int j = 0; j < 10; ++j) {
        sum += B2jOverFact[j] * rising * term;
        rising *= (s + 2 * j + 1) * (s + 2 * j + 2);
        term /= static_cast<Real>(N) * N;
    }
    return sum;
}

} // namespace portable

template <std::floating_point Real>
inline Real besselI(std::type_identity_t<Real> nu, Real x)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::cyl_bessel_i(nu, x);
#else
    return portable::besselI(nu, x);
#endif
}

template <std::floating_point Real>
inline Real besselK(std::type_identity_t<Real> nu, Real x)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::cyl_bessel_k(nu, x);
#else
    return portable::besselK(nu, x);
#endif
}

template <std::floating_point Real>
inline Real riemannZeta(Real s)
{
#ifdef RANDLIB_HAS_STD_SPECIAL_FUNCTIONS
    return std::riemann_zeta(s);
#else
    return portable::riemannZeta(s);
#endif
}

template <std::floating_point Real>
inline Real lgammaCorrection(Real x)
{
    if (x < 10.0)
        throw std::domain_error("Stirling correction: argument should be at least 10, but it's equal to " +
                                std::to_string(x));
    static constexpr Real c[] = {Real{1} / 12,   -Real{1} / 360,      Real{1} / 1260, -Real{1} / 1680,
                                 Real{1} / 1188, -Real{691} / 360360, Real{1} / 156};
    Real xInv = 1 / x, xInvSq = xInv * xInv;
    Real sum = c[6];
    for (int i = 5; i >= 0; --i)
        sum = sum * xInvSq + c[i];
    return sum * xInv;
}

} // namespace math

} // namespace randlib
