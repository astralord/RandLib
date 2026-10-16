#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <functional>
#include <complex>
#include <vector>
#include <utility>
#include <type_traits>
#include <algorithm>
#include <numeric>

#include <randlib/math/constants.hpp>

namespace randlib
{

namespace math
{

template <int TABLE_SIZE>
constexpr std::array<LongDoublePair, TABLE_SIZE>
createZiggurat(const std::function<LongDoublePair(const LongDoublePair &)> &funPtr, LongDoublePair first,
               LongDoublePair second);

namespace detail
{

/**
 * @param a
 * @param b
 * @param eps
 * @return |a - b| < eps * max(a, b)
 */
template <typename RealType>
bool areClose(RealType a, RealType b, RealType eps = 1e-6)
{
    if (a == b)
        return true;
    RealType fa = std::fabs(a);
    RealType fb = std::fabs(b);
    return std::fabs(b - a) < eps * std::max(fa, fb);
}

} // namespace detail

/**
 * @param x
 * @return sign of x
 */
template <std::floating_point Real>
int sign(Real x);

/**
 * @param x
 * @return safe atan(x)
 */
template <std::floating_point Real>
Real atan(Real x);

/**
 * @param x
 * @return log(1 + exp(x))
 */
template <std::floating_point Real>
Real softplus(Real x);

/**
 * @param x
 * @return log(1 - exp(x))
 */
template <std::floating_point Real>
Real log1mexp(Real x);

/**
 * @param x
 * @return log(exp(x) - 1)
 */
template <std::floating_point Real>
Real logexpm1l(Real x);

/**
 * @param x
 * @return log(2 - exp(x))
 */
template <std::floating_point Real>
Real log2mexp(Real x);

/**
 * @param p
 * @return inverse error function: such x that erf(x) = p
 */
template <std::floating_point Real>
Real erfinv(Real p);

/**
 * @param p
 * @return inverse complementary error function: such x that erfc(x) = p
 */
template <std::floating_point Real>
Real erfcinv(Real p);

/**
 * @param x
 * @return scaled complementary error function exp(x^2) erfc(x)
 */
template <std::floating_point Real>
Real erfcx(Real x);

/**
 * @brief xexpxsqerfc
 * @param x
 * @return x * exp(x^2) * erfc(x)
 */
template <std::floating_point Real>
Real xexpxsqerfc(Real x);

/**
 * @param exponent
 * @param number
 * @return sum_{i=1}^{number} i^{-exponent}
 */
template <std::floating_point Real>
Real harmonicNumber(Real exponent, int number);

/**
 * Calculates logarithm of modified Bessel function of the 1st kind
 * @param nu
 * @param x
 * @return log(I_ν(x))
 */
template <std::floating_point Real>
Real logBesselI(std::type_identity_t<Real> nu, Real x);

/**
 * Calculates logarithm of modified Bessel function of the 2nd kind
 * @param nu
 * @param x
 * @return log(K_ν(x))
 */
template <std::floating_point Real>
Real logBesselK(std::type_identity_t<Real> nu, Real x);

/**
 * @param x
 * @param epsilon
 * @return W0 branch of Lambert W function
 */
template <std::floating_point Real>
Real W0Lambert(Real x, std::type_identity_t<Real> epsilon = Real{1e-11});

/**
 * @param x
 * @param epsilon
 * @return W-1 branch of Lambert W function
 */
template <std::floating_point Real>
Real Wm1Lambert(Real x, std::type_identity_t<Real> epsilon = Real{1e-11});

/**
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @param logX log(x)
 * @param logY log(y)
 * @return 1 - Marcum Q-function
 */
template <std::floating_point Real>
Real MarcumP(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y, std::type_identity_t<Real> sqrtX,
             std::type_identity_t<Real> sqrtY, std::type_identity_t<Real> logX, std::type_identity_t<Real> logY);

/**
 * @param mu
 * @param x
 * @param y
 * @return 1 - Marcum Q-function
 */
template <std::floating_point Real>
Real MarcumP(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y);

/**
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @param logX log(x)
 * @param logY log(y)
 * @return Marcum Q-function
 */
template <std::floating_point Real>
Real MarcumQ(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y, std::type_identity_t<Real> sqrtX,
             std::type_identity_t<Real> sqrtY, std::type_identity_t<Real> logX, std::type_identity_t<Real> logY);

/**
 * @param mu
 * @param x
 * @param y
 * @return Marcum Q-function
 */
template <std::floating_point Real>
Real MarcumQ(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y);
} // namespace math

} // namespace randlib

#include <randlib/math/gamma.hpp>
#include <randlib/math/special_functions.hpp>
#include <randlib/math/beta.hpp>
#include <randlib/math/quadrature.hpp>
#include <randlib/math/roots.hpp>
#include <randlib/math/optimize.hpp>

namespace randlib
{

namespace math
{

template <std::floating_point Real>
inline int sign(Real x)
{
    return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

template <std::floating_point Real>
inline Real atan(Real x)
{
    /// For small absolute values we use standard technique
    /// Otherwise we use relation
    /// atan(x) = +/-π/2 - atan(1/x)
    /// to avoid numeric problems
    if (x == 0.0)
        return 0.0;
    if (x > 1.0)
        return randlib::consts::half_pi - std::atan(1.0 / x);
    return (x < -1.0) ? -randlib::consts::half_pi - std::atan(1.0 / x) : std::atan(x);
}

template <std::floating_point Real>
inline Real softplus(Real x)
{
    if (x < 20.0)
        return std::log1p(std::exp(x));
    return (x < 35.0) ? x + std::exp(-x) : x;
}

template <std::floating_point Real>
inline Real log1mexp(Real x)
{
    return (x < -randlib::consts::ln2) ? std::log1p(-std::exp(x)) : std::log(-std::expm1(x));
}

template <std::floating_point Real>
inline Real logexpm1l(Real x)
{
    if (x < 20.0)
        return std::log(std::expm1(x));
    return (x < 35.0) ? x - std::exp(-x) : x;
}

template <std::floating_point Real>
inline Real log2mexp(Real x)
{
    return std::log1p(-std::expm1(x));
}

namespace detail
{

template <std::floating_point Real>
inline Real erfinvChebyshevSeries(Real x, long double t, const long double *array, int size)
{
    /// We approximate inverse erf via Chebyshev polynomials
    long double Tn = t, Tnm1 = 1, sum = 0.0;
    for (int i = 1; i != size; ++i) {
        sum += array[i] * Tn;
        /// calculate next Chebyshev polynomial
        long double temp = Tn;
        Tn *= 2 * t;
        Tn -= Tnm1;
        Tnm1 = temp;
    }
    return static_cast<Real>(static_cast<long double>(x) * (array[0] + sum));
}

template <std::floating_point Real>
inline Real erfinvAux1(Real beta)
{
    /// |1 - p| < 5e-16
    static constexpr long double D5 = -9.199992358830151031278420l, D6 = 2.794990820124599493768426l;
    static constexpr int muSize = 26;
    static constexpr long double mu[muSize] = {.9885750640661893136460358l,
                                               .0108577051845994776160281l,
                                               -.0017511651027627952594825l,
                                               .0000211969932065633437984l,
                                               .0000156648714042435087911l,
                                               -.05190416869103124261e-5l,
                                               -.00371357897426717780e-5l,
                                               .00012174308662357429e-5l,
                                               -.00001768115526613442e-5l,
                                               -.119372182556161e-10l,
                                               .003802505358299e-10l,
                                               -.000660188322362e-10l,
                                               -.000087917055170e-10l,
                                               -.3506869329e-15l,
                                               -.0697221497e-15l,
                                               -.0109567941e-15l,
                                               -.0011536390e-15l,
                                               -.0000263938e-15l,
                                               .05341e-20l,
                                               -.22610e-20l,
                                               .09552e-20l,
                                               -.05250e-20l,
                                               .02487e-20l,
                                               -.01134e-20l,
                                               .00420e-20l};
    return erfinvChebyshevSeries(beta, D5 / std::sqrt(beta) + D6, mu, muSize);
}

template <std::floating_point Real>
inline Real erfinvAux2(Real beta)
{
    /// 5e-16 < |1 - p| < 0.025
    static constexpr long double D3 = -0.5594576313298323225436913l, D4 = 2.287915716263357638965891l;
    static constexpr int deltaSize = 38;
    static constexpr long double delta[deltaSize] = {.9566797090204925274526373l,
                                                     -.0231070043090649036999908l,
                                                     -.0043742360975084077333218l,
                                                     -.0005765034226511854809364l,
                                                     -.0000109610223070923931242l,
                                                     .0000251085470246442787982l,
                                                     .0000105623360679477511955l,
                                                     .27544123300306391503e-5l,
                                                     .04324844983283380689e-5l,
                                                     -.00205303366552086916e-5l,
                                                     -.00438915366654316784e-5l,
                                                     -.00176840095080881795e-5l,
                                                     -.00039912890280463420e-5l,
                                                     -.00001869324124559212e-5l,
                                                     .00002729227396746077e-5l,
                                                     .00001328172131565497e-5l,
                                                     .318342484482286e-10l,
                                                     .016700607751926e-10l,
                                                     -.020364649611537e-10l,
                                                     -.009648468127965e-10l,
                                                     -.002195672778128e-10l,
                                                     -.000095689813014e-10l,
                                                     .000137032572230e-10l,
                                                     .000062538505417e-10l,
                                                     .000014584615266e-10l,
                                                     .1078123993e-15l,
                                                     -.0709229988e-15l,
                                                     -.0391411775e-15l,
                                                     -.0111659209e-15l,
                                                     -.0015770366e-15l,
                                                     .0002853149e-15l,
                                                     .0002716662e-15l,
                                                     .0000176835e-15l,
                                                     .09828e-20l,
                                                     .20464e-20l,
                                                     .08020e-20l,
                                                     .01650e-20l};
    return erfinvChebyshevSeries(beta, D3 * beta + D4, delta, deltaSize);
}

template <std::floating_point Real>
inline Real erfinvAux3(Real beta)
{
    /// 0.8 < p < 0.975
    static constexpr long double D1 = -1.548813042373261659512742l, D2 = 2.565490123147816151928163l;
    static constexpr int lambdaSize = 27;
    static constexpr long double lambda[lambdaSize] = {.9121588034175537733059200l,
                                                       -.0162662818676636958546661l,
                                                       .0004335564729494453650589l,
                                                       .0002144385700744592065205l,
                                                       .26257510757648130176e-5l,
                                                       -.30210910501037969912e-5l,
                                                       -.00124060618367572157e-5l,
                                                       .00624066092999917380e-5l,
                                                       -.00005401247900957858e-5l,
                                                       -.00014232078975315910e-5l,
                                                       .343840281955305e-10l,
                                                       .335848703900138e-10l,
                                                       -.014584288516512e-10l,
                                                       -.008102174258833e-10l,
                                                       .000525324085874e-10l,
                                                       .000197115408612e-10l,
                                                       -.000017494333828e-10l,
                                                       -.4800596619e-15l,
                                                       .0557302987e-15l,
                                                       .0116326054e-15l,
                                                       -.0017262489e-15l,
                                                       -.0002784973e-15l,
                                                       .0000524481e-15l,
                                                       .65270e-20l,
                                                       -.15707e-20l,
                                                       -.01475e-20l,
                                                       .00450e-20l};
    return erfinvChebyshevSeries(beta, D1 * beta + D2, lambda, lambdaSize);
}

template <std::floating_point Real>
inline Real erfinvAux4(Real p)
{
    /// 0 < p < 0.8
    static constexpr int xiSize = 39;
    static constexpr long double xi[xiSize] = {.9928853766189408231495800l,
                                               .1204675161431044864647846l,
                                               .0160781993420999447267039l,
                                               .0026867044371623158279591l,
                                               .0004996347302357262947170l,
                                               .0000988982185991204409911l,
                                               .0000203918127639944337340l,
                                               .43272716177354218758e-5l,
                                               .09380814128593406758e-5l,
                                               .02067347208683427411e-5l,
                                               .00461596991054300078e-5l,
                                               .00104166797027146217e-5l,
                                               .00023715009995921222e-5l,
                                               .00005439284068471390e-5l,
                                               .00001255489864097987e-5l,
                                               .291381803663201e-10l,
                                               .067949421808797e-10l,
                                               .015912343331569e-10l,
                                               .003740250585245e-10l,
                                               .000882087762421e-10l,
                                               .000208650897725e-10l,
                                               .000049488041039e-10l,
                                               .000011766394740e-10l,
                                               .2803855725e-15l,
                                               .0669506638e-15l,
                                               .0160165495e-15l,
                                               .0038382583e-15l,
                                               .0009212851e-15l,
                                               .0002214615e-15l,
                                               .0000533091e-15l,
                                               .0000128488e-15l,
                                               .31006e-20l,
                                               .07491e-20l,
                                               .01812e-20l,
                                               .00439e-20l,
                                               .00106e-20l,
                                               .00026e-20l,
                                               .00006e-20l,
                                               .00002e-20l};
    return erfinvChebyshevSeries(p, p * p / 0.32 - 1.0, xi, xiSize);
}

} // namespace detail

template <std::floating_point Real>
inline Real erfinv(Real p)
{
    /// Consider special cases
    if (p < 0.0)
        return -erfinv(-p);
    if (p > 1.0)
        throw std::invalid_argument("Argument p should be in interval [-1, 1]");
    if (p == 1.0)
        return INFINITY;
    if (p == 0.0)
        return 0.0;
    if (p < 0.8)
        return detail::erfinvAux4(p);
    /// Handle tails
    Real beta = std::sqrt(-std::log1p(-p * p));
    if (p < 0.9975)
        return detail::erfinvAux3(beta);
    return (1.0 - p < 5e-16) ? detail::erfinvAux1(beta) : detail::erfinvAux2(beta);
}

template <std::floating_point Real>
inline Real erfcinv(Real p)
{
    /// Consider special cases
    if (p > 1.0)
        return -erfcinv(2.0 - p);
    if (p == 0.0)
        return INFINITY;
    if (p == 1.0)
        return 0.0;
    if (p > 0.2)
        return detail::erfinvAux4(1.0 - p);
    Real pSq = p * p, p2 = 2 * p;
    Real beta = std::sqrt(-std::log(p2 - pSq));
    if (p > 0.0025)
        return detail::erfinvAux3(beta);
    return (p > 5e-16) ? detail::erfinvAux2(beta) : detail::erfinvAux1(beta);
}

template <std::floating_point Real>
inline Real xexpxsqerfc(Real x)
{
    static constexpr int MAX_X = 10;
    static constexpr int N = 10;
    if (x < MAX_X) {
        Real y = x * x;
        y += std::log(std::erfc(x));
        return x * std::exp(y);
    }
    Real log2xSq = randlib::consts::ln2 + 2 * std::log(x);
    Real sum = 0.0;
    for (int n = 1; n != N; ++n) {
        Real add = math::ldfact(2 * n - 1);
        add -= n * log2xSq;
        add = std::exp(add);
        sum += (n & 1) ? -add : add;
    }
    return (1.0 + sum) / randlib::consts::sqrtpi;
}

template <std::floating_point Real>
inline Real erfcx(Real x)
{
    if (x < 1.0)
        return std::exp(x * x) * std::erfc(x);
    return xexpxsqerfc(x) / x;
}

template <std::floating_point Real>
inline Real harmonicNumber(Real exponent, int number)
{
    if (number < 1)
        return 0;
    if (exponent == 1)
        return static_cast<Real>(randlib::consts::egamma) + digamma(static_cast<Real>(number + 1));
    if (exponent == 2)
        return static_cast<Real>(randlib::consts::pi_sq) / 6 - trigamma(static_cast<Real>(number + 1));
    Real res = 1.0;
    for (int i = 2; i <= number; ++i)
        res += std::pow(i, -exponent);
    return res;
}

template <std::floating_point Real>
inline Real logBesselI(std::type_identity_t<Real> nu, Real x)
{
    if (x < 0) {
        Real roundNu = std::round(nu);
        bool nuIsInt = detail::areClose(nu, roundNu);
        if (nuIsInt) {
            int nuInt = roundNu;
            return (nuInt % 2) ? NAN : logBesselI(nu, -x);
        }
        return -INFINITY;
    }

    if (x == 0) {
        if (nu == 0)
            return 0.0;
        Real roundNu = std::round(nu);
        bool nuIsInt = detail::areClose(nu, roundNu);
        return (nu > 0 || nuIsInt) ? -INFINITY : INFINITY;
    }

    if (std::fabs(nu) == 0.5) {
        /// log(sinh(x)) or log(cosh(x))
        Real y = x - 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + std::log(x));
        y += (nu > 0) ? math::log1mexp(-2 * x) : math::softplus(-2 * x);
        return y;
    }

    if (nu < 0) {
        /// I(ν, x) = I(-ν, x) - 2 / π sin(πν) K(-ν, x), and sin(πν) = -sin(π|ν|) for ν < 0
        Real besseli = math::besselI(-nu, x);
        Real sinPiAbsNu = std::sin(-randlib::consts::pi * nu);
        Real y = 0;
        if (sinPiAbsNu == 0 || detail::areClose(nu, std::round(nu)))
            y = besseli;
        else {
            Real besselk = math::besselK(-nu, x);
            y = besseli + randlib::consts::two_inv_pi * sinPiAbsNu * besselk;
        }
        return (y <= 0) ? -INFINITY : std::log(y);
    }

    Real besseli = math::besselI(nu, x); // TODO: expand Hankel asymptotic expansions
    return std::isfinite(besseli) ? std::log(besseli)
                                  : x - 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + std::log(x));
}

template <std::floating_point Real>
inline Real logBesselK(std::type_identity_t<Real> nu, Real x)
{
    if (nu < 0.0)
        return NAN; /// K(-ν, x) = -K(ν, x) < 0

    if (x == 0.0)
        return INFINITY;

    Real besselk = 0;
    if (nu == 0.5 || (besselk = math::besselK(nu, x)) == 0)
        return 0.5 * (randlib::consts::lnpi - randlib::consts::ln2 - std::log(x)) - x;

    if (!std::isfinite(besselk)) // TODO: expand Hankel asymptotic expansions
        return (nu == 0) ? std::log(-std::log(x))
                         : std::lgamma(nu) - static_cast<Real>(randlib::consts::ln2) - nu * std::log(Real{0.5} * x);

    return std::log(besselk);
}

/**
 * @param x
 * @param w0
 * @param epsilon
 * @return
 */
namespace detail
{

template <std::floating_point Real>
inline Real WLambert(Real x, Real w0, Real epsilon)
{
    Real w = w0;
    Real step = 0;
    do {
        Real ew = std::exp(w);
        Real wew = w * ew;
        Real numerator1 = wew - x;
        Real wp1 = w + 1;
        Real denominator1 = ew * wp1;
        Real numerator2 = (w + 2) * numerator1;
        Real denominator2 = 2 * wp1;
        step = numerator2 / denominator2;
        step = numerator1 / (denominator1 - step);
        w -= step;
    } while (std::fabs(step) > epsilon);
    return w;
}

} // namespace detail

template <std::floating_point Real>
inline Real W0Lambert(Real x, std::type_identity_t<Real> epsilon)
{
    Real w = 0;
    if (x < -randlib::consts::inv_e)
        throw std::invalid_argument("Argument x should be greater than -1/e, but it's equal to " + std::to_string(x));
    if (x > 10) {
        Real logX = std::log(x);
        Real loglogX = std::log(logX);
        w = logX - loglogX;
    }
    return detail::WLambert(x, w, epsilon);
}

template <std::floating_point Real>
inline Real Wm1Lambert(Real x, std::type_identity_t<Real> epsilon)
{
    Real w = -2;
    if (x < -randlib::consts::inv_e || x > 0)
        throw std::invalid_argument(
            "Argument x should be greater than -1/e and smaller or equal to 0, but it's equal to " + std::to_string(x));
    if (x > -0.1) {
        Real logmX = std::log(-x);
        Real logmlogmX = std::log(-logmX);
        w = logmX - logmlogmX;
    }
    return detail::WLambert(x, w, epsilon);
}

namespace detail
{

/**
 * @param mu
 * @param x
 * @param y
 * @param logX log(x)
 * @param logY log(y)
 * @return series expansion for Marcum-P function
 */
template <std::floating_point Real>
inline Real MarcumPSeries(Real mu, Real x, Real y, Real logX, Real logY)
{
    /// ~log(2πε) for ε = 1e-16
    static constexpr Real ln2piEps = -35.0;
    Real lgammamu = std::lgamma(mu);
    Real C = lgammamu - ln2piEps + mu;

    /// solving equation f(n) = 0
    /// to find first negleted term
    Real root = std::max(Real{0.5} * (mu * mu + 4 * x * y - mu), Real{1});
    Real logXY = logX + logY;
    if (!math::findRootNewtonSecondOrder<Real>(
            [C, mu, logXY](Real n) {
                Real npmu = n + mu;
                Real logn = std::log(n), lognpmu = std::log(npmu);
                Real first = logn - 2 - logXY;
                first *= n;
                first += npmu * lognpmu;
                first -= C;
                Real second = logn + lognpmu - logXY;
                Real third = 1.0 / n + 1.0 / npmu;
                return Triplet<Real>(first, second, third);
            },
            root))
        /// unexpected return
        throw std::runtime_error("Marcum P function: failure in numerical procedure");

    /// series expansion
    Real sum = 0.0;
    /// sanity check
    int n0 = static_cast<int>(std::max(std::ceil(root), Real{5}));
    Real mpn0 = mu + n0;
    Real P = pgamma(mpn0, y, logY);
    Real diffP = (mpn0 - 1) * logY - y - std::lgamma(mpn0);
    diffP = std::exp(diffP);
    for (int n = n0; n > 0; --n) {
        Real term = n * logX - x - lfact(n);
        Real mupnm1 = mu + n - 1;
        term = std::exp(term) * P;
        sum += term;
        if (n % 5 == 0) {
            /// every 5 iterations we recalculate P and diffP
            /// in order to achieve enough accuracy
            P = pgamma(mupnm1, y, logY);
            diffP = (mupnm1 - 1) * logY - y - std::lgamma(mupnm1);
            diffP = std::exp(diffP);
        } else {
            /// otherwise we use recurrent relations
            P += diffP;
            diffP *= mupnm1 / y;
        }
    }
    /// add the last 0-term
    Real lastTerm = lpgamma(mu, y, logY) - x;
    lastTerm = std::exp(lastTerm);
    sum += lastTerm;
    return sum;
}

template <std::floating_point Real>
inline Real MarcumQSeries(Real mu, Real x, Real y, Real logX, Real logY)
{
    if (x <= 0)
        return qgamma(mu, y, logY);
    const Real saddle = Real{0.5} * (std::sqrt(mu * mu + 4 * x * y) - mu);
    const int center = static_cast<int>(std::floor(std::max(x, saddle)));
    Real shape = mu + center;
    Real weight = std::exp(-x + center * logX - std::lgamma(center + 1.0));
    Real q = qgamma(shape, y, logY);
    Real delta = std::exp(shape * logY - y - std::lgamma(shape + 1));
    Real sum = weight * q;

    Real lowerWeight = weight;
    Real previousTerm = sum;
    int decreasing = 0;
    for (int n = center; n > 0; --n) {
        lowerWeight *= n / x;
        const Real lowerQ = qgamma(mu + n - 1, y, logY);
        const Real term = lowerWeight * lowerQ;
        sum += term;
        decreasing = (term <= previousTerm) ? decreasing + 1 : 0;
        if (decreasing >= 8 && term <= 4 * std::numeric_limits<Real>::epsilon() * sum)
            break;
        previousTerm = term;
    }

    decreasing = 0;
    previousTerm = weight * q;
    for (int n = center + 1; n < 1000000; ++n) {
        weight *= x / n;
        q += delta;
        shape += 1;
        delta *= y / shape;
        const Real term = weight * std::min(Real{1}, q);
        sum += term;
        decreasing = (term <= previousTerm) ? decreasing + 1 : 0;
        if (decreasing >= 8 && term <= 4 * std::numeric_limits<Real>::epsilon() * sum)
            break;
        previousTerm = term;
    }
    return std::clamp(sum, Real{0}, Real{1});
}

/**
 * @param mu
 * @param x
 * @param y
 * @param sqrtX √x
 * @param sqrtY √y
 * @return asymptotic expansion for Marcum-P function for large x*y
 */
template <std::floating_point Real>
inline Real MarcumPAsymptoticForLargeXY(Real mu, Real x, Real y, Real sqrtX, Real sqrtY)
{
    Real xi = 2 * sqrtX * sqrtY;
    Real sigma = (y + x) / xi - 1;
    Real sigmaXi = sigma * xi;
    Real rho = sqrtY / sqrtX;
    Real aux = std::erfc(sqrtX - sqrtY);
    Real Phi = (sigma == 0) ? Real{0} : sign(x - y) * std::sqrt(static_cast<Real>(randlib::consts::pi) / sigma) * aux;
    Real Psi0 = aux / std::sqrt(rho);
    Real logXi = randlib::consts::ln2 + 0.5 * std::log(x * y);
    /// sanity check
    int n0 = static_cast<int>(std::max(std::ceil(sigmaXi), Real{7}));
    Real sum = 0.0;
    Real A1 = 1.0, A2 = 1.0;
    for (int n = 1; n <= n0; ++n) {
        /// change φ
        Real nmHalf = n - 0.5;
        Phi *= -sigma;
        Real temp = -sigmaXi - nmHalf * logXi;
        Phi += std::exp(temp);
        Phi /= nmHalf;
        /// calculate A(μ-1) and A(μ)
        Real coef1 = (nmHalf - mu) * (nmHalf + mu);
        Real coef2 = (nmHalf - mu + 1) * (nmHalf + mu - 1);
        Real denom = -2 * n;
        A1 *= coef1 / denom;
        A2 *= coef2 / denom;
        /// compute term ψ and add it to the sum
        Real Psi = Phi * (A2 - A1 / rho);
        sum += (n & 1) ? Psi : -Psi;
    }
    sum /= randlib::consts::sqrt2pi;
    sum += Psi0;
    sum *= 0.5 * std::pow(rho, mu);
    return sum;
}

/**
 * @param mu
 * @param x
 * @param y
 * @param logX precomputed log(x)
 * @param logY precomputed log(y)
 * @return
 */
template <std::floating_point Real>
inline Real MarcumPForMuSmallerThanOne(Real mu, Real x, Real y, Real logX, Real logY)
{
    // TODO: check Krishnamoorthy paper for alternative representation

    /// in this case we use numerical integration,
    /// however we have singularity point at 0,
    /// so we get rid of it by subtracting the function
    /// which has the same behaviour at this point
    Real aux = x + mu * static_cast<Real>(randlib::consts::ln2) + std::lgamma(mu);
    Real log2x = randlib::consts::ln2 + logX;
    Real I = (randlib::consts::ln2 + logY) * mu - aux;
    I = std::exp(I) / mu;

    Real mum1 = mu - 1.0;
    I += math::integrate(
        [x, log2x, mum1, aux](Real t) {
            if (t <= 0)
                return Real{0};
            /// Calculate log of leveling factor
            Real log2T = std::log(2 * t);
            Real exponent = mum1 * log2T;
            exponent -= aux;

            /// Calculate log(2f(t))
            Real logBessel = math::logBesselI(mum1, 2 * std::sqrt(x * t));
            Real z = mum1 * (log2T - log2x);
            Real log2F = 0.5 * z - t - x + logBessel;

            /// Return difference f(t) - factor
            return std::exp(log2F) - 2 * std::exp(exponent);
        },
        0, y);
    return I;
}

template <std::floating_point Real>
inline Real MarcumQIntegrand(Real theta, Real xi, Real sqrt1pXiSq, Real mu, Real y)
{
    Real sinTheta = std::sin(theta);
    Real theta_sinTheta = theta / sinTheta;
    Real rho = std::hypot(theta_sinTheta, xi);
    Real theta_sinThetapRho = theta_sinTheta + rho;
    Real r = 0.5 * theta_sinThetapRho / y;
    Real cosTheta = std::cos(theta);
    Real psi = cosTheta * rho - sqrt1pXiSq;
    Real frac = theta_sinThetapRho / (1.0 + sqrt1pXiSq);
    psi -= std::log(frac);
    Real numerator = (sinTheta - theta * cosTheta) / (sinTheta * rho);
    numerator += cosTheta - r;
    numerator *= r;
    Real denominator = r * r - 2.0 * r * cosTheta + 1.0;
    Real f = numerator / denominator;
    return std::exp(mu * psi) * f;
}

template <std::floating_point Real>
inline Real MarcumQIntergralRepresentation(Real mu, Real x, Real y, Real sqrtX, Real sqrtY)
{
    Real xi = 2.0 * sqrtX * sqrtY / mu;
    Real sqrt1pXiSq = std::sqrt(1.0 + xi * xi);
    Real yPrime = y / mu;
    Real s0 = 0.5 * (1.0 + sqrt1pXiSq) / yPrime;
    Real phi = x / s0 + y * s0 - std::log(s0) * mu;
    auto integrand = [=](Real theta) { return MarcumQIntegrand(theta, xi, sqrt1pXiSq, mu, yPrime); };
    Real integral =
        math::integrate(integrand, -static_cast<Real>(randlib::consts::pi), static_cast<Real>(randlib::consts::pi));
    return 0.5 * std::exp(-x - y + phi) / randlib::consts::pi * integral;
}

} // namespace detail

template <std::floating_point Real>
inline Real MarcumP(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y,
                    std::type_identity_t<Real> sqrtX, std::type_identity_t<Real> sqrtY, std::type_identity_t<Real> logX,
                    std::type_identity_t<Real> logY)
{
    /* 1 - ok
     * 2 - ok
     * 3 - no
     * 4 - no
     * 5 - not yet
     * */
    if (x < 0.0 || y <= 0.0)
        return 0.0;
    if (std::isinf(y))
        return 1.0;
    if (std::isinf(x))
        return 0.0;

    if (x < 30)
        return detail::MarcumPSeries(mu, x, y, logX, logY);

    Real xi = 2 * sqrtX * sqrtY;
    if (xi > 30 && mu * mu < 2 * xi)
        return detail::MarcumPAsymptoticForLargeXY(mu, x, y, sqrtX, sqrtY);

    Real temp = std::sqrt(4 * x + 2 * mu);
    Real f1 = x + mu - temp, f2 = x + mu + temp;
    if (y > f1 && y < f2)
        // IF mu > 135
        return 1.0 - detail::MarcumQIntergralRepresentation(mu, x, y, sqrtX, sqrtY);

    // TODO: implement the rest techniques

    Real mum1 = mu - 1;
    return math::integrate(
        [mum1, logX, x](Real t) {
            if (t < 0.0)
                return Real{0};
            if (t == 0.0)
                return (mum1 == 0) ? Real{0.5} * std::exp(-x) : Real{0};
            Real logBesseli = math::logBesselI(mum1, 2 * std::sqrt(x * t));
            Real z = 0.5 * mum1 * (std::log(t) - logX);
            Real h = t + x;
            return std::exp(logBesseli + z - h);
        },
        0, y);
}

template <std::floating_point Real>
inline Real MarcumP(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y)
{
    Real sqrtX = std::sqrt(x), sqrtY = std::sqrt(y);
    Real logX = std::log(x), logY = std::log(y);
    return MarcumP(mu, x, y, sqrtX, sqrtY, logX, logY);
}

template <std::floating_point Real>
inline Real MarcumQ(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y,
                    std::type_identity_t<Real> sqrtX, std::type_identity_t<Real> sqrtY, std::type_identity_t<Real> logX,
                    std::type_identity_t<Real> logY)
{
    if (x < 0 || y <= 0)
        return 1;
    if (std::isinf(y))
        return 0;
    if (std::isinf(x))
        return 1;
    const Real saddle = Real{0.5} * (std::hypot(mu, 2 * sqrtX * sqrtY) - mu);
    if (!std::isfinite(saddle) || std::max(x, saddle) >= std::numeric_limits<int>::max() - 1)
        return 1 - MarcumP(mu, x, y, sqrtX, sqrtY, logX, logY);
    return detail::MarcumQSeries(mu, x, y, logX, logY);
}

template <std::floating_point Real>
inline Real MarcumQ(std::type_identity_t<Real> mu, Real x, std::type_identity_t<Real> y)
{
    Real sqrtX = std::sqrt(x), sqrtY = std::sqrt(y);
    Real logX = std::log(x), logY = std::log(y);
    return MarcumQ(mu, x, y, sqrtX, sqrtY, logX, logY);
}

} // namespace math

} // namespace randlib
