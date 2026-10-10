#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include <algorithm>
#include <functional>
#include <vector>

template < typename RealType >
StableDistribution<RealType>::StableDistribution(double exponent, double skewness, double scale, double location)
{
    SetParameters(exponent, skewness, scale, location);
}

template < typename RealType >
SUPPORT_TYPE StableDistribution<RealType>::SupportType() const
{
    if (alpha < 1) {
        if (beta == 1)
            return RIGHTSEMIFINITE_T;
        if (beta == -1)
            return LEFTSEMIFINITE_T;
    }
    return INFINITE_T;
}

template < typename RealType >
RealType StableDistribution<RealType>::MinValue() const
{
    return (alpha < 1 && beta == 1) ? mu : -INFINITY;
}

template < typename RealType >
RealType StableDistribution<RealType>::MaxValue() const
{
    return (alpha < 1 && beta == -1) ? mu : INFINITY;
}

template < typename RealType >
void StableDistribution<RealType>::parametersVerification(double exponent, double skewness, double scale)
{
    if (!(exponent >= 0.1 && exponent <= 2.0))
        throw std::invalid_argument("Stable distribution: exponent should be in the interval [0.1, 2], but it's equal to "
                                    + std::to_string(exponent));
    if (!(std::fabs(skewness) <= 1.0))
        throw std::invalid_argument("Stable distribution: skewness should be in the interval [-1, 1], but it's equal to "
                                    + std::to_string(skewness));
    if (!(scale > 0.0))
        throw std::invalid_argument("Stable distribution: scale should be positive, but it's equal to "
                                    + std::to_string(scale));
}

template < typename RealType >
void StableDistribution<RealType>::checkDensitySupport() const
{
    bool exponentNearOne = alpha != 1.0 && std::fabs(alpha - 1.0) < 0.01 && beta != 0.0;
    bool skewnessNearZero = alpha == 1.0 && beta != 0.0 && std::fabs(beta) < 0.01;
    if (exponentNearOne || skewnessNearZero)
        throw randlib::not_implemented("Stable distribution: the density and distribution function are not yet implemented "
                                       "for exponent within 0.01 of 1 and non-zero skewness, or for exponent 1 and skewness "
                                       "within 0.01 of 0; sampling and the characteristic function are available");
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForNormal()
{   
    distributionType = NORMAL;
    pdfCoef = randlib::consts::ln2 + logGamma + 0.5 * randlib::consts::lnpi;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForCauchy()
{   
    distributionType = CAUCHY;
    pdfCoef = -logGamma - randlib::consts::lnpi;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForLevy()
{   
    distributionType = LEVY;
    pdfCoef = logGamma - randlib::consts::ln2 - randlib::consts::lnpi;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForUnityExponent()
{   
    distributionType = UNITY_EXPONENT;
    pdfCoef = 0.5 / (gamma * std::fabs(beta));
    pdftailBound = 0; // not in the use for now
    /// In S1, X = γZ + (2/π)βγ log γ + μ for α = 1, which is where log γ enters
    logGammaPi_2 = randlib::consts::lnpi - randlib::consts::ln2 - logGamma;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForGeneralExponent()
{   
    distributionType = GENERAL;
    if (beta != 0.0) {
        zeta = -beta * std::tan(randlib::consts::half_pi * alpha);
        omega = 0.5 * alphaInv * std::log1pl(zeta * zeta);
        xi = alphaInv * RandMath::atan(-zeta);
    }
    else {
        zeta = omega = xi = 0.0;
    }
    pdfCoef = randlib::consts::inv_pi * std::fabs(alpha_alpham1) / gamma;
    pdftailBound = 3.0 / (1.0 + alpha) * randlib::consts::ln10;
    cdftailBound = 3.0 / alpha * randlib::consts::ln10;
    /// define boundaries of region near 0, where we use series expansion
    if (alpha <= ALMOST_TWO) {
        seriesZeroParams.first = std::round(std::min(alpha * alpha * 40 + 1, 10.0));
        seriesZeroParams.second = -(alphaInv * 1.5 + 0.5) * randlib::consts::ln10; /// corresponds to interval [10^(-15.5), ~0.056]
    }
    else {
        seriesZeroParams.first = 85;
        seriesZeroParams.second = randlib::consts::ln2 + randlib::consts::ln3;///< corresponds to 6
    }
}

template < typename RealType >
void StableDistribution<RealType>::SetParameters(double exponent, double skewness, double scale, double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Stable distribution: location should be finite");
    parametersVerification(exponent, skewness, scale);
    alpha = exponent;
    alphaInv = 1.0 / alpha;
    beta = skewness;
    mu = location;
    gamma = scale;
    logGamma = std::log(gamma);
    alpha_alpham1 = alpha / (alpha - 1.0);

    /// Set id of distribution
    if (alpha == 2.0)
        return setParametersForNormal();
    if (alpha == 1.0)
        return (beta == 0.0) ? setParametersForCauchy() : setParametersForUnityExponent();
    if (alpha == 0.5 && std::fabs(beta) == 1.0)
        return setParametersForLevy();
    return setParametersForGeneralExponent();
}

template < typename RealType >
void StableDistribution<RealType>::SetLocation(double location)
{
    SetParameters(alpha, beta, gamma, location);
}

template < typename RealType >
void StableDistribution<RealType>::SetScale(double scale)
{
    SetParameters(alpha, beta, scale, mu);
}

template < typename RealType >
double StableDistribution<RealType>::pdfNormal(RealType x) const
{
    return std::exp(logpdfNormal(x));
}

template < typename RealType >
double StableDistribution<RealType>::logpdfNormal(RealType x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    y *= y;
    y += pdfCoef;
    return -y;
}

template < typename RealType >
double StableDistribution<RealType>::pdfCauchy(RealType x) const
{
    double y = x - mu;
    y *= y;
    y /= gamma;
    y += gamma;
    return randlib::consts::inv_pi / y;
}

template < typename RealType >
double StableDistribution<RealType>::logpdfCauchy(RealType x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    double xSq = x0 * x0;
    return pdfCoef - std::log1pl(xSq);
}

template < typename RealType >
double StableDistribution<RealType>::pdfLevy(RealType x) const
{
    return (x <= mu) ? 0.0 : std::exp(logpdfLevy(x));
}

template < typename RealType >
double StableDistribution<RealType>::logpdfLevy(RealType x) const
{
    double x0 = x - mu;
    if (x0 <= 0.0)
        return -INFINITY;
    double y = gamma / x0;
    y += 3 * std::log(x0);
    y -= pdfCoef;
    return -0.5 * y;
}

template < typename RealType >
double StableDistribution<RealType>::fastpdfExponentiation(double u)
{
    /// exp(u - exp(u)); for very negative u, exp(u) alone, and for large u it underflows to 0 by itself
    return (u < -700) ? std::exp(u) : std::exp(u - std::exp(u));
}

template < typename RealType >
double StableDistribution<RealType>::pdfShortTailExpansionForUnityExponent(double x) const
{
    if (x > 10)
        return 0.0;
    double xm1 = x - 1.0;
    double y = 0.5 * xm1 - std::exp(xm1);
    y -= 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi);
    return std::exp(y);
}

template < typename RealType >
double StableDistribution<RealType>::limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const
{
    if (theta > 0.0) {
        if (beta > 0.0)
            return BIG_NUMBER;
        return (beta == -1) ? xAdj - 1.0 : -BIG_NUMBER;
    }
    if (beta < 0.0)
        return BIG_NUMBER;
    return (beta == 1) ? xAdj - 1.0 : -BIG_NUMBER;
}

template < typename RealType >
double StableDistribution<RealType>::integrandAuxForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
    if (theta == 0.0)
        return xAdj + randlib::consts::lnpi - randlib::consts::ln2;
    double thetaAdj = (randlib::consts::half_pi + beta * theta) / std::cos(theta);
    double u = std::log(thetaAdj);
    u += thetaAdj * std::sin(theta) / beta;
    return std::isfinite(u) ? u + xAdj : limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
}

template < typename RealType >
double StableDistribution<RealType>::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    return fastpdfExponentiation(u);
}

template < typename RealType >
double StableDistribution<RealType>::rootOfMonotone(const std::function<double (double)> &u, double a, double b)
{
    /// Bisection on the sign of u, which changes exactly once on (a, b)
    bool increasing = u(b) > u(a);
    for (int i = 0; i < 200; ++i) {
        double middle = 0.5 * (a + b);
        if (!(middle > a && middle < b))
            break;
        bool belowRoot = (u(middle) < 0) == increasing;
        (belowRoot ? a : b) = middle;
    }
    return 0.5 * (a + b);
}

template < typename RealType >
double StableDistribution<RealType>::integrateAroundRoot(const std::function<double (double)> &integrand,
                                                         const std::function<double (double)> &u,
                                                         double a, double b, double root)
{
    /// In the tails the integrand is concentrated in a region of width ~1/|u'(root)| around the
    /// root of u, which can be far narrower than [a, b]: a quadrature rule over the whole
    /// interval may not place a single node there. Break the interval at log-spaced distances
    /// from the root, scaled by that width.
    double h = 1e-7 * std::max(std::min(root - a, b - root), 1e-300);
    double slope = (u(root + h) - u(root - h)) / (2 * h);
    double width = (std::isfinite(slope) && slope != 0.0) ? 1.0 / std::fabs(slope) : b - a;
    std::vector<double> points = {a, root, b};
    for (double scale : {1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0, 1000.0}) {
        points.push_back(root - scale * width);
        points.push_back(root + scale * width);
    }
    std::sort(points.begin(), points.end());
    double sum = 0.0, left = a;
    for (double point : points) {
        if (point <= left || point > b)
            continue;
        sum += RandMath::integral(integrand, left, point, 1e-11, 20);
        left = point;
    }
    return sum;
}

template < typename RealType >
double StableDistribution<RealType>::peakForUnityExponent(double xAdj) const
{
    return rootOfMonotone([this, xAdj] (double theta) { return integrandAuxForUnityExponent(theta, xAdj); },
                          -randlib::consts::half_pi, randlib::consts::half_pi);
}

template < typename RealType >
double StableDistribution<RealType>::pdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / gamma;
    double xAdj = -randlib::consts::half_pi * xSt / beta - logGammaPi_2;
    /// The integrand g exp(-g) with g = exp(u) peaks where u = 0, which approaches ±π/2 in the
    /// tails; integrating on both sides of the peak lets the quadrature resolve it
    double theta0 = peakForUnityExponent(xAdj);
    auto integrand = [this, xAdj] (double theta) { return integrandForUnityExponent(theta, xAdj); };
    auto u = [this, xAdj] (double theta) { return integrandAuxForUnityExponent(theta, xAdj); };
    return pdfCoef * integrateAroundRoot(integrand, u, -randlib::consts::half_pi, randlib::consts::half_pi, theta0);
}

template < typename RealType >
double StableDistribution<RealType>::pdfShortTailExpansionForGeneralExponent(double logX) const
{
    double logAlpha = std::log(alpha);
    double log1mAlpha = (alpha < 1) ? std::log1pl(-alpha) : std::log(alpha - 1);
    double temp = logX - logAlpha;
    double y = std::exp(alpha_alpham1 * temp);
    y *= -std::fabs(1.0 - alpha);
    double z = (1.0 - 0.5 * alpha) / (alpha - 1) * temp;
    z -= 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + logAlpha + log1mAlpha);
    return std::exp(y + z);
}

template < typename RealType >
double StableDistribution<RealType>::cdfShortTailExpansion(double logX) const
{
    /// Mass between a totally skewed law (α < 1) and its support boundary at distance x:
    /// exp(-(1 - α) u) / sqrt(2π α (1 - α) u) with u = (x/α)^(-α/(1-α)), to leading order
    double logU = alpha / (alpha - 1.0) * (logX - std::log(alpha));
    double u = std::exp(logU);
    double y = -(1.0 - alpha) * u - 0.5 * (randlib::consts::ln2 + randlib::consts::lnpi + std::log(alpha) + std::log1p(-alpha) + logU);
    return std::exp(y);
}

template < typename RealType >
double StableDistribution<RealType>::pdfAtZero() const
{
    double y0 = 0.0;
    if (beta == 0.0)
        y0 = std::tgammal(alphaInv);
    else {
        y0 = std::lgammal(alphaInv) - omega;
        y0 = std::exp(y0) * std::cos(xi);
    }
    return y0 * randlib::consts::inv_pi / alpha;
}

template < typename RealType >
double StableDistribution<RealType>::pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x is 0, only this term is non-zero)
    double y0 = pdfAtZero();
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int n = 1; n <= k; ++n)
        {
            int n2 = n + n;
            double term = std::lgammal((n2 + 1) / alpha);
            term += n2 * logX;
            term -= RandMath::lfact(n2);
            term = std::exp(term);
            sum += (n & 1) ? -term : term;
        }
    }
    else {
        /// Asymmetric distribution
        double rhoPi_alpha = randlib::consts::half_pi + xiAdj;
        for (int n = 1; n <= k; ++n) {
            int np1 = n + 1;
            double term = std::lgammal(np1 * alphaInv);
            term += n * logX;
            term -= RandMath::lfact(n);
            term = std::exp(term - omega);
            term *= std::sin(np1 * rhoPi_alpha);
            sum += (n & 1) ? -term : term;
        }
    }
    return y0 + sum * randlib::consts::inv_pi / alpha;
}

template < typename RealType >
double StableDistribution<RealType>::pdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; ///< number of elements in the series
    double rhoPi = randlib::consts::half_pi + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha + 1.0;
        double term = std::lgammal(aux);
        term -= aux * logX;
        term -= RandMath::lfact(n);
        term = std::exp(term - omega);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return randlib::consts::inv_pi * sum;
}

template < typename RealType >
double StableDistribution<RealType>::pdfTaylorExpansionTailNearCauchy(double x) const
{
    double xSq = x * x;
    double y = 1.0 + xSq;
    double ySq = y * y;
    double z = RandMath::atan(x);
    double zSq = z * z;
    double logY = std::log1pl(xSq);
    double alpham1 = alpha - 1.0;
    double temp = 1.0 - randlib::consts::egamma - 0.5 * logY;
    /// first derivative
    double f_a = temp;
    f_a *= xSq - 1.0;
    f_a += 2 * x * z;
    f_a /= ySq;
    static constexpr long double M_PI_SQ_6 = 1.64493406684822643647l; /// π^2 / 6
    /// second derivative
    double f_aa1 = M_PI_SQ_6;
    f_aa1 += temp * temp;
    f_aa1 -= 1.0 + z * z;
    f_aa1 *= xSq * xSq - 6.0 * xSq + 1.0;
    double f_aa2 = 0.5 + temp;
    f_aa2 *= z;
    f_aa2 *= 8 * x * (xSq - 1.0);
    double f_aa3 = (1.0 - 3 * xSq) * temp;
    f_aa3 -= x * y * z;
    f_aa3 += f_aa3;
    double f_aa = f_aa1 + f_aa2 + f_aa3;
    f_aa /= std::pow(y, 3);
    /// Hashed values of special functions for x = 2, 3, 4
    /// Gamma(x)
    static constexpr int gammaTable[] = {1, 2, 6};
    /// Gamma'(x)
    static constexpr long double gammaDerTable[] = {1.0 - randlib::consts::egamma, 3.0 - 2.0 * randlib::consts::egamma, 11.0 - 6.0 * randlib::consts::egamma};
    /// Gamma''(x)
    static constexpr long double gammaSecDerTable[] = {0.82368066085287938958l, 2.49292999190269305794l, 11.1699273161019477314l};
    /// Digamma(x)
    static constexpr long double digammaTable[] = {1.0 - randlib::consts::egamma, 1.5 - randlib::consts::egamma, 11.0 / 6 - randlib::consts::egamma};
    /// Digamma'(x)
    static constexpr long double digammaDerTable[] = {M_PI_SQ_6 - 1.0, M_PI_SQ_6 - 1.25, M_PI_SQ_6 - 49.0 / 36};
    /// Digamma''(x)
    static constexpr long double digammaSecDerTable[] = {-0.40411380631918857080l, -0.15411380631918857080l, -0.08003973224511449673l};
    /// third derivative
    double gTable[] = {0, 0, 0};
    for (int i = 0; i < 3; ++i) {
        double g_11 = 0.25 * gammaTable[i] * logY * logY;
        g_11 -= gammaDerTable[i] * logY;
        g_11 += gammaSecDerTable[i];
        double aux = digammaTable[i] - 0.5 * logY;
        double g_12 = aux;
        double zip2 = z * (i + 2);
        double cosZNu = std::cos(zip2), zSinZNu = z * std::sin(zip2);
        g_12 *= cosZNu;
        g_12 -= zSinZNu;
        double g_1 = g_11 * g_12;
        double g_21 = -gammaTable[i] * logY + 2 * gammaDerTable[i];
        double g_22 = -zSinZNu * aux;
        g_22 -= zSq * cosZNu;
        g_22 += cosZNu * digammaDerTable[i];
        double g_2 = g_21 * g_22;
        double g_3 = -zSq * cosZNu * aux;
        g_3 -= 2 * zSinZNu * digammaDerTable[i];
        g_3 += zSq * zSinZNu;
        g_3 += cosZNu * digammaSecDerTable[i];
        g_3 *= gammaTable[i];
        double g = g_1 + g_2 + g_3;
        g *= std::pow(y, -0.5 * i - 1);
        gTable[i] = g;
    }
    double f_aaa = -gTable[0] + 3 * gTable[1] - gTable[2];
    /// sum all three derivatives
    double tail = f_a * alpham1;
    tail += 0.5 * f_aa * alpham1 * alpham1;
    tail += std::pow(alpham1, 3) * f_aaa / 6.0;
    tail /= randlib::consts::pi;
    return tail;
}

template < typename RealType >
double StableDistribution<RealType>::limitCaseForIntegrandAuxForGeneralExponent(double theta, double xiAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta < 0.5 * (randlib::consts::half_pi - xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
}

template < typename RealType >
double StableDistribution<RealType>::integrandAuxForGeneralExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi || theta <= -xiAdj)
        return limitCaseForIntegrandAuxForGeneralExponent(theta, xiAdj);
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    double y = std::log(std::cos(theta));
    y -= alpha * std::log(sinThetaAdj);
    y /= alpha - 1.0;
    y += std::log(std::cos(thetaAdj - theta));
    y += xAdj;
    return std::isfinite(y) ? y : limitCaseForIntegrandAuxForGeneralExponent(theta, xiAdj);
}

template < typename RealType >
double StableDistribution<RealType>::integrandFoGeneralExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= randlib::consts::half_pi)
        return 0.0;
    if (theta <= -xiAdj)
        return 0.0;
    double u = integrandAuxForGeneralExponent(theta, xAdj, xiAdj);
    return fastpdfExponentiation(u);
}

template < typename RealType >
double StableDistribution<RealType>::pdfForGeneralExponent(double x) const
{
    /// Standardize
    double xSt = (x - mu) / gamma;
    double absXSt = xSt;
    /// +- xi
    double xiAdj = xi;
    if (xSt > 0) {
        if (alpha < 1 && beta == -1)
            return 0.0;
    }
    else {
        if (alpha < 1 && beta == 1)
            return 0.0;
        absXSt = -xSt;
        xiAdj = -xi;
    }

    /// If α is too close to 1 and distribution is symmetric, then we approximate using Taylor series
    if (beta == 0.0 && std::fabs(alpha - 1.0) < 0.01)
        return pdfCauchy(x) + pdfTaylorExpansionTailNearCauchy(absXSt) / gamma;

    /// If x = 0, we know the analytic solution
    if (xSt == 0.0)
        return pdfAtZero() / gamma;

    double logAbsX = std::log(absXSt) - omega;

    /// If x is too close to 0, we do series expansion avoiding numerical problems
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1 && std::fabs(beta) == 1)
            return pdfShortTailExpansionForGeneralExponent(logAbsX);
        return pdfSeriesExpansionAtZero(logAbsX, xiAdj, seriesZeroParams.first) / gamma;
    }

    /// If x is large enough we use tail approximation
    if (logAbsX > pdftailBound && alpha <= ALMOST_TWO) {
        /// for α > 1 and |β| = 1 only the tail on the side of -β is short; the other one is heavy
        double skewnessTowardsX = (xSt > 0) ? beta : -beta;
        if (alpha > 1 && skewnessTowardsX == -1)
            return pdfShortTailExpansionForGeneralExponent(logAbsX);
        return pdfSeriesExpansionAtInf(logAbsX, xiAdj) / gamma;
    }

    double xAdj = alpha_alpham1 * logAbsX;

    /// Integrate around the peak of the integrand, where its exponent u vanishes
    auto u = [this, xAdj, xiAdj] (double theta) { return integrandAuxForGeneralExponent(theta, xAdj, xiAdj); };
    auto integrand = [this, xAdj, xiAdj] (double theta) { return integrandFoGeneralExponent(theta, xAdj, xiAdj); };
    double theta0 = rootOfMonotone(u, -xiAdj, randlib::consts::half_pi);
    double integral = integrateAroundRoot(integrand, u, -xiAdj, randlib::consts::half_pi, theta0);
    double res = pdfCoef * integral / absXSt;

    /// Finally we check if α is not too close to 2
    if (alpha <= ALMOST_TWO)
        return res;

    /// If α is near 2, we use tail aprroximation for large x
    /// and compare it with integral representation
    double alphap1 = alpha + 1.0;
    double tail = std::lgammal(alphap1);
    tail -= alphap1 * logAbsX;
    tail = std::exp(tail);
    tail *= (1.0 - 0.5 * alpha) / gamma;
    return std::max(tail, res);
}

template < typename RealType >
double StableDistribution<RealType>::f(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return pdfNormal(x);
    case CAUCHY:
        return pdfCauchy(x);
    case LEVY:
        return (beta > 0) ? pdfLevy(x) : pdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return pdfForUnityExponent(x);
    case GENERAL:
        return pdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::logf(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return logpdfNormal(x);
    case CAUCHY:
        return logpdfCauchy(x);
    case LEVY:
        return (beta > 0) ? logpdfLevy(x) : logpdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return std::log(pdfForUnityExponent(x));
    case GENERAL:
        return std::log(pdfForGeneralExponent(x));
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::cdfNormal(RealType x) const
{
    double y = mu - x;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

template < typename RealType >
double StableDistribution<RealType>::cdfNormalCompl(RealType x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

template < typename RealType >
double StableDistribution<RealType>::cdfCauchy(RealType x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    /// In the lower tail 1/2 + atan(x)/π cancels; atan(-1/x)/π does not
    return (x0 < -1.0) ? randlib::consts::inv_pi * std::atan(-1.0 / x0) : 0.5 + randlib::consts::inv_pi * std::atan(x0);
}

template < typename RealType >
double StableDistribution<RealType>::cdfCauchyCompl(RealType x) const
{
    double x0 = mu - x;
    x0 /= gamma;
    return (x0 < -1.0) ? randlib::consts::inv_pi * std::atan(-1.0 / x0) : 0.5 + randlib::consts::inv_pi * std::atan(x0);
}

template < typename RealType >
double StableDistribution<RealType>::cdfLevy(RealType x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

template < typename RealType >
double StableDistribution<RealType>::cdfLevyCompl(RealType x) const
{
    if (x <= mu)
        return 1.0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erf(y);
}

template < typename RealType >
double StableDistribution<RealType>::fastcdfExponentiation(double u)
{
    /// exp(-exp(u)), which underflows to 0 by itself for large u
    return std::exp(-std::exp(u));
}

template < typename RealType >
double StableDistribution<RealType>::fastccdfExponentiation(double u)
{
    return (u > 5.0) ? 1.0 : -std::expm1(-std::exp(u));
}

template < typename RealType >
double StableDistribution<RealType>::cdfAtZero(double xiAdj) const
{
    return 0.5 - randlib::consts::inv_pi * xiAdj;
}

template < typename RealType >
double StableDistribution<RealType>::cdfForUnityExponent(double x, bool complement) const
{
    /// y = (1/π) ∫ exp(-exp(u(θ))) dθ is F for β > 0 and S for β < 0; the complementary
    /// integrand 1 - exp(-exp(u)) gives the other one without cancellation. The integrand
    /// switches between 0 and 1 around the root θ0 of u, so integrate on both sides of it.
    double xSt = (x - mu) / gamma;
    double xAdj = -randlib::consts::half_pi * xSt / beta - logGammaPi_2;
    bool useComplement = complement == (beta > 0);
    auto integrand = [this, xAdj, useComplement] (double theta)
    {
        double u = integrandAuxForUnityExponent(theta, xAdj);
        return useComplement ? fastccdfExponentiation(u) : fastcdfExponentiation(u);
    };
    double theta0 = peakForUnityExponent(xAdj);
    auto u = [this, xAdj] (double theta) { return integrandAuxForUnityExponent(theta, xAdj); };
    return randlib::consts::inv_pi * integrateAroundRoot(integrand, u, -randlib::consts::half_pi, randlib::consts::half_pi, theta0);
}

template < typename RealType >
double StableDistribution<RealType>::cdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x = 0, only this term is non-zero)
    double y0 = cdfAtZero(xiAdj);
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int m = 0; m <= k; ++m) {
            int m2p1 = 2 * m + 1;
            double term = std::lgammal(m2p1 * alphaInv);
            term += m2p1 * logX;
            term -= RandMath::lfact(m2p1);
            term = std::exp(term);
            sum += (m & 1) ? -term : term;
        }
    }
    else {
        /// Asymmetric distribution
        double rhoPi_alpha = randlib::consts::half_pi + xiAdj;
        for (int n = 1; n <= k; ++n) {
            double term = std::lgammal(n * alphaInv);
            term += n * logX;
            term -= RandMath::lfact(n);
            term = std::exp(term);
            term *= std::sin(n * rhoPi_alpha);
            sum += (n & 1) ? term : -term;
        }
    }
    return y0 + sum * randlib::consts::inv_pi * alphaInv;
}

template < typename RealType >
double StableDistribution<RealType>::cdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; /// number of elements in the series
    double rhoPi = randlib::consts::half_pi + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha;
        double term = std::lgammal(aux);
        term -= aux * logX;
        term -= RandMath::lfact(n);
        term = std::exp(term);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return randlib::consts::inv_pi * sum;
}

template < typename RealType >
double StableDistribution<RealType>::cdfIntegralRepresentation(double logX, double xiAdj, bool complement) const
{
    /// (1/π) ∫_{-ξ}^{π/2} exp(-exp(u(θ))) dθ, or the same integral of 1 - exp(-exp(u(θ)))
    double xAdj = alpha_alpham1 * logX;
    auto u = [this, xAdj, xiAdj] (double theta) { return integrandAuxForGeneralExponent(theta, xAdj, xiAdj); };
    auto integrand = [&u, complement] (double theta) {
        return complement ? fastccdfExponentiation(u(theta)) : fastcdfExponentiation(u(theta));
    };
    double theta0 = rootOfMonotone(u, -xiAdj, randlib::consts::half_pi);
    return randlib::consts::inv_pi * integrateAroundRoot(integrand, u, -xiAdj, randlib::consts::half_pi, theta0);
}

template < typename RealType >
double StableDistribution<RealType>::cdfForGeneralExponent(double x) const
{
    double xSt = (x - mu) / gamma; /// Standardize
    if (xSt == 0) {
        /// totally skewed laws with α < 1 start (or end) exactly at μ
        if (alpha < 1 && std::fabs(beta) == 1)
            return (beta > 0) ? 0.0 : 1.0;
        return cdfAtZero(xi);
    }
    if (xSt > 0.0) {
        if (alpha > 1.0)
            return 1.0 - ccdfForGeneralExponent(x);
        /// α < 1: F is small near 0, so compute it directly rather than as 1 - S
        if (beta == -1.0)
            return 1.0;
        double logAbsX = std::log(xSt) - omega;
        if (logAbsX > cdftailBound)
            return 1.0 - cdfSeriesExpansionAtInf(logAbsX, xi);
        if (logAbsX < seriesZeroParams.second)
            return (beta == 1.0) ? cdfShortTailExpansion(logAbsX) : cdfSeriesExpansionAtZero(logAbsX, xi, seriesZeroParams.first);
        double atZero = (beta == 1.0) ? 0.0 : cdfAtZero(xi);
        return atZero + cdfIntegralRepresentation(logAbsX, xi, false);
    }
    /// For x < 0 we use relation F(x, ξ) = S(-x, -ξ)
    double logAbsX = std::log(-xSt) - omega;
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1.0 && beta == -1.0)
            return 1.0 - cdfShortTailExpansion(logAbsX);
        return 1.0 - cdfSeriesExpansionAtZero(logAbsX, -xi, seriesZeroParams.first);
    }
    if (logAbsX > cdftailBound)
        return cdfSeriesExpansionAtInf(logAbsX, -xi);
    if (alpha > 1.0)
        return cdfIntegralRepresentation(logAbsX, -xi, false);
    return (beta == 1.0) ? 0.0 : cdfIntegralRepresentation(logAbsX, -xi, true);
}

template < typename RealType >
double StableDistribution<RealType>::ccdfForGeneralExponent(double x) const
{
    double xSt = (x - mu) / gamma; /// Standardize
    if (xSt <= 0)
        return 1.0 - cdfForGeneralExponent(x);
    double logAbsX = std::log(xSt) - omega;
    /// If x is too close to 0, we do series expansion avoiding numerical problems
    if (logAbsX < seriesZeroParams.second) {
        if (alpha < 1.0 && beta == 1.0)
            return 1.0 - cdfShortTailExpansion(logAbsX);
        if (alpha < 1.0 && beta == -1.0)
            return cdfShortTailExpansion(logAbsX);
        return 1.0 - cdfSeriesExpansionAtZero(logAbsX, xi, seriesZeroParams.first);
    }
    /// If x is large enough we use tail approximation
    if (logAbsX > cdftailBound)
        return cdfSeriesExpansionAtInf(logAbsX, xi);
    if (alpha > 1.0)
        return cdfIntegralRepresentation(logAbsX, xi, false);
    return (beta == -1.0) ? 0.0 : cdfIntegralRepresentation(logAbsX, xi, true);
}

template < typename RealType >
double StableDistribution<RealType>::F(const RealType &x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return cdfNormal(x);
    case CAUCHY:
        return cdfCauchy(x);
    case LEVY:
        return (beta > 0) ? cdfLevy(x) : cdfLevyCompl(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x, false);
    case GENERAL:
        return cdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::S(const RealType & x) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return cdfNormalCompl(x);
    case CAUCHY:
        return cdfCauchyCompl(x);
    case LEVY:
        return (beta > 0) ? cdfLevyCompl(x) : cdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x, true);
    case GENERAL:
        return ccdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::variateForUnityExponent() const
{
    RealType U = randlib::consts::pi * UniformRand<RealType>::StandardVariate(this->localRandGenerator) - randlib::consts::half_pi;
    RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType pi_2pBetaU = randlib::consts::half_pi + beta * U;
    RealType Y = W * std::cos(U) / pi_2pBetaU;
    RealType X = std::log(Y);
    X += logGammaPi_2;
    X *= -beta;
    X += pi_2pBetaU * std::tan(U);
    X *= randlib::consts::two_inv_pi;
    return mu + gamma * X;
}

template < typename RealType >
double StableDistribution<RealType>::variateForGeneralExponent() const
{
    RealType U = randlib::consts::pi * UniformRand<RealType>::StandardVariate(this->localRandGenerator) - randlib::consts::half_pi;
    RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType alphaUpxi = alpha * (U + xi);
    RealType X = std::sin(alphaUpxi);
    RealType W_adj = W / std::cos(U - alphaUpxi);
    X *= W_adj;
    RealType R = omega - alphaInv * std::log(W_adj * std::cos(U));
    X *= std::exp(R);
    return mu + gamma * X;
}

template < typename RealType >
double StableDistribution<RealType>::variateForExponentEqualOneHalf() const
{
    RealType Z1 = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType Z2 = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType temp1 = (1.0 + beta) / Z1, temp2 = (1.0 - beta) / Z2;
    RealType var = temp1 - temp2;
    var *= temp1 + temp2;
    var *= 0.25;
    return mu + gamma * var;
}

template < typename RealType >
RealType StableDistribution<RealType>::Variate() const
{
    switch (distributionType) {
    case NORMAL:
        return mu + randlib::consts::sqrt2 * gamma * NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    case CAUCHY:
        return mu + gamma * CauchyRand<RealType>::StandardVariate(this->localRandGenerator);
    case LEVY:
        return mu + RandMath::sign(beta) * gamma * LevyRand<RealType>::StandardVariate(this->localRandGenerator);
    case UNITY_EXPONENT:
        return variateForUnityExponent();
    case GENERAL:
        return (alpha == 0.5) ? variateForExponentEqualOneHalf() : variateForGeneralExponent();
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
void StableDistribution<RealType>::Sample(std::vector<RealType> &outputData) const
{
    switch (distributionType) {
    case NORMAL: {
        double stdev = randlib::consts::sqrt2 * gamma;
        for (RealType &var : outputData)
            var = mu + stdev * NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    }
        break;
    case CAUCHY: {
        for (RealType &var : outputData)
            var = mu + gamma * CauchyRand<RealType>::StandardVariate(this->localRandGenerator);
    }
        break;
    case LEVY: {
        if (beta > 0) {
            for (RealType &var : outputData)
                var = mu + gamma * LevyRand<RealType>::StandardVariate(this->localRandGenerator);
        }
        else {
            for (RealType &var : outputData)
                var = mu - gamma * LevyRand<RealType>::StandardVariate(this->localRandGenerator);
        }
    }
        break;
    case UNITY_EXPONENT: {
        for (RealType &var : outputData)
            var = variateForUnityExponent();
    }
        break;
    case GENERAL: {
        if (alpha == 0.5) {
            for (RealType &var : outputData)
                var = variateForExponentEqualOneHalf();
        }
        else {
            for (RealType &var : outputData)
                var = variateForGeneralExponent();
        }
    }
        break;
    default:
        break;
    }
}

template < typename RealType >
long double StableDistribution<RealType>::Mean() const
{
    if (alpha > 1)
        return mu;
    if (beta == 1)
        return INFINITY;
    return (beta == -1) ? -INFINITY : NAN;
}

template < typename RealType >
long double StableDistribution<RealType>::Variance() const
{
    return (distributionType == NORMAL) ? 2 * gamma * gamma : INFINITY;
}

template < typename RealType >
RealType StableDistribution<RealType>::Mode() const
{
    checkDensitySupport();
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionType == NORMAL)
        return mu;
    if (distributionType == LEVY)
        return mu + beta * gamma / 3.0;
    return ContinuousDistribution<RealType>::Mode();
}

template < typename RealType >
RealType StableDistribution<RealType>::Median() const
{
    checkDensitySupport();
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionType == NORMAL)
        return mu;
    return ContinuousDistribution<RealType>::Median();
}

template < typename RealType >
long double StableDistribution<RealType>::Skewness() const
{
    return (distributionType == NORMAL) ? 0 : NAN;
}

template < typename RealType >
long double StableDistribution<RealType>::ExcessKurtosis() const
{
    return (distributionType == NORMAL) ? 0 : NAN;
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileNormal(double p) const
{
    return mu - 2 * gamma * RandMath::erfcinv(2 * p);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileNormal1m(double p) const
{
    return mu + 2 * gamma * RandMath::erfcinv(2 * p);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileCauchy(double p) const
{
    return mu - gamma / std::tan(randlib::consts::pi * p);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileCauchy1m(double p) const
{
    return mu + gamma / std::tan(randlib::consts::pi * p);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileLevy(double p) const
{
    double y = RandMath::erfcinv(p);
    return mu + 0.5 * gamma / (y * y);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileLevy1m(double p) const
{
    double y = RandMath::erfinv(p);
    return mu + 0.5 * gamma / (y * y);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileImpl(double p) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return quantileNormal(p);
    case CAUCHY:
        return quantileCauchy(p);
    case LEVY:
        return (beta > 0) ? quantileLevy(p) : 2 * mu - quantileLevy1m(p);
    default:
        return ContinuousDistribution<RealType>::quantileImpl(p);
    }
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileImpl1m(double p) const
{
    checkDensitySupport();
    switch (distributionType) {
    case NORMAL:
        return quantileNormal1m(p);
    case CAUCHY:
        return quantileCauchy1m(p);
    case LEVY:
        return (beta > 0) ? quantileLevy1m(p) : 2 * mu - quantileLevy(p);
    default:
        return ContinuousDistribution<RealType>::quantileImpl1m(p);
    }
}

template < typename RealType >
std::complex<double> StableDistribution<RealType>::cfNormal(double t) const
{
    double gammaT = gamma * t;
    std::complex<double> y(-gammaT * gammaT, mu * t);
    return std::exp(y);
}

template < typename RealType >
std::complex<double> StableDistribution<RealType>::cfCauchy(double t) const
{
    std::complex<double> y(-gamma * t, mu * t);
    return std::exp(y);
}

template < typename RealType >
std::complex<double> StableDistribution<RealType>::cfLevy(double t) const
{
    std::complex<double> y(0.0, -2 * gamma * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

template < typename RealType >
std::complex<double> StableDistribution<RealType>::CFImpl(double t) const
{
    double x = 0;
    switch (distributionType) {
    case NORMAL:
        return cfNormal(t);
    case CAUCHY:
        return cfCauchy(t);
    case LEVY: {
        std::complex<double> phi = cfLevy(t);
        return (beta > 0) ? phi : std::conj(phi);
    }
    case UNITY_EXPONENT:
        x = beta * randlib::consts::two_inv_pi * std::log(t);
        break;
    default:
        x = zeta;
    }
    double re = std::pow(gamma * t, alpha);
    std::complex<double> psi = std::complex<double>(re, re * x - mu * t);
    return std::exp(-psi);
}

template class StableDistribution<float>;
template class StableDistribution<double>;
template class StableDistribution<long double>;

template < typename RealType >
String StableRand<RealType>::Name() const
{
    return "Stable("
            + this->toStringWithPrecision(this->GetExponent()) + ", "
            + this->toStringWithPrecision(this->GetSkewness()) + ", "
            + this->toStringWithPrecision(this->GetScale()) + ", "
            + this->toStringWithPrecision(this->GetLocation()) + ")";
}

template < typename RealType >
void StableRand<RealType>::SetExponent(double exponent)
{
    this->SetParameters(exponent, this->GetSkewness(), this->GetScale(), this->GetLocation());
}

template < typename RealType >
void StableRand<RealType>::SetSkewness(double skewness)
{
    this->SetParameters(this->GetExponent(), skewness, this->GetScale(), this->GetLocation());
}

template class StableRand<float>;
template class StableRand<double>;
template class StableRand<long double>;

template < typename RealType >
String HoltsmarkRand<RealType>::Name() const
{
    return "Holtsmark("
            + this->toStringWithPrecision(this->GetScale()) + ", "
            + this->toStringWithPrecision(this->GetLocation()) + ")";
}

template class HoltsmarkRand<float>;
template class HoltsmarkRand<double>;
template class HoltsmarkRand<long double>;

template < typename RealType >
String LandauRand<RealType>::Name() const
{
    return "Landau("
            + this->toStringWithPrecision(this->GetScale()) + ", "
            + this->toStringWithPrecision(this->GetLocation()) + ")";
}

template class LandauRand<float>;
template class LandauRand<double>;
template class LandauRand<long double>;
