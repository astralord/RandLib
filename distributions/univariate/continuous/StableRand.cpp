#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include <functional>

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
    if (exponent < 0.1 || exponent > 2.0)
        throw std::invalid_argument("Stable distribution: exponent should be in the interval [0.1, 2], but it's equal to "
                                    + std::to_string(exponent));
    if (std::fabs(skewness) > 1.0)
        throw std::invalid_argument("Stable distribution: skewness should be in the interval [-1, 1], but it's equal to "
                                    + std::to_string(skewness));
    if (scale <= 0.0)
        throw std::invalid_argument("Stable distribution: scale should be positive, but it's equal to "
                                    + std::to_string(scale));

    /// the following errors should be removed in the future
    if (exponent != 1.0 && std::fabs(exponent - 1.0) < 0.01 && skewness != 0.0)
        throw std::invalid_argument("Stable distribution: exponent close to 1 with non-zero skewness is not yet supported");
    if (exponent == 1.0 && skewness != 0.0 && std::fabs(skewness) < 0.01)
        throw std::invalid_argument("Stable distribution: skewness close to 0 with exponent equal to 1 is not yet supported");
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForNormal()
{   
    distributionType = NORMAL;
    pdfCoef = M_LN2 + logGamma + 0.5 * M_LNPI;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForCauchy()
{   
    distributionType = CAUCHY;
    pdfCoef = -logGamma - M_LNPI;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForLevy()
{   
    distributionType = LEVY;
    pdfCoef = logGamma - M_LN2 - M_LNPI;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForUnityExponent()
{   
    distributionType = UNITY_EXPONENT;
    pdfCoef = 0.5 / (gamma * std::fabs(beta));
    pdftailBound = 0; // not in the use for now
    logGammaPi_2 = logGamma + M_LNPI - M_LN2;
}

template < typename RealType >
void StableDistribution<RealType>::setParametersForGeneralExponent()
{   
    distributionType = GENERAL;
    if (beta != 0.0) {
        zeta = -beta * std::tan(M_PI_2 * alpha);
        omega = 0.5 * alphaInv * std::log1pl(zeta * zeta);
        xi = alphaInv * RandMath::atan(-zeta);
    }
    else {
        zeta = omega = xi = 0.0;
    }
    pdfCoef = M_1_PI * std::fabs(alpha_alpham1) / gamma;
    pdftailBound = 3.0 / (1.0 + alpha) * M_LN10;
    cdftailBound = 3.0 / alpha * M_LN10;
    /// define boundaries of region near 0, where we use series expansion
    if (alpha <= ALMOST_TWO) {
        seriesZeroParams.first = std::round(std::min(alpha * alpha * 40 + 1, 10.0));
        seriesZeroParams.second = -(alphaInv * 1.5 + 0.5) * M_LN10; /// corresponds to interval [10^(-15.5), ~0.056]
    }
    else {
        seriesZeroParams.first = 85;
        seriesZeroParams.second = M_LN2 + M_LN3;///< corresponds to 6
    }
}

template < typename RealType >
void StableDistribution<RealType>::SetParameters(double exponent, double skewness, double scale, double location)
{
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
    return M_1_PI / y;
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
    if (u > 5 || u < -50)
        return 0.0;
    return (u < -25) ? std::exp(u) : std::exp(u - std::exp(u));
}

template < typename RealType >
double StableDistribution<RealType>::pdfShortTailExpansionForUnityExponent(double x) const
{
    if (x > 10)
        return 0.0;
    double xm1 = x - 1.0;
    double y = 0.5 * xm1 - std::exp(xm1);
    y -= 0.5 * (M_LN2 + M_LNPI);
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
    if (std::fabs(theta) >= M_PI_2)
        return limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
    if (theta == 0.0)
        return xAdj + M_LNPI - M_LN2;
    double thetaAdj = (M_PI_2 + beta * theta) / std::cos(theta);
    double u = std::log(thetaAdj);
    u += thetaAdj * std::sin(theta) / beta;
    return std::isfinite(u) ? u + xAdj : limitCaseForIntegrandAuxForUnityExponent(theta, xAdj);
}

template < typename RealType >
double StableDistribution<RealType>::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    return fastpdfExponentiation(u);
}

template < typename RealType >
double StableDistribution<RealType>::pdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / gamma;
    double xAdj = -M_PI_2 * xSt / beta - logGammaPi_2;

    /// We squeeze boudaries for too peaked integrands
    double boundary = RandMath::atan(M_2_PI * beta * (5.0 - xAdj));
    double upperBoundary = (beta > 0.0) ? boundary : M_PI_2;
    double lowerBoundary = (beta < 0.0) ? boundary : -M_PI_2;

    /// Find peak of the integrand
    double theta0 = 0;
    std::function<double (double)> funPtr = std::bind(&StableDistribution<RealType>::integrandAuxForUnityExponent, this, std::placeholders::_1, xAdj);
    RandMath::findRootNewtonFirstOrder(funPtr, lowerBoundary, upperBoundary, theta0);

    /// Sanity check
    /// if we failed while looking for the peak position
    /// we set it in the middle between boundaries
    if (theta0 >= upperBoundary || theta0 <= lowerBoundary)
        theta0 = 0.5 * (upperBoundary + lowerBoundary);

    std::function<double (double)> integrandPtr = std::bind(&StableDistribution<RealType>::integrandForUnityExponent, this, std::placeholders::_1, xAdj);

    /// If theta0 is too close to +/-π/2 then we can still underestimate the integral
    int maxRecursionDepth = 11;
    double closeness = M_PI_2 - std::fabs(theta0);
    if (closeness < 0.1)
        maxRecursionDepth = 20;
    else if (closeness < 0.2)
        maxRecursionDepth = 15;

    double int1 = RandMath::integral(integrandPtr, lowerBoundary, theta0, 1e-11, maxRecursionDepth);
    double int2 = RandMath::integral(integrandPtr, theta0, upperBoundary, 1e-11, maxRecursionDepth);
    return pdfCoef * (int1 + int2);
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
    z -= 0.5 * (M_LN2 + M_LNPI + logAlpha + log1mAlpha);
    return std::exp(y + z);
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
    return y0 * M_1_PI / alpha;
}

template < typename RealType >
double StableDistribution<RealType>::pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x = 0, only this term is non-zero)
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
        double rhoPi_alpha = M_PI_2 + xiAdj;
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
    return y0 + sum * M_1_PI / alpha;
}

template < typename RealType >
double StableDistribution<RealType>::pdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; ///< number of elements in the series
    double rhoPi = M_PI_2 + xiAdj;
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
    return M_1_PI * sum;
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
    double temp = 1.0 - M_EULER - 0.5 * logY;
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
    static constexpr long double gammaDerTable[] = {1.0 - M_EULER, 3.0 - 2.0 * M_EULER, 11.0 - 6.0 * M_EULER};
    /// Gamma''(x)
    static constexpr long double gammaSecDerTable[] = {0.82368066085287938958l, 2.49292999190269305794l, 11.1699273161019477314l};
    /// Digamma(x)
    static constexpr long double digammaTable[] = {1.0 - M_EULER, 1.5 - M_EULER, 11.0 / 6 - M_EULER};
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
    tail /= M_PI;
    return tail;
}

template < typename RealType >
double StableDistribution<RealType>::limitCaseForIntegrandAuxForGeneralExponent(double theta, double xiAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta < 0.5 * (M_PI_2 - xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
}

template < typename RealType >
double StableDistribution<RealType>::integrandAuxForGeneralExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= M_PI_2 || theta <= -xiAdj)
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
    if (std::fabs(theta) >= M_PI_2)
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
        if (alpha > 1 && std::fabs(beta) == 1)
            return pdfShortTailExpansionForGeneralExponent(logAbsX);
        return pdfSeriesExpansionAtInf(logAbsX, xiAdj) / gamma;
    }

    double xAdj = alpha_alpham1 * logAbsX;

    /// Search for the peak of the integrand
    double theta0;
    std::function<double (double)> funPtr = std::bind(&StableDistribution<RealType>::integrandAuxForGeneralExponent, this, std::placeholders::_1, xAdj, xiAdj);
    RandMath::findRootNewtonFirstOrder(funPtr, -xiAdj, M_PI_2, theta0);

    /// If theta0 is too close to π/2 or -xiAdj then we can still underestimate the integral
    int maxRecursionDepth = 11;
    double closeness = std::min(M_PI_2 - theta0, theta0 + xiAdj);
    if (closeness < 0.1)
        maxRecursionDepth = 20;
    else if (closeness < 0.2)
        maxRecursionDepth = 15;

    /// Calculate sum of two integrals
    std::function<double (double)> integrandPtr = std::bind(&StableDistribution<RealType>::integrandFoGeneralExponent, this, std::placeholders::_1, xAdj, xiAdj);
    double int1 = RandMath::integral(integrandPtr, -xiAdj, theta0, 1e-11, maxRecursionDepth);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2, 1e-11, maxRecursionDepth);
    double res = pdfCoef * (int1 + int2) / absXSt;

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
    return 0.5 + M_1_PI * RandMath::atan(x0);
}

template < typename RealType >
double StableDistribution<RealType>::cdfCauchyCompl(RealType x) const
{
    double x0 = mu - x;
    x0 /= gamma;
    return 0.5 + M_1_PI * RandMath::atan(x0);
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
    if (u > 5.0)
        return 0.0;
    else if (u < -50.0)
        return 1.0;
    double y = std::exp(u);
    return std::exp(-y);
}

template < typename RealType >
double StableDistribution<RealType>::cdfAtZero(double xiAdj) const
{
    return 0.5 - M_1_PI * xiAdj;
}

template < typename RealType >
double StableDistribution<RealType>::cdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / gamma;
    double xAdj = -M_PI_2 * xSt / beta - logGammaPi_2;
    double y = M_1_PI * RandMath::integral([this, xAdj] (double theta)
    {
        double u = integrandAuxForUnityExponent(theta, xAdj);
        return fastcdfExponentiation(u);
    },
    -M_PI_2, M_PI_2);
    return (beta > 0) ? y : 1.0 - y;
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
        double rhoPi_alpha = M_PI_2 + xiAdj;
        for (int n = 1; n <= k; ++n) {
            double term = std::lgammal(n * alphaInv);
            term += n * logX;
            term -= RandMath::lfact(n);
            term = std::exp(term);
            term *= std::sin(n * rhoPi_alpha);
            sum += (n & 1) ? term : -term;
        }
    }
    return y0 + sum * M_1_PI * alphaInv;
}

template < typename RealType >
double StableDistribution<RealType>::cdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; /// number of elements in the series
    double rhoPi = M_PI_2 + xiAdj;
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
    return M_1_PI * sum;
}

template < typename RealType >
double StableDistribution<RealType>::cdfIntegralRepresentation(double logX, double xiAdj) const
{
    double xAdj = alpha_alpham1 * logX;
    return M_1_PI * RandMath::integral([this, xAdj, xiAdj] (double theta)
    {
        double u = integrandAuxForGeneralExponent(theta, xAdj, xiAdj);
        return fastcdfExponentiation(u);
    },
    -xiAdj, M_PI_2);
}

template < typename RealType >
double StableDistribution<RealType>::cdfForGeneralExponent(double x) const
{
    double xSt = (x - mu) / gamma; /// Standardize
    if (xSt == 0)
        return cdfAtZero(xi);
    if (xSt > 0.0) {
        double logAbsX = std::log(xSt) - omega;
        /// If x is too close to 0, we do series expansion avoiding numerical problems
        if (logAbsX < seriesZeroParams.second)
            return cdfSeriesExpansionAtZero(logAbsX, xi, seriesZeroParams.first);
        /// If x is large enough we use tail approximation
        if (logAbsX > cdftailBound)
            return 1.0 - cdfSeriesExpansionAtInf(logAbsX, xi);
        if (alpha > 1.0)
            return 1.0 - cdfIntegralRepresentation(logAbsX, xi);
        return (beta == -1.0) ? 1.0 : cdfAtZero(xi) + cdfIntegralRepresentation(logAbsX, xi);
    }
    /// For x < 0 we use relation F(-x, xi) + F(x, -xi) = 1
    double logAbsX = std::log(-xSt) - omega;
    if (logAbsX < seriesZeroParams.second)
        return 1.0 - cdfSeriesExpansionAtZero(logAbsX, -xi, seriesZeroParams.first);
    if (logAbsX > cdftailBound)
        return cdfSeriesExpansionAtInf(logAbsX, -xi);
    if (alpha > 1.0)
        return cdfIntegralRepresentation(logAbsX, -xi);
    return (beta == 1.0) ? 0.0 : cdfAtZero(xi) - cdfIntegralRepresentation(logAbsX, -xi);
}

template < typename RealType >
double StableDistribution<RealType>::F(const RealType &x) const
{
    switch (distributionType) {
    case NORMAL:
        return cdfNormal(x);
    case CAUCHY:
        return cdfCauchy(x);
    case LEVY:
        return (beta > 0) ? cdfLevy(x) : cdfLevyCompl(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x);
    case GENERAL:
        return cdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::S(const RealType & x) const
{
    switch (distributionType) {
    case NORMAL:
        return cdfNormalCompl(x);
    case CAUCHY:
        return cdfCauchyCompl(x);
    case LEVY:
        return (beta > 0) ? cdfLevyCompl(x) : cdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return 1.0 - cdfForUnityExponent(x);
    case GENERAL:
        return 1.0 - cdfForGeneralExponent(x);
    default:
        throw std::runtime_error("Stable distribution: invalid distribution type");
    }
}

template < typename RealType >
double StableDistribution<RealType>::variateForUnityExponent() const
{
    RealType U = M_PI * UniformRand<RealType>::StandardVariate(this->localRandGenerator) - M_PI_2;
    RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType pi_2pBetaU = M_PI_2 + beta * U;
    RealType Y = W * std::cos(U) / pi_2pBetaU;
    RealType X = std::log(Y);
    X += logGammaPi_2;
    X *= -beta;
    X += pi_2pBetaU * std::tan(U);
    X *= M_2_PI;
    return mu + gamma * X;
}

template < typename RealType >
double StableDistribution<RealType>::variateForGeneralExponent() const
{
    RealType U = M_PI * UniformRand<RealType>::StandardVariate(this->localRandGenerator) - M_PI_2;
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
        return mu + M_SQRT2 * gamma * NormalRand<RealType>::StandardVariate(this->localRandGenerator);
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
        double stdev = M_SQRT2 * gamma;
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
    return mu - gamma / std::tan(M_PI * p);
}

template < typename RealType >
RealType StableDistribution<RealType>::quantileCauchy1m(double p) const
{
    return mu + gamma / std::tan(M_PI * p);
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
        x = beta * M_2_PI * std::log(t);
        break;
    default:
        x = -zeta;
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
