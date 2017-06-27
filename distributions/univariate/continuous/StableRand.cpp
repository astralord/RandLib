#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include <functional>

StableDistribution::StableDistribution(double exponent, double skewness, double scale, double location)
{
    SetParameters(exponent, skewness, scale, location);
}

SUPPORT_TYPE StableDistribution::SupportType() const
{
    if (alpha < 1) {
        if (beta == 1)
            return RIGHTSEMIFINITE_T;
        if (beta == -1)
            return LEFTSEMIFINITE_T;
    }
    return INFINITE_T;
}

double StableDistribution::MinValue() const
{
    return (alpha < 1 && beta == 1) ? mu : -INFINITY;
}

double StableDistribution::MaxValue() const
{
    return (alpha < 1 && beta == -1) ? mu : INFINITY;
}

void StableDistribution::SetParameters(double exponent, double skewness, double scale, double location)
{
    alpha = std::min(std::max(exponent, 0.1), 2.0);
    alphaInv = 1.0 / alpha;
    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);
    mu = location;
    gamma = scale > 0 ? scale : M_SQRT2;
    logGamma = std::log(gamma);

    /// Set id of distribution
    if (alpha == 2.0)
        distributionId = NORMAL;
    else if (alpha == 1.0)
        distributionId = (beta == 0.0) ? CAUCHY : UNITY_EXPONENT;
    else if (alpha == 0.5 && std::fabs(beta) == 1.0)
        distributionId = LEVY;
    else
        distributionId = COMMON;

    alpha_alpham1 = alpha / (alpha - 1.0);

    if (distributionId == NORMAL)
        pdfCoef = M_LN2 + logGamma + 0.5 * M_LNPI;
    else if (distributionId == LEVY)
        pdfCoef = logGamma - M_LN2 - M_LNPI;
    else if (distributionId == CAUCHY)
        pdfCoef = -logGamma - M_LNPI;
    else if (distributionId == UNITY_EXPONENT) {
        pdfCoef = 0.5 / (gamma * std::fabs(beta));
        pdftailBound = std::sqrt(2e4 / M_PI * M_E);
        logGammaPi_2 = logGamma + M_LNPI - M_LN2;
    }
    else if (distributionId == COMMON) {
        if (beta != 0.0) {
            zeta = -beta * std::tan(M_PI_2 * alpha);
            omega = 0.5 * alphaInv * std::log1p(zeta * zeta);
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
            /// corresponds to boundaries from 10^(-15.5) to ~ 0.056
            seriesZeroParams.second = -(alphaInv * 1.5 + 0.5) * M_LN10;
        }
        else {
            seriesZeroParams.first = 85;
            /// corresponds to 6
            seriesZeroParams.second = M_LN2 + M_LN3;
        }
    }
}

void StableDistribution::SetLocation(double location)
{
    mu = location;
}

void StableDistribution::SetScale(double scale)
{
    gamma = scale > 0 ? scale : M_SQRT2;
    logGamma = std::log(gamma);
    if (distributionId == NORMAL)
        pdfCoef = M_LN2 + logGamma + 0.5 * M_LNPI;
    else if (distributionId == LEVY)
        pdfCoef = logGamma - M_LN2 - M_LNPI;
    else if (distributionId == CAUCHY)
        pdfCoef = -logGamma - M_LNPI;
    else if (distributionId == UNITY_EXPONENT) {
        pdfCoef = 0.5 / (gamma * std::fabs(beta));
        logGammaPi_2 = logGamma + M_LNPI - M_LN2;
    }
    else if (distributionId == COMMON)
        pdfCoef = M_1_PI * std::fabs(alpha_alpham1) / gamma;
}

double StableDistribution::pdfNormal(double x) const
{
    return std::exp(logpdfNormal(x));
}

double StableDistribution::logpdfNormal(double x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    y *= y;
    y += pdfCoef;
    return -y;
}

double StableDistribution::pdfCauchy(double x) const
{
    double y = x - mu;
    y *= y;
    y /= gamma;
    y += gamma;
    return M_1_PI / y;
}

double StableDistribution::logpdfCauchy(double x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    double xSq = x0 * x0;
    return pdfCoef - std::log1p(xSq);
}

double StableDistribution::pdfLevy(double x) const
{
    return (x <= mu) ? 0.0 : std::exp(logpdfLevy(x));
}

double StableDistribution::logpdfLevy(double x) const
{
    double x0 = x - mu;
    if (x0 <= 0.0)
        return -INFINITY;
    double y = gamma / x0;
    y += 3 * std::log(x0);
    y -= pdfCoef;
    return -0.5 * y;
}

double StableDistribution::fastpdfExponentiation(double u)
{
    if (u > 5 || u < -50)
        return 0.0;
    return (u < -25) ? std::exp(u) : std::exp(u - std::exp(u));
}

double StableDistribution::limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta > 0.0) {
        if (beta > 0.0)
            return BIG_NUMBER;
        return (beta == -1) ? xAdj - 1.0 : -BIG_NUMBER;
    }
    if (beta < 0.0)
        return BIG_NUMBER;
    return (beta == 1) ? xAdj - 1.0 : -BIG_NUMBER;
}

double StableDistribution::integrandAuxForUnityExponent(double theta, double xAdj) const
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

double StableDistribution::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    return fastpdfExponentiation(u);
}

double StableDistribution::pdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / gamma;
    double xAdj = -M_PI_2 * xSt / beta - logGammaPi_2;

    // TODO: this limit is not available anymore, re-check
    if (std::fabs(xSt) > pdftailBound) {
        double y = 1.0 / (M_PI * xSt * xSt);
        y *= (xSt > 0) ? (1 + beta) : (1 - beta);
        return y;
    }

    /// We squeeze boudaries for too peaked integrands
    double boundary = RandMath::atan(M_2_PI * beta * (5.0 - xAdj));
    double upperBoundary = (beta > 0.0) ? boundary : M_PI_2;
    double lowerBoundary = (beta < 0.0) ? boundary : -M_PI_2;

    /// Find peak of the integrand
    double theta0 = 0;
    std::function<double (double)> funPtr = std::bind(&StableDistribution::integrandAuxForUnityExponent, this, std::placeholders::_1, xAdj);
    RandMath::findRoot(funPtr, lowerBoundary, upperBoundary, theta0);

    /// Sanity check
    /// if we failed while looking for the peak position
    /// we set it in the middle between boundaries
    if (theta0 >= upperBoundary || theta0 <= lowerBoundary)
        theta0 = 0.5 * (upperBoundary + lowerBoundary);

    std::function<double (double)> integrandPtr = std::bind(&StableDistribution::integrandForUnityExponent, this, std::placeholders::_1, xAdj);

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

double StableDistribution::pdfAtZero() const
{
    double y0 = 0.0;
    if (beta == 0.0)
        y0 = std::tgamma(alphaInv);
    else {
        y0 = std::lgamma(alphaInv) - omega;
        y0 = std::exp(y0) * std::cos(xi);
    }
    return y0 * M_1_PI / alpha;
}

double StableDistribution::pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
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
            double term = std::lgamma((n2 + 1) / alpha);
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
            double term = std::lgamma(np1 * alphaInv);
            term += n * logX;
            term -= RandMath::lfact(n);
            term = std::exp(term - omega);
            term *= std::sin(np1 * rhoPi_alpha);
            sum += (n & 1) ? -term : term;
        }
    }
    return y0 + sum * M_1_PI / alpha;
}

double StableDistribution::pdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; /// number of elements in the series
    double rhoPi = M_PI_2 + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha + 1.0;
        double term = std::lgamma(aux);
        term -= aux * logX;
        term -= RandMath::lfact(n);
        term = std::exp(term - omega);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return M_1_PI * sum;
}

double StableDistribution::pdfTaylorExpansionTailNearCauchy(double x) const
{
    // TODO: recheck derivatives (there are some typos in the paper)
    double xSq = x * x;
    double y = 1.0 + xSq;
    double ySq = y * y;
    double z = RandMath::atan(x);
    double zSq = z * z;
    double logY = std::log1p(xSq);
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
    /// summarize all three derivatives
    double tail = f_a * alpham1;
    tail += 0.5 * f_aa * alpham1 * alpham1;
    tail += std::pow(alpham1, 3) * f_aaa / 6.0;
    tail /= M_PI;
    return tail;
}

double StableDistribution::limitCaseForIntegrandAuxForCommonExponent(double theta, double xiAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta < 0.5 * (M_PI_2 - xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
}

double StableDistribution::integrandAuxForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= M_PI_2 || theta <= -xiAdj)
        return limitCaseForIntegrandAuxForCommonExponent(theta, xiAdj);
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    double y = std::log(std::cos(theta));
    y -= alpha * std::log(sinThetaAdj);
    y /= alpha - 1.0;
    y += std::log(std::cos(thetaAdj - theta));
    y += xAdj;
    return std::isfinite(y) ? y : limitCaseForIntegrandAuxForCommonExponent(theta, xiAdj);
}

double StableDistribution::integrandForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    if (theta <= -xiAdj)
        return 0.0;
    double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
    return fastpdfExponentiation(u);
}

double StableDistribution::pdfForCommonExponent(double x) const
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
    if (logAbsX < seriesZeroParams.second)
        return pdfSeriesExpansionAtZero(logAbsX, xiAdj, seriesZeroParams.first) / gamma;

    /// If x is large enough we use tail approximation
    if (logAbsX > pdftailBound && alpha <= ALMOST_TWO)
        return pdfSeriesExpansionAtInf(logAbsX, xiAdj) / gamma;

    double xAdj = alpha_alpham1 * logAbsX;

    /// Search for the peak of the integrand
    double theta0;
    std::function<double (double)> funPtr = std::bind(&StableDistribution::integrandAuxForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
    RandMath::findRoot(funPtr, -xiAdj, M_PI_2, theta0);

    /// If theta0 is too close to π/2 or -xiAdj then we can still underestimate the integral
    int maxRecursionDepth = 11;
    double closeness = M_PI_2 - theta0;
    closeness = std::min(closeness, theta0 + xiAdj);
    if (closeness < 0.1)
        maxRecursionDepth = 20;
    else if (closeness < 0.2)
        maxRecursionDepth = 15;

    /// Calculate sum of two integrals
    std::function<double (double)> integrandPtr = std::bind(&StableDistribution::integrandForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
    double int1 = RandMath::integral(integrandPtr, -xiAdj, theta0, 1e-11, maxRecursionDepth);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2, 1e-11, maxRecursionDepth);
    double res = pdfCoef * (int1 + int2) / absXSt;

    /// Finally we check if α is not too close to 2
    if (alpha <= ALMOST_TWO)
        return res;

    /// If α is near 2, we use tail aprroximation for large x
    /// and compare it with integral representation
    double alphap1 = alpha + 1.0;
    double tail = std::lgamma(alphap1);
    tail -= alphap1 * logAbsX;
    tail = std::exp(tail);
    tail *= (1.0 - 0.5 * alpha) / gamma;
    return std::max(tail, res);
}

double StableDistribution::f(const double & x) const
{
    switch (distributionId) {
    case NORMAL:
        return pdfNormal(x);
    case CAUCHY:
        return pdfCauchy(x);
    case LEVY:
        return (beta > 0) ? pdfLevy(x) : pdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return pdfForUnityExponent(x);
    case COMMON:
        return pdfForCommonExponent(x);
    default:
        return NAN; /// unexpected return
    }
}

double StableDistribution::logf(const double & x) const
{
    switch (distributionId) {
    case NORMAL:
        return logpdfNormal(x);
    case CAUCHY:
        return logpdfCauchy(x);
    case LEVY:
        return (beta > 0) ? logpdfLevy(x) : logpdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return std::log(pdfForUnityExponent(x));
    case COMMON:
        return std::log(pdfForCommonExponent(x));
    default:
        return NAN; /// unexpected return
    }
}

double StableDistribution::cdfNormal(double x) const
{
    double y = mu - x;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

double StableDistribution::cdfNormalCompl(double x) const
{
    double y = x - mu;
    y *= 0.5 / gamma;
    return 0.5 * std::erfc(y);
}

double StableDistribution::cdfCauchy(double x) const
{
    double x0 = x - mu;
    x0 /= gamma;
    return 0.5 + M_1_PI * RandMath::atan(x0);
}

double StableDistribution::cdfCauchyCompl(double x) const
{
    double x0 = mu - x;
    x0 /= gamma;
    return 0.5 + M_1_PI * RandMath::atan(x0);
}

double StableDistribution::cdfLevy(double x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

double StableDistribution::cdfLevyCompl(double x) const
{
    if (x <= mu)
        return 1.0;
    double y = x - mu;
    y += y;
    y = gamma / y;
    y = std::sqrt(y);
    return std::erf(y);
}

double StableDistribution::fastcdfExponentiation(double u)
{
    if (u > 5.0)
        return 0.0;
    else if (u < -50.0)
        return 1.0;
    double y = std::exp(u);
    return std::exp(-y);
}

double StableDistribution::cdfAtZero(double xiAdj) const
{
    return 0.5 - M_1_PI * xiAdj;
}

double StableDistribution::cdfForUnityExponent(double x) const
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

double StableDistribution::cdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x = 0, only this term is non-zero)
    double y0 = cdfAtZero(xiAdj);
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int m = 0; m <= k; ++m) {
            int m2p1 = 2 * m + 1;
            double term = std::lgamma(m2p1 * alphaInv);
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
            double term = std::lgamma(n * alphaInv);
            term += n * logX;
            term -= RandMath::lfact(n);
            term = std::exp(term);
            term *= std::sin(n * rhoPi_alpha);
            sum += (n & 1) ? term : -term;
        }
    }
    return y0 + sum * M_1_PI * alphaInv;
}

double StableDistribution::cdfSeriesExpansionAtInf(double logX, double xiAdj) const
{
    static constexpr int k = 10; /// number of elements in the series
    double rhoPi = M_PI_2 + xiAdj;
    rhoPi *= alpha;
    double sum = 0.0;
    for (int n = 1; n <= k; ++n) {
        double aux = n * alpha;
        double term = std::lgamma(aux);
        term -= aux * logX;
        term -= RandMath::lfact(n);
        term = std::exp(term);
        term *= std::sin(rhoPi * n);
        sum += (n & 1) ? term : -term;
    }
    return M_1_PI * sum;
}

double StableDistribution::cdfIntegralRepresentation(double logX, double xiAdj) const
{
    double xAdj = alpha_alpham1 * logX;
    return M_1_PI * RandMath::integral([this, xAdj, xiAdj] (double theta)
    {
        double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
        return fastcdfExponentiation(u);
    },
    -xiAdj, M_PI_2);
}

double StableDistribution::cdfForCommonExponent(double x) const
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

double StableDistribution::F(const double & x) const
{
    switch (distributionId) {
    case NORMAL:
        return cdfNormal(x);
    case CAUCHY:
        return cdfCauchy(x);
    case LEVY:
        return (beta > 0) ? cdfLevy(x) : cdfLevyCompl(2 * mu - x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x);
    case COMMON:
        return cdfForCommonExponent(x);
    default:
        return NAN; /// unexpected return
    }
}

double StableDistribution::S(const double & x) const
{
    switch (distributionId) {
    case NORMAL:
        return cdfNormalCompl(x);
    case CAUCHY:
        return cdfCauchyCompl(x);
    case LEVY:
        return (beta > 0) ? cdfLevyCompl(x) : cdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return 1.0 - cdfForUnityExponent(x);
    case COMMON:
        return 1.0 - cdfForCommonExponent(x);
    default:
        return NAN; /// unexpected return
    }
}

double StableDistribution::variateForUnityExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::StandardVariate();
    double pi_2pBetaU = M_PI_2 + beta * U;
    double Y = W * std::cos(U) / pi_2pBetaU;
    double X = std::log(Y);
    X += logGammaPi_2;
    X *= -beta;
    X += pi_2pBetaU * std::tan(U);
    X *= M_2_PI;
    return mu + gamma * X;
}

double StableDistribution::variateForCommonExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::StandardVariate();
    double alphaUpxi = alpha * (U + xi);
    double X = std::sin(alphaUpxi);
    double W_adj = W / std::cos(U - alphaUpxi);
    X *= W_adj;
    double R = omega - alphaInv * std::log(W_adj * std::cos(U));
    X *= std::exp(R);
    return mu + gamma * X;
}

double StableDistribution::variateForExponentEqualOneHalf() const
{
    double Z1 = NormalRand::StandardVariate(), Z2 = NormalRand::StandardVariate();
    double temp1 = (1.0 + beta) / Z1, temp2 = (1.0 - beta) / Z2;
    double var = temp1 - temp2;
    var *= temp1 + temp2;
    var *= 0.25;
    return mu + gamma * var;
}

double StableDistribution::Variate() const
{
    switch (distributionId) {
    case NORMAL:
        return mu + M_SQRT2 * gamma * NormalRand::StandardVariate();
    case CAUCHY:
        return CauchyRand::Variate(mu, gamma);
    case LEVY:
        return mu + RandMath::sign(beta) * LevyRand::Variate(0, gamma);
    case UNITY_EXPONENT:
        return variateForUnityExponent();
    case COMMON:
        return (alpha == 0.5) ? variateForExponentEqualOneHalf() : variateForCommonExponent();
    default:
        return NAN; /// unexpected return
    }
}

void StableDistribution::Sample(std::vector<double> &outputData) const
{
    switch (distributionId) {
    case NORMAL: {
        double stdev = M_SQRT2 * gamma;
        for (double &var : outputData)
            var = mu + stdev * NormalRand::StandardVariate();
    }
        break;
    case CAUCHY: {
        for (double &var : outputData)
            var = CauchyRand::Variate(mu, gamma);
    }
        break;
    case LEVY: {
        if (beta > 0) {
            for (double &var : outputData)
                var = LevyRand::Variate(mu, gamma);
        }
        else {
            for (double &var : outputData)
                var = mu - gamma * LevyRand::StandardVariate();
        }
    }
        break;
    case UNITY_EXPONENT: {
        for (double &var : outputData)
            var = variateForUnityExponent();
    }
        break;
    case COMMON: {
        if (alpha == 0.5) {
            for (double &var : outputData)
                var = variateForExponentEqualOneHalf();
        }
        else {
            for (double &var : outputData)
                var = variateForCommonExponent();
        }
    }
        break;
    default:
        break;
    }
}

double StableDistribution::Mean() const
{
    if (alpha > 1)
        return mu;
    if (beta == 1)
        return INFINITY;
    return (beta == -1) ? -INFINITY : NAN;
}

double StableDistribution::Variance() const
{
    return (distributionId == NORMAL) ? 2 * gamma * gamma : INFINITY;
}

double StableDistribution::Mode() const
{
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionId == NORMAL)
        return mu;
    if (distributionId == LEVY)
        return mu + beta * gamma / 3.0;
    return ContinuousDistribution::Mode();
}

double StableDistribution::Median() const
{
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionId == NORMAL)
        return mu;
    return ContinuousDistribution::Median();
}

double StableDistribution::Skewness() const
{
    return (distributionId == NORMAL) ? 0 : NAN;
}

double StableDistribution::ExcessKurtosis() const
{
    return (distributionId == NORMAL) ? 0 : NAN;
}

double StableDistribution::quantileNormal(double p) const
{
    return mu - 2 * gamma * RandMath::erfcinv(2 * p);
}

double StableDistribution::quantileNormal1m(double p) const
{
    return mu + 2 * gamma * RandMath::erfcinv(2 * p);
}

double StableDistribution::quantileCauchy(double p) const
{
    return mu - gamma / std::tan(M_PI * p);
}

double StableDistribution::quantileCauchy1m(double p) const
{
    return mu + gamma / std::tan(M_PI * p);
}

double StableDistribution::quantileLevy(double p) const
{
    double y = RandMath::erfcinv(p);
    return mu + 0.5 * gamma / (y * y);
}

double StableDistribution::quantileLevy1m(double p) const
{
    double y = RandMath::erfinv(p);
    return mu + 0.5 * gamma / (y * y);
}

double StableDistribution::quantileImpl(double p) const
{
    switch (distributionId) {
    case NORMAL:
        return quantileNormal(p);
    case CAUCHY:
        return quantileCauchy(p);
    case LEVY:
        return (beta > 0) ? quantileLevy(p) : 2 * mu - quantileLevy1m(p);
    default:
        return ContinuousDistribution::quantileImpl(p);
    }
}

double StableDistribution::quantileImpl1m(double p) const
{
    switch (distributionId) {
    case NORMAL:
        return quantileNormal1m(p);
    case CAUCHY:
        return quantileCauchy1m(p);
    case LEVY:
        return (beta > 0) ? quantileLevy1m(p) : 2 * mu - quantileLevy(p);
    default:
        return ContinuousDistribution::quantileImpl1m(p);
    }
}

std::complex<double> StableDistribution::cfNormal(double t) const
{
    double gammaT = gamma * t;
    std::complex<double> y(-gammaT * gammaT, mu * t);
    return std::exp(y);
}

std::complex<double> StableDistribution::cfCauchy(double t) const
{
    std::complex<double> y(-gamma * t, mu * t);
    return std::exp(y);
}

std::complex<double> StableDistribution::cfLevy(double t) const
{
    std::complex<double> y(0.0, -2 * gamma * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

std::complex<double> StableDistribution::CFImpl(double t) const
{
    double x = 0;
    switch (distributionId) {
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

std::string StableRand::Name() const
{
    return "Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

std::string HoltsmarkRand::Name() const
{
    return "Holtsmark("
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}


std::string LandauRand::Name() const
{
    return "Landau("
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}
