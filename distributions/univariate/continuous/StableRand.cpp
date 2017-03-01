#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"


StableRand::StableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location)
{
    SetParameters(exponent, skewness);
    SetScale(scale);
    SetLocation(location);
}

std::string StableRand::Name() const
{
    return "Stable("
            + toStringWithPrecision(GetExponent()) + ", "
            + toStringWithPrecision(GetSkewness()) + ", "
            + toStringWithPrecision(GetScale()) + ", "
            + toStringWithPrecision(GetLocation()) + ")";
}

void StableRand::SetParameters(double exponent, double skewness)
{
    LimitingDistribution::SetParameters(exponent, skewness);

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
    alpham1Inv = alpha_alpham1 - 1.0;

    if (distributionId == NORMAL) {
        pdfCoef = std::log(2 * sigma) + 0.5 * M_LNPI;
    }
    else if (distributionId == LEVY) {
        pdfCoef = std::log(sigma) - M_LN2 - M_LNPI;
    }
    else if (distributionId == UNITY_EXPONENT) {
        pdfCoef = 0.5 / (sigma * std::fabs(beta));
        /// pdfXLimit is such k that f(x) < 1e-4 for |x| > k
        pdfXLimit = std::sqrt(2e4 / M_PI * M_E);
    }
    else if (distributionId == COMMON) {
        xi = beta * std::tan(M_PI_2 * alpha);
        zeta = -xi;
        omega = 0.5 * alphaInv * std::log1p(zeta * zeta);
        xi = alphaInv * RandMath::atan(xi);
        pdfCoef = M_1_PI * std::fabs(alpha_alpham1) / sigma;
        /// pdfXLimit is such k that for x > k we use asymptotic expansion
        pdfXLimit = 3.0 / (1.0 + alpha) * M_LN10;
        /// define boundaries of region near 0, where we use series expansion
        if (alpha <= ALMOST_TWO) {
            seriesZeroParams.first = std::round(std::min(alpha * alpha * 40 + 1, 10.0));
            seriesZeroParams.second = -(1.5 / alpha + 0.5) * M_LN10; /// corresponds to boundaries from 10^(-15.5) to ~ 0.056
        }
        else {
            seriesZeroParams.first = 85;
            seriesZeroParams.second = M_LN2 + M_LN3; /// corresponds to 6
        }
    }
}

void StableRand::SetScale(double scale)
{
    LimitingDistribution::SetScale(scale);
    if (distributionId == NORMAL)
        pdfCoef = std::log(2 * sigma) + 0.5 * M_LNPI;
    else if (distributionId == LEVY)
        pdfCoef = std::log(sigma) - M_LN2 - M_LNPI;
    else if (distributionId == COMMON)
        pdfCoef = M_1_PI * std::fabs(alpha_alpham1) / sigma;
}

double StableRand::pdfNormal(double x) const
{
    return std::exp(logpdfNormal(x));
}

double StableRand::logpdfNormal(double x) const
{
    double y = x - mu;
    y *= 0.5 / sigma;
    y *= y;
    y += pdfCoef;
    return -y;
}

double StableRand::pdfCauchy(double x) const
{
    double y = x - mu;
    y *= y;
    y /= sigma;
    y += sigma;
    return M_1_PI / y;
}

double StableRand::logpdfCauchy(double x) const
{
    return std::log(pdfCauchy(x));
}

double StableRand::pdfLevy(double x) const
{
    return (x <= mu) ? 0.0 : std::exp(logpdfLevy(x));
}

double StableRand::logpdfLevy(double x) const
{
    double x0 = x - mu;
    if (x0 <= 0.0)
        return -INFINITY;
    double y = sigma / x0;
    y += 3 * std::log(x0);
    y -= pdfCoef;
    return -0.5 * y;
}

double StableRand::fastpdfExponentiation(double u)
{
    if (u > 5 || u < -150)
        return 0.0;
    return (u < -50) ? std::exp(u) : std::exp(u - std::exp(u));
}

double StableRand::limitCaseForIntegrandAuxForUnityExponent(double theta, double xAdj) const
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

double StableRand::integrandAuxForUnityExponent(double theta, double xAdj) const
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

double StableRand::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    return fastpdfExponentiation(u);
}

double StableRand::pdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / sigma;
    double xAdj = -M_PI_2 * xSt / beta - logsigmaPi_2;

    // TODO: this limit is not available anymore, re-check
    if (std::fabs(xSt) > pdfXLimit) {
        double y = 1.0 / (M_PI * xSt * xSt);
        y *= (xSt > 0) ? (1 + beta) : (1 - beta);
        return y;
    }

    /// We squize boudaries for too peaked integrands
    double boundary = RandMath::atan(M_2_PI * beta * (5.0 - xAdj));
    double upperBoundary = (beta > 0.0) ? boundary : M_PI_2;
    double lowerBoundary = (beta < 0.0) ? boundary : -M_PI_2;

    /// Find peak of the integrand
    double theta0 = 0;
    std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForUnityExponent, this, std::placeholders::_1, xAdj);
    RandMath::findRoot(funPtr, lowerBoundary, upperBoundary, theta0);

    /// Sanity check
    /// if we failed while looking for the peak position
    /// we set it in the middle between boundaries
    if (theta0 >= upperBoundary || theta0 <= lowerBoundary)
        theta0 = 0.5 * (upperBoundary + lowerBoundary);

    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForUnityExponent, this, std::placeholders::_1, xAdj);

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

double StableRand::pdfAtZero() const
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

double StableRand::pdfSeriesExpansionAtZero(double logX, double xiAdj, int k) const
{
    /// Calculate first term of the sum
    /// (if x = 0, only this term is non-zero)
    double y0 = pdfAtZero();
    double sum = 0.0;
    if (beta == 0.0) {
        /// Symmetric distribution
        for (int i = 1; i <= k; ++i)
        {
            int i2 = i + i;
            double term = std::lgamma((i2 + 1) / alpha);
            term += i2 * logX;
            term = std::exp(term);
            term /= RandMath::factorial(i2);
            sum += (i & 1) ? -term : term;
        }
    }
    else {
        /// Asymmetric distribution
        double coef = M_PI_2 + xiAdj;
        for (int i = 1; i <= k; ++i) {
            int ip1 = i + 1;
            double term = std::lgamma(ip1 * alphaInv);
            term += i * logX;
            term = std::exp(term - omega);
            term /= RandMath::factorial(i);
            term *= std::sin(ip1 * coef);
            sum += (i & 1) ? -term : term;
        }
    }
    return y0 + sum * M_1_PI / alpha;
}

double StableRand::pdfSeriesExpansionAtInf(double logX, double xiAdj, int k) const
{
    double coef = M_PI_2;
    if (beta != 0.0)
        coef += xiAdj;
    coef *= alpha;
    double sum = 0.0;
    for (int i = 1; i <= k; ++i) {
        double aux = i * alpha + 1.0;
        double term = std::lgamma(aux);
        term -= aux * logX;
        term = std::exp(term - omega);
        term /= RandMath::factorial(i);
        term *= std::sin(coef * i);
        sum += (i & 1) ? term : -term;
    }
    return M_1_PI * sum;
}

double StableRand::pdfTaylorExpansionTailNearCauchy(double x) const
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

double StableRand::limitCaseForIntegrandAuxForCommonExponent(double theta, double xiAdj) const
{
    /// We got numerical error, need to investigate to which extreme point we are closer
    if (theta < 0.5 * (M_PI_2 - xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
}

double StableRand::integrandAuxForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= M_PI_2 || theta <= -xiAdj)
        return limitCaseForIntegrandAuxForCommonExponent(theta, xiAdj);
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    double y = std::log(std::cos(theta));
    y -= alpha * std::log(sinThetaAdj);
    y *= alpham1Inv;
    y += std::log(std::cos(thetaAdj - theta));
    return std::isfinite(y) ? y + xAdj : limitCaseForIntegrandAuxForCommonExponent(theta, xiAdj);
}

double StableRand::integrandForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    if (theta <= -xiAdj)
        return 0.0;
    double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
    return fastpdfExponentiation(u);
}

double StableRand::pdfForCommonExponent(double x) const
{
    double xSt = (x - mu) / sigma; /// Standardize

    double absXSt = xSt;
    double xiAdj = xi; /// +- xi
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
    if (beta == 0.0 && alpha > 0.99 && alpha < 1.01)
        return pdfCauchy(x) + pdfTaylorExpansionTailNearCauchy(absXSt) / sigma;

    /// If x = 0, we know the analytic solution
    if (xSt == 0.0)
        return pdfAtZero() / sigma;

    double logAbsX = std::log(absXSt) - omega;

    /// If x is too close to 0, we do series expansion avoiding numerical problems
    if (logAbsX < seriesZeroParams.second)
        return pdfSeriesExpansionAtZero(logAbsX, xiAdj, seriesZeroParams.first) / sigma;

    /// If x is large enough we use tail approximation
    if (logAbsX > pdfXLimit && alpha <= ALMOST_TWO)
        return pdfSeriesExpansionAtInf(logAbsX, xiAdj, 10) / sigma;

    double xAdj = alpha_alpham1 * logAbsX;

    /// Search for the peak of the integrand
    double theta0;
    std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
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
    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
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
    tail *= (1.0 - 0.5 * alpha) / sigma;
    return std::max(tail, res);
}

double StableRand::f(const double & x) const
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

double StableRand::logf(const double & x) const
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

double StableRand::cdfNormal(double x) const
{
    double y = mu - x;
    y *= 0.5 / sigma;
    return 0.5 * std::erfc(y);
}

double StableRand::cdfNormalCompl(double x) const
{
    double y = x - mu;
    y *= 0.5 / sigma;
    return 0.5 * std::erfc(y);
}

double StableRand::cdfCauchy(double x) const
{
    double x0 = x - mu;
    x0 /= sigma;
    return 0.5 + M_1_PI * RandMath::atan(x0);
}

double StableRand::cdfCauchyCompl(double x) const
{
    double x0 = mu - x;
    x0 /= sigma;
    return 0.5 + M_1_PI * RandMath::atan(x0);
}

double StableRand::cdfLevy(double x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y += y;
    y = sigma / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

double StableRand::cdfLevyCompl(double x) const
{
    if (x <= mu)
        return 1.0;
    double y = x - mu;
    y += y;
    y = sigma / y;
    y = std::sqrt(y);
    return std::erf(y);
}

double StableRand::fastcdfExponentiation(double u)
{
    if (u > 5.0)
        return 0.0;
    else if (u < -150.0)
        return 1.0;
    double y = std::exp(u);
    return std::exp(-y);
}

double StableRand::cdfForUnityExponent(double x) const
{
    double xSt = (x - mu) / sigma;
    double xAdj = -M_PI_2 * xSt / beta - logsigmaPi_2;
    double y = M_1_PI * RandMath::integral([this, xAdj] (double theta)
    {
        double u = integrandAuxForUnityExponent(theta, xAdj);
        return fastcdfExponentiation(u);
    },
    -M_PI_2, M_PI_2);
    return (beta > 0) ? y : 1.0 - y;
}

double StableRand::cdfIntegralRepresentation(double absXSt, double xiAdj) const
{
    /// Here we assume that absXSt > 0
    double xAdj = alpha_alpham1 * (std::log(absXSt) - omega);
    return M_1_PI * RandMath::integral([this, xAdj, xiAdj] (double theta)
    {
        double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
        return fastcdfExponentiation(u);
    },
    -xiAdj, M_PI_2);
}

double StableRand::cdfForCommonExponent(double x) const
{
    double xSt = (x - mu) / sigma; /// Standardize
    
    if (std::fabs(xSt) < 1e-4) /// If we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double y0 = 0.5 - M_1_PI * xi; /// f(0)
        if (std::fabs(xSt) < MIN_POSITIVE)
            return y0;
        double b = (xSt > 0) ? 1.1e-4 : -1.1e-4;
        double y1 = cdfForCommonExponent(mu + sigma * b);
        return RandMath::linearInterpolation(0, b, y0, y1, xSt);
    }

    if (alpha > 1.0) {
        if (xSt > 0.0)
            return 1.0 - cdfIntegralRepresentation(xSt, xi);
        /// We use relation F(-x, xi) + F(x, -xi) = 1
        return cdfIntegralRepresentation(-xSt, -xi);
    }

    /// α < 1
    double temp = 0.5 - xi / M_PI;
    if (xSt > 0.0)
        return (beta == -1.0) ? 1.0 : temp + cdfIntegralRepresentation(xSt, xi);
    /// We use relation F(-x, xi) + F(x, -xi) = 1
    return (beta == 1.0) ? 0.0 : temp - cdfIntegralRepresentation(-xSt, -xi);
}

double StableRand::F(const double & x) const
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

double StableRand::S(const double & x) const
{
    switch (distributionId) {
    case NORMAL:
        return cdfNormalCompl(x);
    case CAUCHY:
        return cdfCauchyCompl(x);
    case LEVY:
        return (beta > 0) ? cdfLevyCompl(x) : cdfLevy(2 * mu - x);
    case UNITY_EXPONENT:
        return 1.0 - cdfForUnityExponent(x); // TODO: implement
    case COMMON:
        return 1.0 - cdfForCommonExponent(x); // TODO: implement
    default:
        return NAN; /// unexpected return
    }
}

double StableRand::variateForUnityExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::StandardVariate();
    double pi_2pBetaU = M_PI_2 + beta * U;
    double Y = W * std::cos(U) / pi_2pBetaU;
    double X = std::log(Y);
    X += logsigmaPi_2;
    X *= -beta;
    X += pi_2pBetaU * std::tan(U);
    X *= M_2_PI;
    return mu + sigma * X;
}

double StableRand::variateForCommonExponent() const
{
    double U = UniformRand::Variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::StandardVariate();
    double alphaUpxi = alpha * (U + xi);
    double X = std::sin(alphaUpxi);
    double W_adj = W / std::cos(U - alphaUpxi);
    X *= W_adj;
    double R = omega - alphaInv * std::log(W_adj * std::cos(U));
    X *= std::exp(R);
    return mu + sigma * X;
}

double StableRand::Variate() const
{
    switch (distributionId) {
    case NORMAL:
        return mu + M_SQRT2 * sigma* NormalRand::StandardVariate();
    case CAUCHY:
        return CauchyRand::Variate(mu, sigma);
    case LEVY:
        return mu + RandMath::sign(beta) * LevyRand::Variate(0, sigma);
    case UNITY_EXPONENT:
        return variateForUnityExponent();
    case COMMON:
        return variateForCommonExponent();
    default:
        return NAN; /// unexpected return
    }
}

void StableRand::Sample(std::vector<double> &outputData) const
{
    switch (distributionId) {
    case NORMAL: {
        double stdev = M_SQRT2 * sigma;
        for (double &var : outputData)
            var = mu + stdev * NormalRand::StandardVariate();
    }
        break;
    case CAUCHY: {
        for (double &var : outputData)
            var = CauchyRand::Variate(mu, sigma);
    }
        break;
    case LEVY: {
        if (beta > 0) {
            for (double &var : outputData)
                var = LevyRand::Variate(mu, sigma);
        }
        else {
            for (double &var : outputData)
                var = mu - sigma * LevyRand::StandardVariate();
        }
    }
        break;
    case UNITY_EXPONENT: {
        for (double &var : outputData)
            var = variateForUnityExponent();
    }
        break;
    case COMMON: {
        for (double &var : outputData)
            var = variateForCommonExponent();
    }
        break;
    default:
        break;
    }
}

double StableRand::Variance() const
{
    return (distributionId == NORMAL) ? 2 * sigma * sigma : INFINITY;
}

double StableRand::Mode() const
{
    /// For symmetric and normal distributions mode is μ
    if (beta == 0 || distributionId == NORMAL)
        return mu;
    if (distributionId == LEVY)
        return mu + beta * sigma / 3.0;
    return ContinuousDistribution::Mode();
}

double StableRand::Median() const
{
    return (beta == 0) ? mu : ContinuousDistribution::Median();
}

double StableRand::Skewness() const
{
    return (distributionId == NORMAL) ? 0 : NAN;
}

double StableRand::ExcessKurtosis() const
{
    return (distributionId == NORMAL)  ? 0 : NAN;
}

std::complex<double> StableRand::CFImpl(double t) const
{
    return std::exp(-psi(t));
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
