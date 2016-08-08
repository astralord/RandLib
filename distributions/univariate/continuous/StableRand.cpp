#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"


StableRand::StableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location)
{
    setParameters(exponent, skewness);
    setScale(scale);
    setLocation(location);
}

std::string StableRand::name() const
{
    return "Stable("
            + toStringWithPrecision(getExponent()) + ", "
            + toStringWithPrecision(getSkewness()) + ", "
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}

void StableRand::setParameters(double exponent, double skewness)
{
    LimitingDistribution::setParameters(exponent, skewness);

    /// set id of distribution
    if (RandMath::areClose(alpha, 2.0))
        distributionId = NORMAL;
    else if (RandMath::areClose(alpha, 1.0))
        distributionId = (beta == 0.0) ? CAUCHY : UNITY_EXPONENT;
    else if (RandMath::areClose(alpha, 0.5) && RandMath::areClose(std::fabs(beta), 1.0))
        distributionId = LEVY;
    else
        distributionId = COMMON;

    alpha_alpham1 = alpha / (alpha - 1.0);
    alpham1Inv = alpha_alpham1 - 1.0;

    if (distributionId == UNITY_EXPONENT)
        pdfCoef = 0.5 / beta;
    else if (distributionId == COMMON) {
        xi = alphaInv * B;
        integrandCoef = alpham1Inv * std::log(std::cos(B));
    }
}

void StableRand::setScale(double scale)
{
    LimitingDistribution::setScale(scale);
    if (distributionId == NORMAL)
        pdfCoef = 0.5 / sigma;
    else if (distributionId == LEVY)
        pdfCoef = M_1_SQRT2PI * std::sqrt(sigma);
    else if (distributionId == COMMON)
        pdfCoef = M_1_PI * alpha / (std::fabs(1 - alpha) * sigma);
}

double StableRand::pdfNormal(double x) const
{
    double y = x - mu;
    y *= pdfCoef;
    y *= y;
    y = std::exp(-y);
    return M_1_SQRTPI * pdfCoef * y;
}

double StableRand::pdfCauchy(double x) const
{
    double y = x - mu;
    y *= y;
    y /= sigma;
    y += sigma;
    return M_1_PI / y;
}

double StableRand::pdfLevy(double x) const
{
    if (x <= mu)
        return 0;
    double xInv = 1.0 / (x - mu);
    double y = -0.5 * sigma * xInv;
    y = std::exp(y);
    y *= xInv;
    y *= std::sqrt(xInv);
    return pdfCoef * y;
}

double StableRand::integrandAuxForUnityExponent(double theta, double xAdj) const
{
    if (theta > M_PI_2 || RandMath::areClose(theta, M_PI_2))
        return (beta > 0) ? BIG_NUMBER : -BIG_NUMBER;
    if (theta < -M_PI_2 || RandMath::areClose(theta, -M_PI_2))
        return (beta > 0) ? -BIG_NUMBER : BIG_NUMBER;
    if (theta == 0.0)
        return xAdj;
    double cosTheta = std::cos(theta);
    double thetaAdj = (M_PI_2 + beta * theta) / cosTheta;
    double u = std::log(M_2_PI * thetaAdj);
    u += thetaAdj * std::sin(theta) / beta;
    if (std::isinf(u) || std::isnan(u))
    {
        if (theta < 0.0)
            return (beta > 0) ? -BIG_NUMBER : BIG_NUMBER;
        return (beta > 0) ? BIG_NUMBER : -BIG_NUMBER;
    }
    return u + xAdj;
}

double StableRand::integrandForUnityExponent(double theta, double xAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    double u = integrandAuxForUnityExponent(theta, xAdj);
    if (std::fabs(u) >= BIG_NUMBER)
        return 0.0;
    u = std::exp(u);
    u *= std::exp(-u);
    if (std::isinf(u) || std::isnan(u))
        return 0.0;
    return u;
}

double StableRand::pdfForUnityExponent(double x) const
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize

    double xAdj = -M_PI * x * pdfCoef;

    /// find peak of integrand
    double theta0 = 0;

    std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForUnityExponent, this, std::placeholders::_1, xAdj);
    RandMath::findRoot(funPtr, -M_PI_2, M_PI_2, theta0);

    /// Sanity check
    if (std::fabs(theta0) >= M_PI_2)
        theta0 = 0.0;

    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForUnityExponent, this, std::placeholders::_1, xAdj);

    /// if theta0 is too close to +/- pi/2 then we can still underestimate the integral
    int maxRecursionDepth = 10;
    double closeness = M_PI_2 - std::fabs(theta0);
    if (closeness < 0.1)
        maxRecursionDepth = 20;
    else if (closeness < 0.2)
        maxRecursionDepth = 15;

    double int1 = RandMath::integral(integrandPtr, -M_PI_2, theta0, 1e-11, maxRecursionDepth);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2, 1e-11, maxRecursionDepth);

    return std::fabs(pdfCoef) * (int1 + int2) / sigma;
}

double StableRand::integrandAuxForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    if (theta > M_PI_2 || RandMath::areClose(theta, M_PI_2))
        return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
    if (theta < -M_PI_2 || RandMath::areClose(theta, -M_PI_2) || theta < -xiAdj || RandMath::areClose(theta, -xiAdj))
        return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    double cosTheta = std::cos(theta);
    double y = std::cos(thetaAdj - theta) / cosTheta;
    y = std::log(y);
    y += alpha_alpham1 * std::log(cosTheta / sinThetaAdj);
    if (std::isinf(y) || std::isnan(y))
    {
        /// we got numerical error, need to investigate to which extreme point we are closer
        if (theta < 0.5 * (M_PI_2 - xiAdj))
            return alpha < 1 ? -BIG_NUMBER : BIG_NUMBER;
        return alpha < 1 ? BIG_NUMBER : -BIG_NUMBER;
    }
    return integrandCoef + xAdj + y;
}

double StableRand::integrandForCommonExponent(double theta, double xAdj, double xiAdj) const
{
    double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
    if (std::fabs(u) >= BIG_NUMBER)
        return 0.0;
    u = std::exp(u);
    u *= std::exp(-u);
    if (std::isinf(u) || std::isnan(u))
        return 0.0;
    return u;
}

double StableRand::pdfForCommonExponent(double x) const
{
    x = (x - mu) / sigma; /// Standardize

    // TODO: elaborate dangerous values - alpha close to 1 and zero
    if (std::fabs(x) < 1e-4) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double y0 = std::lgamma(1 + alphaInv);
        y0 -= 0.5 * alphaInv * std::log1p(zeta * zeta);
        y0 = std::exp(y0);
        y0 *= std::cos(xi);
        y0 /= (M_PI * sigma);
        /// Now y0 = f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 1.1e-4 : -1.1e-4;
        double y1 = pdfForCommonExponent(mu + sigma * b);
        return RandMath::linearInterpolation(0, b, y0, y1, x);
    }

    double xiAdj = xi; /// +- xi
    if (x > 0)
    {
        if (alpha < 1 && beta == -1)
            return 0;
    }
    else
    {
        if (alpha < 1 && beta == 1)
            return 0;
        x = -x;
        xiAdj = -xi;
    }

    double xAdj = alpha_alpham1 * std::log(x);

    /// find the peak of the integrand
    double theta0 = 0.5 * (M_PI_2 - xiAdj);
    if (-xiAdj < M_PI_2)
    {
        std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
        RandMath::findRoot(funPtr, -xiAdj, M_PI_2, theta0);
    }

    /// calculate two integrals
    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForCommonExponent, this, std::placeholders::_1, xAdj, xiAdj);
    double int1 = RandMath::integral(integrandPtr, -xiAdj, theta0);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2);
    return pdfCoef * (int1 + int2) / x;
}

double StableRand::f(double x) const
{
    switch (distributionId) {
    case NORMAL:
        return pdfNormal(x);
    case CAUCHY:
        return pdfCauchy(x);
    case LEVY:
        return (beta > 0) ? pdfLevy(x) : pdfLevy(-x);
    case UNITY_EXPONENT:
        return pdfForUnityExponent(x);
    case COMMON:
        return pdfForCommonExponent(x);
    default:
        return NAN; /// unexpected return
    }
}

double StableRand::cdfNormal(double x) const
{
    double y = x - mu;
    y *= pdfCoef;
    y = std::erf(y);
    ++y;
    return 0.5 * y;
}

double StableRand::cdfCauchy(double x) const
{
    double y = x - mu;
    y /= sigma;
    y = std::atan(y);
    y *= M_1_PI;
    return y + 0.5;
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

double StableRand::cdfForCommonExponent(double x) const
{
    x = (x - mu) / sigma; /// Standardize
    
    if (std::fabs(x) < 1e-4) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double y0 = 0.5 - M_1_PI * xi; /// f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 1.1e-4 : -1.1e-4;
        double y1 = cdfForCommonExponent(mu + sigma * b);
        return RandMath::linearInterpolation(0, b, y0, y1, x);
    }
    
    double xiAdj = xi; /// +- xi
    double xAdj;
    if (x > 0)
    {
        if (alpha < 1 && beta == -1)
            return 1.0;
        xAdj = alpha_alpham1 * std::log(x);
    }
    else
    {
        if (alpha < 1 && beta == 1)
            return 0.0;
        xAdj = alpha_alpham1 * std::log(-x);
        xiAdj = -xi;
    }

    double y = RandMath::integral([this, xAdj, xiAdj] (double theta)
    {
        double u = integrandAuxForCommonExponent(theta, xAdj, xiAdj);
        if (u >= BIG_NUMBER)
            return 0.0;
        else if (u <= -BIG_NUMBER)
            return 1.0;
        u = std::exp(u);
        return std::exp(-u);
    },
    -xiAdj, M_PI_2);

    if (alpha > 1)
        y = 1.0 - y * M_1_PI;
    else
        y = 0.5 + (y - xiAdj) * M_1_PI;

    return (x > 0) ? y : 1 - y;
}
    
double StableRand::cdfForUnityExponent(double x) const
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize
    double xAdj = -M_PI * x * pdfCoef;
    double y = RandMath::integral([this, xAdj] (double theta)
    {
        double u = integrandAuxForUnityExponent(theta, xAdj);
        if (u >= BIG_NUMBER)
            return 0.0;
        else if (u <= -BIG_NUMBER)
            return 1.0;
        u = std::exp(u);
        return std::exp(-u);
    },
    -M_PI_2, M_PI_2);
    y *= M_1_PI;
    return (beta > 0) ? y : 1 - y;
}

double StableRand::F(double x) const
{
    switch (distributionId) {
    case NORMAL:
        return cdfNormal(x);
    case CAUCHY:
        return cdfCauchy(x);
    case LEVY:
        return (beta > 0) ? cdfLevy(x) : 1.0 - cdfLevy(-x);
    case UNITY_EXPONENT:
        return cdfForUnityExponent(x);
    case COMMON:
        return cdfForCommonExponent(x);
    default:
        return NAN; /// unexpected return
    }
}

double StableRand::variateForCommonExponent() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::standardVariate();
    double alphaUB = alpha * U + B;
    double X = S * std::sin(alphaUB); /// S * sin(alpha * V + B)
    double W_adj = W / std::cos(U - alphaUB);
    X *= W_adj; /// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B)
    X *= std::pow(W_adj * std::cos(U), -alphaInv); /// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B) /
                                                   /// ((W * cos(V) / cos((1 - alpha) * V - B)) ^ (1 / alpha))
    return mu + sigma * X;
}

double StableRand::variateForUnityExponent() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::standardVariate();
    double pi_2pBetaU = M_PI_2 + beta * U;
    double X = logSigma;
    X -= std::log(M_PI_2 * W * std::cos(U) / pi_2pBetaU);
    X *= beta;
    X += pi_2pBetaU * std::tan(U);
    X *= M_2_PI;
    return mu + sigma * X;
}

double StableRand::variate() const
{
    switch (distributionId) {
    case NORMAL:
        return NormalRand::variate(mu, M_SQRT2 * sigma);
    case CAUCHY:
        return CauchyRand::variate(mu, sigma);
    case LEVY:
        return RandMath::sign(beta) * LevyRand::variate(mu, sigma);
    case UNITY_EXPONENT:
        return variateForUnityExponent();
    case COMMON:
        return variateForCommonExponent();
    default:
        return NAN; /// unexpected return
    }
}

void StableRand::sample(std::vector<double> &outputData) const
{
    switch (distributionId) {
    case NORMAL: {
        for (double &var : outputData)
            var = NormalRand::variate(mu, M_SQRT2 * sigma);
    }
        break;
    case CAUCHY: {
        for (double &var : outputData)
            var = CauchyRand::variate(mu, sigma);
    }
        break;
    case LEVY: {
        if (beta > 0) {
            for (double &var : outputData)
                var = LevyRand::variate(mu, sigma);
        }
        else {
            for (double &var : outputData)
                var = -LevyRand::variate(mu, sigma);
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
    return (distributionId == NORMAL)  ? 2 * sigma * sigma : INFINITY;
}

std::complex<double> StableRand::CF(double t) const
{
    return (t == 0) ? 1.0 : std::exp(-psi(t));
}

double StableRand::Mode() const
{
    switch (distributionId) {
    case NORMAL:
    case CAUCHY:
        return mu;
    case LEVY:
        return mu + RandMath::sign(beta) * sigma / 3.0;
    default:
        return ContinuousDistribution::Mode();
    }
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

std::string HoltsmarkRand::name() const
{
    return "Holtsmark("
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}


std::string LandauRand::name() const
{
    return "Landau("
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}
