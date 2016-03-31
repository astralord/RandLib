#include "StableRand.h"
#include "LevyRand.h"
#include "CauchyRand.h"
#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"

StableRand::StableRand(double exponent, double skewness, double scale, double location)
{
    setParameters(exponent, skewness, scale, location);
}

std::string StableRand::name()
{
    return "Stable("
            + toStringWithPrecision(getExponent()) + ", "
            + toStringWithPrecision(getSkewness()) + ", "
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}

void StableRand::setParameters(double exponent, double skewness, double scale, double location)
{
    alpha = std::min(exponent, 2.0);
    if (alpha <= 0)
        alpha = 2.0;
    alphaInv = 1.0 / alpha;
    alpha_alpham1 = alpha / (alpha - 1.0);
    alpham1Inv = alpha_alpham1 - 1.0;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    /// Should be cautious, known distributions in priority
    if (RandMath::areClose(alpha, 2))
        alpha = 2;
    else if (RandMath::areClose(alpha, 1))
        alpha = 1;
    else if (RandMath::areClose(alpha, 0.5))
        alpha = 0.5;

    if (RandMath::areClose(beta, 1))
        beta = 1;
    else if (RandMath::areClose(beta, -1))
        beta = -1;

    if (alpha == 1 && beta != 0)
        pdfCoef = 0.5 / beta;
    if (alpha != 1 && alpha != 2 &&
        !(alpha == 0.5 && std::fabs(beta) == 1)) /// Common case
    {
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = std::pow(1 + B * B, .5 * alphaInv);
        B = std::atan(B);
        xi = alphaInv * B;
        integrandCoef = alpham1Inv * std::log(std::cos(B));
    }

    setScale(scale);
    setLocation(location);
}

void StableRand::setLocation(double location)
{
    mu = location;
}

void StableRand::setScale(double scale)
{
    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    if (alpha == 1)
    {
        if (beta != 0) /// not Cauchy
            logSigma = std::log(sigma);
    }
    else if (alpha == 2) /// Normal
        pdfCoef = M_SQRT1_2 / sigma;
    else if (alpha == 0.5 && std::fabs(beta) == 1) /// +/- Levy
        pdfCoef = M_1_SQRT2PI * std::sqrt(sigma);
    else
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

double StableRand::integrandAuxForAlphaEqualOne(double theta, double xAdj) const
{
    if (theta >= M_PI_2)
        return (beta > 0) ? 1e6 : -1e6;
    if (theta <= -M_PI_2)
        return (beta > 0) ? -1e6 : 1e6;
    if (theta == 0.0)
        return xAdj;
    double cosTheta = std::cos(theta);
    /// if theta ~ +-pi / 2
    if (std::fabs(cosTheta) < MIN_POSITIVE)
        return -1.0; // -INF
    double thetaAdj = (M_PI_2 + beta * theta) / cosTheta;
    double u = std::log(M_2_PI * thetaAdj);
    u += thetaAdj * std::sin(theta) / beta;
    if (std::isinf(u) || std::isnan(u))
    {
        if (theta < 0.0)
            return (beta > 0) ? -1e6 : 1e6;
        return (beta > 0) ? 1e6 : -1e6;
    }
    return u + xAdj;
}

double StableRand::integrandForAlphaEqualOne(double theta, double xAdj) const
{
    if (std::fabs(theta) >= M_PI_2)
        return 0.0;
    double u = integrandAuxForAlphaEqualOne(theta, xAdj);
    u = std::exp(u);
    u *= std::exp(-u);
    if (std::isinf(u) || std::isnan(u))
        return 0.0;
    return u;
}

double StableRand::pdfForAlphaEqualOne(double x) const
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize

    double xAdj = -M_PI * x * pdfCoef;

    /// find peak of integrand
    double theta0 = 0;

    std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForAlphaEqualOne, this, std::placeholders::_1, xAdj);
    RandMath::findRoot(funPtr, -M_PI_2, M_PI_2, theta0);

    /// Sanity check
    if (std::fabs(theta0) >= M_PI_2)
        theta0 = 0.0;

    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForAlphaEqualOne, this, std::placeholders::_1, xAdj);

    /// if theta0 is too close to +/- pi/2 then we can still underestimate the integral
    int maxRecursionDepth = 10;
    if (M_PI_2 - std::fabs(theta0) < 0.1)
        maxRecursionDepth = 20;
    else if (M_PI_2 - std::fabs(theta0) < 0.2)
        maxRecursionDepth = 15;

    double int1 = RandMath::integral(integrandPtr, -M_PI_2, theta0, 1e-11, maxRecursionDepth);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2, 1e-11, maxRecursionDepth);

    return std::fabs(pdfCoef) * (int1 + int2) / sigma;
}

double StableRand::integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const
{
    if (theta >= M_PI_2)
        return alpha < 1 ? 1e6 : -1e6; // +/-INF
    if (theta <= -M_PI_2 || theta <= -xiAdj)
        return alpha < 1 ? -1e6 : +1e6; // +/-INF
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
            return alpha < 1 ? -1e6 : +1e6; // +/-INF
        return alpha < 1 ? 1e6 : -1e6; // +/-INF
    }
    return integrandCoef + xAdj + y;
}

double StableRand::integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const
{
    double u = integrandAuxForCommonAlpha(theta, xAdj, xiAdj);
    u = std::exp(u);
    u *= std::exp(-u);
    if (std::isinf(u) || std::isnan(u))
        return 0.0;
    return u;
}

double StableRand::pdfForCommonAlpha(double x) const
{
    x = (x - mu) / sigma; /// Standardize

    // TODO: elaborate dangerous values - alpha close to 1 and zero
    if (std::fabs(x) < 1e-4) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double numerator = std::tgamma(1 + alphaInv);
        numerator *= std::cos(xi);
        double denominator = sigma * std::pow(1 + zeta * zeta, .5 * alphaInv);
        double y0 =  M_1_PI * numerator / denominator; /// f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 1.1e-4 : -1.1e-4;
        double y1 = pdfForCommonAlpha(mu + sigma * b);
        return RandMath::linearInterpolation(0, b, y0, y1, x); // TODO: find a way to do better interpolation than linear
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
        std::function<double (double)> funPtr = std::bind(&StableRand::integrandAuxForCommonAlpha, this, std::placeholders::_1, xAdj, xiAdj);
        RandMath::findRoot(funPtr, -xiAdj, M_PI_2, theta0);
    }

    /// calculate two integrals
    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForCommonAlpha, this, std::placeholders::_1, xAdj, xiAdj);
    double int1 = RandMath::integral(integrandPtr, -xiAdj, theta0);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2);
    return pdfCoef * (int1 + int2) / x;
}

double StableRand::f(double x) const
{
    /// Check all 'good' cases
    if (alpha == 2)
        return pdfNormal(x);
    if (alpha == 1 && beta == 0)
        return pdfCauchy(x);
    if (alpha == .5 && beta == 1)
        return pdfLevy(x);
    if (alpha == .5 && beta == -1)
        return pdfLevy(-x);

    /// Now check the others
    return (alpha == 1) ? pdfForAlphaEqualOne(x) : pdfForCommonAlpha(x);
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

double StableRand::cdfForCommonAlpha(double x) const
{
    x = (x - mu) / sigma; /// Standardize
    
    if (std::fabs(x) < 1e-4) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double y0 = 0.5 - M_1_PI * xi; /// f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 1.1e-4 : -1.1e-4;
        double y1 = cdfForCommonAlpha(mu + sigma * b);
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
        double u = integrandAuxForCommonAlpha(theta, xAdj, xiAdj);
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
    
double StableRand::cdfForAlphaEqualOne(double x) const
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize

    double xAdj = -M_PI * x * pdfCoef;

    double y = RandMath::integral([this, xAdj] (double theta)
    {
        double u = integrandAuxForAlphaEqualOne(theta, xAdj);
        u = std::exp(u);
        return std::exp(-u);
    },
    -M_PI_2, M_PI_2);

    y *= M_1_PI;
    return (beta > 0) ? y : 1 - y;
}

double StableRand::F(double x) const
{
    /// Check all 'good' cases    
    if (alpha == 2)
        return cdfNormal(x);
    if (alpha == 1 && beta == 0)
        return cdfCauchy(x);
    if (alpha == .5 && beta == 1)
        return cdfLevy(x);
    if (alpha == .5 && beta == -1)
        return 1.0 - cdfLevy(-x);

    /// Now check the others
    return (alpha == 1) ? cdfForAlphaEqualOne(x) : cdfForCommonAlpha(x);
}

double StableRand::variateForCommonAlpha() const
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

double StableRand::variateForAlphaEqualOne() const
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
    /// Check all 'good' cases
    if (alpha == 2)
        return NormalRand::variate(mu, sigma);
    if (alpha == 1 && beta == 0)
        return CauchyRand::variate(mu, sigma);
    if (alpha == .5 && beta == 1)
        return LevyRand::variate(mu, sigma);
    if (alpha == .5 && beta == -1)
        return -LevyRand::variate(mu, sigma);

    /// Now check the others
    return (alpha == 1) ? variateForAlphaEqualOne() : variateForCommonAlpha();
}

void StableRand::sample(std::vector<double> &outputData) const
{
    /// Check all 'good' cases
    if (alpha == 2) {
        for (double &var : outputData)
            var = NormalRand::variate(mu, sigma);
    }
    else if (alpha == 1 && beta == 0) {
        for (double &var : outputData)
            var = CauchyRand::variate(mu, sigma);
    }
    else if (alpha == .5 && beta == 1) {
        for (double &var : outputData)
            var = LevyRand::variate(mu, sigma);
    }
    else if (alpha == .5 && beta == -1) {
        for (double &var : outputData)
            var = -LevyRand::variate(mu, sigma);
    }
    else if (alpha == 1) {
        for (double &var : outputData)
            var = variateForAlphaEqualOne();
    }
    else {
        for (double &var : outputData)
            var = variateForCommonAlpha();
    }
}

std::complex<double> StableRand::psi(double t) const
{
    double x = (alpha == 1) ? beta * M_2_PI * log(std::fabs(t)) : zeta;
    if (t > 0)
        x = -x;
    double re = std::pow(std::fabs(sigma * t), alpha);
    return std::complex<double>(re, re * x + mu * t);
}

std::complex<double> StableRand::CF(double t) const
{
    return std::exp(-psi(t));
}

double StableRand::Mean() const
{
    if (alpha > 1)
        return mu;
    return (alpha == 1.0) ? NAN : INFINITY;
}

double StableRand::Variance() const
{
    return (alpha == 2) ? 2 * sigma * sigma : INFINITY;
}

double StableRand::Median() const
{
    return (beta == 0) ? mu : ContinuousDistribution::Median();
}

double StableRand::Mode() const
{
    return (beta == 0) ? mu : ContinuousDistribution::Mode();
}

double StableRand::Skewness() const
{
    return (alpha == 2) ? 0 : NAN;
}

double StableRand::ExcessKurtosis() const
{
    return (alpha == 2) ? 0 : NAN;
}


std::string HoltsmarkRand::name()
{
    return "Holtsmark("
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}


std::string LandauRand::name()
{
    return "Landau("
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}
