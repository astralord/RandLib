#include "StableRand.h"
#include <QDebug>

StableRand::StableRand(double exponent, double skewness, double scale, double location)
{
    setParameters(exponent, skewness, scale, location);
}

std::string StableRand::name()
{
    return "Stable("
            + toStringWithPrecision(getAlpha()) + ", "
            + toStringWithPrecision(getBeta()) + ", "
            + toStringWithPrecision(getSigma()) + ", "
            + toStringWithPrecision(getMu()) + ")";
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

    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    mu = location;

    /// Should be cautious, known distributions in priority
    if (RandMath::areEqual(alpha, 2))
        alpha = 2;
    else if (RandMath::areEqual(alpha, 1))
        alpha = 1;
    else if (RandMath::areEqual(alpha, 0.5))
        alpha = 0.5;

    if (RandMath::areEqual(beta, 1))
        beta = 1;
    else if (RandMath::areEqual(beta, -1))
        beta = -1;

    if (alpha == 2) /// X ~ Normal(mu, 2sigma^2)
    {
        N.setMean(mu);
        N.setSigma(sigma * M_SQRT2);
    }
    else if (alpha == 1)
    {
        if (beta == 0) /// X ~ Cauchy(mu, sigma)
        {
            C.setLocation(mu);
            C.setScale(sigma);
        }
        else /// just alpha == 1
        {
            logSigma = std::log(sigma);
            pdfCoef = 0.5 / beta;
        }
    }
    else if (alpha == .5 && std::fabs(beta) == 1) /// +/- X ~ Levy(mu, sigma)
    {
        L.setLocation(mu);
        L.setScale(sigma);
    }
    else /// Common case: alpha != 1
    {
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = std::pow(1 + B * B, .5 * alphaInv);
        B = std::atan(B);
        pdfCoef = M_1_PI * alpha / (std::fabs(1 - alpha) * sigma);
        xi = alphaInv * B;
        integrandCoef = std::pow(std::cos(B), alpham1Inv);
    }
}

double StableRand::integrandAuxForAlphaEqualOne(double theta, double xAdj) const
{
    double cosTheta = std::cos(theta);
    /// if theta ~ +-pi / 2
    if (std::fabs(cosTheta) < MIN_POSITIVE)
        return 0.0;
    double thetaAdj = (M_PI_2 + beta * theta) / cosTheta;
    double u = M_2_PI * thetaAdj;
    u *= std::exp(thetaAdj * std::sin(theta) / beta);
    return u * xAdj;
}

double StableRand::integrandForAlphaEqualOne(double theta, double xAdj) const
{
    double u = integrandAuxForAlphaEqualOne(theta, xAdj);
    return u * std::exp(-u);
}

double StableRand::integrandAuxForCommonAlpha(double theta, double xAdj, double xiAdj) const
{
    double thetaAdj = alpha * (theta + xiAdj);
    double sinThetaAdj = std::sin(thetaAdj);
    /// if theta ~ 0
    if (std::fabs(sinThetaAdj) < MIN_POSITIVE)
        return 0.0;
    double cosTheta = std::cos(theta);
    double y = cosTheta / sinThetaAdj;
    /// if theta ~ pi / 2
    if (std::fabs(y) < MIN_POSITIVE)
        return 0.0;
    y = std::pow(y, alpha_alpham1);
    y *= std::cos(thetaAdj - theta);
    y /= cosTheta;
    return integrandCoef * xAdj * y;
}

double StableRand::integrandForCommonAlpha(double theta, double xAdj, double xiAdj) const
{
    double u = integrandAuxForCommonAlpha(theta, xAdj, xiAdj);
    return u * std::exp(-u);
}

double StableRand::pdfForCommonAlpha(double x) const
{
    x = (x - mu) / sigma; /// Standardize

    if (std::fabs(x) < 0.01) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double numerator = std::tgamma(1 + alphaInv);
        numerator *= std::cos(xi);
        double denominator = sigma * std::pow(1 + zeta * zeta, .5 * alphaInv);
        double y0 =  M_1_PI * numerator / denominator; /// f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 0.011 : -0.011;
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

    double xAdj = std::pow(x, alpha_alpham1);

    /// find the peak of the integrand
    double theta0 = M_PI_4 - .5 * xiAdj;
    if (-xiAdj < almostPI_2)
    {
        RandMath::findRoot([this, xAdj, xiAdj] (double theta)
        {
            return integrandAuxForCommonAlpha(theta, xAdj, xiAdj) - 1.0;
        },
        -xiAdj, almostPI_2, theta0);
    }

    /// calculate two integrals
    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForCommonAlpha, this, std::placeholders::_1, xAdj, xiAdj);
    double int1 = RandMath::integral(integrandPtr, -xiAdj, theta0);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2);
    return pdfCoef * (int1 + int2) / x;
}

double StableRand::pdfForAlphaEqualOne(double x) const
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize

    double xAdj = std::exp(-M_PI * x * pdfCoef);

    /// find peak of integrand
    double theta0 = 0;
    if (beta > 0) {
        RandMath::findRoot([this, xAdj] (double theta)
        {
            return integrandAuxForAlphaEqualOne(theta, xAdj) - 1.0;
        },
        -M_PI_2, almostPI_2, theta0);
    }
    else {
        RandMath::findRoot([this, xAdj] (double theta)
        {
            return integrandAuxForAlphaEqualOne(theta, xAdj) - 1.0;
        },
        -almostPI_2, M_PI_2, theta0);
    }

    std::function<double (double)> integrandPtr = std::bind(&StableRand::integrandForAlphaEqualOne, this, std::placeholders::_1, xAdj);
    double int1 = RandMath::integral(integrandPtr, -M_PI_2, theta0);
    double int2 = RandMath::integral(integrandPtr, theta0, M_PI_2);

    return std::fabs(pdfCoef) * (int1 + int2) / sigma;
}

double StableRand::f(double x) const
{
    /// Check all 'good' cases
    if (alpha == 2)
        return N.f(x);
    if (alpha == 1 && beta == 0)
        return C.f(x);
    if (alpha == .5 && beta == 1)
        return L.f(x);
    if (alpha == .5 && beta == -1)
        return L.f(-x);

    /// Now check the others
    if (alpha == 1)
        return pdfForAlphaEqualOne(x);
    return pdfForCommonAlpha(x);
}

double StableRand::cdfForCommonAlpha(double x) const
{
    x = (x - mu) / sigma; /// Standardize
    
    if (std::fabs(x) < 0.01) /// if we are close to 0 then we do interpolation avoiding dangerous variates
    {
        double y0 = 0.5 - M_1_PI * xi; /// f(0)
        if (std::fabs(x) < MIN_POSITIVE)
            return y0;
        double b = (x > 0) ? 0.011 : -0.011;
        double y1 = cdfForCommonAlpha(mu + sigma * b);
        return RandMath::linearInterpolation(0, b, y0, y1, x); // TODO: find a way to do better interpolation than linear
    }
    
    double xiAdj = xi; /// +- xi
    double xAdj;
    if (x > 0)
    {
        if (alpha < 1 && beta == -1)
            return 1.0;
        xAdj = std::pow(x, alpha_alpham1);
    }
    else
    {
        if (alpha < 1 && beta == 1)
            return 0.0;
        xAdj = std::pow(-x, alpha_alpham1);
        xiAdj = -xi;
    }

    double y = RandMath::integral([this, xAdj, xiAdj] (double theta)
    {
        return std::exp(-integrandAuxForCommonAlpha(theta, xAdj, xiAdj));
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
    return x;
}

double StableRand::F(double x) const
{
    /// Check all 'good' cases
    if (alpha == 2)
        return N.F(x);
    if (alpha == 1 && beta == 0)
        return C.F(x);
    if (alpha == .5 && beta == 1)
        return L.F(x);
    if (alpha == .5 && beta == -1)
        return 1.0 - L.f(-x);

    /// Now check the others
    if (alpha == 1)
        return cdfForAlphaEqualOne(x);
    return cdfForCommonAlpha(x);
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
    double pi_2BetaU = M_PI_2 + beta * U;
    double X = logSigma;
    X -= std::log(M_PI_2 * W * std::cos(U) / pi_2BetaU);
    X *= beta;
    X += pi_2BetaU * std::tan(U);
    X *= M_2_PI;
    return mu + sigma * X;
}

double StableRand::variate() const
{
    /// Check all 'good' cases
    if (alpha == 2)
        return N.variate();
    if (alpha == 1 && beta == 0)
        return C.variate();
    if (alpha == .5 && beta == 1)
        return L.variate();
    if (alpha == .5 && beta == -1)
        return -L.variate();

    /// Now check the others
    if (alpha == 1)
        return variateForAlphaEqualOne();
    return variateForCommonAlpha();
}

void StableRand::sample(QVector<double> &outputData) const
{
    /// Check all 'good' cases
    if (alpha == 2) {
        for (double &var : outputData)
            var = N.variate();
    }
    else if (alpha == 1 && beta == 0) {
        for (double &var : outputData)
            var = C.variate();
    }
    else if (alpha == .5 && beta == 1) {
        for (double &var : outputData)
            var = L.variate();
    }
    else if (alpha == .5 && beta == -1) {
        for (double &var : outputData)
            var = -L.variate();
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
    double x = (alpha == 1) ? -beta * M_2_PI * log(std::fabs(t)) : -zeta;
    if (t > 0)
        x = -x;
    double re = -std::pow(std::fabs(sigma * t), alpha);
    return std::complex<double>(re, re * x + mu * t);
}

std::complex<double> StableRand::CF(double t) const
{
    return std::exp(-psi(t));
}

double StableRand::Mean() const
{
    return (alpha > 1) ? mu : NAN;
}

double StableRand::Variance() const
{
    return (alpha == 2) ? 2 * sigma * sigma : INFINITY;
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
            + toStringWithPrecision(getSigma()) + ", "
            + toStringWithPrecision(getMu()) + ")";
}


std::string LandauRand::name()
{
    return "Landau("
            + toStringWithPrecision(getSigma()) + ", "
            + toStringWithPrecision(getMu()) + ")";
}
