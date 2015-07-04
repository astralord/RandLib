#include "StableRand.h"
#include <QDebug>

StableRand::StableRand(double exponent, double skewness, double scale, double location) :    
    U(-M_PI_2, M_PI_2), N(0, M_SQRT2)
{
    setParameters(exponent, skewness, scale, location);
}

void StableRand::setParameters(double exponent, double skewness, double scale, double location)
{
    alpha = std::min(exponent, 2.0);
    alpha = std::max(alpha, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    alpha_alpham1 = alpha / (alpha - 1.0);
    alpham1Inv = alpha_alpham1 - 1.0;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    sigma = std::max(scale, MIN_POSITIVE);
    mu = location;

    /// Should be cautious, known distributions in priority
    if (std::fabs(alpha - 2) < MIN_POSITIVE)
        alpha = 2;
    else if (std::fabs(alpha - 1) < MIN_POSITIVE)
        alpha = 1;
    else if (std::fabs(alpha - .5) < MIN_POSITIVE)
        alpha = .5;

    if (std::fabs(beta) < MIN_POSITIVE)
        beta = 0;
    else if (std::fabs(beta - 1) < MIN_POSITIVE)
        beta = 1;
    else if (std::fabs(beta + 1) < MIN_POSITIVE)
        beta = -1;

    if (alpha == 2) /// X ~ Normal(mu, 2sigma^2)
    {
        N.setMean(mu);
        N.setSigma(sigma * M_SQRT2);
        valuePtr = std::bind(&StableRand::valueNormal, this);
        pdfPtr = std::bind(&StableRand::pdfNormal, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfNormal, this, std::placeholders::_1);
    }
    else if (alpha == 1)
    {
        if (beta == 0) /// X ~ Cauchy(mu, sigma)
        {
            C.setLocation(mu);
            C.setScale(sigma);
            valuePtr = std::bind(&StableRand::valueCauchy, this);
            pdfPtr = std::bind(&StableRand::pdfCauchy, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfCauchy, this, std::placeholders::_1);
        }
        else /// just alpha == 1
        {
            logSigma = std::log(sigma);
            pdfCoef = 0.5 / beta;

            valuePtr = std::bind(&StableRand::valueForAlphaEqualOne, this);
            pdfPtr = std::bind(&StableRand::pdfForAlphaEqualOne, this, std::placeholders::_1);
            cdfPtr = std::bind(&StableRand::cdfForAlphaEqualOne, this, std::placeholders::_1);
        }
    }
    else if (alpha == .5 && beta == 1) /// X ~ Levy(mu, sigma)
    {
        L.setLocation(mu);
        L.setScale(sigma);
        valuePtr = std::bind(&StableRand::valueLevy, this);
        pdfPtr = std::bind(&StableRand::pdfLevy, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfLevy, this, std::placeholders::_1);
    }
    else if (alpha == .5 && beta == -1) /// -X ~ Levy(mu, sigma)
    {
        L.setLocation(mu);
        L.setScale(sigma);
        valuePtr = std::bind(&StableRand::valueLevyNegative, this);
        pdfPtr = std::bind(&StableRand::pdfLevyNegative, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfLevyNegative, this, std::placeholders::_1);
    }
    else /// Common case: alpha != 1
    {
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = std::pow(1 + B * B, .5 * alphaInv);
        B = std::atan(B);
        pdfCoef = M_1_PI * alpha / (std::fabs(1 - alpha) * sigma);
        xi = alphaInv * B;
        integrandCoef = std::pow(qFastCos(B), alpham1Inv);

        valuePtr = std::bind(&StableRand::valueForCommonAlpha, this);
        pdfPtr = std::bind(&StableRand::pdfForCommonAlpha, this, std::placeholders::_1);
        cdfPtr = std::bind(&StableRand::cdfForCommonAlpha, this, std::placeholders::_1);
    }
}

double StableRand::value()
{
    /// Get standard value
    double rv = 0;
    int iter = 0;
    do {
        rv = valuePtr();
        ++iter;
    } while ((std::isnan(rv) || std::isinf(rv)) &&
             iter < 10); /// if we got nan 10 times - we have a problem, get out
    return rv;
}

double StableRand::valueForCommonAlpha()
{
    double V = U.value();
    double W = E.value();
    double alphaVB = alpha * V + B;
    double rv = S * qFastSin(alphaVB); /// S * sin(alpha * V + B)
    double W_adj = W / qFastCos(V - alphaVB);
    rv *= W_adj; /// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B)
    rv *= std::pow(W_adj * qFastCos(V), -alphaInv);/// S * sin(alpha * V + B) * W / cos((1 - alpha) * V - B) /
                                                   /// ((W * cos(V) / cos((1 - alpha) * V - B)) ^ (1 / alpha))
    return mu + sigma * rv;
}

double StableRand::valueForAlphaEqualOne()
{
    double V = U.value();
    double W = E.value();
    double pi_2BetaV = M_PI_2 + beta * V;

    double rv = logSigma;
    rv -= std::log(M_PI_2 * W * qFastCos(V) / pi_2BetaV);
    rv *= beta;
    rv += pi_2BetaV * std::tan(V);
    rv *= M_2_PI;
    return mu + sigma * rv;
}

double StableRand::f(double x) const
{
    return pdfPtr(x);
}

double StableRand::pdfForCommonAlpha(double x)
{
    x = (x - mu) / sigma; /// Standardize

    if (std::fabs(x) < MIN_POSITIVE)
    {
        double numen = std::tgamma(1 + alphaInv);
        numen *= qFastCos(xi);
        double denom = std::pow(1 + xi * xi, .5 * alphaInv);
        return M_1_PI * numen / denom;
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
    integrandPtr = std::bind(&StableRand::integrandForCommonAlpha, this, std::placeholders::_1, xAdj, xiAdj);
    return pdfCoef * RandMath::integral(integrandPtr, -xiAdj, M_PI_2) / x;
}

double StableRand::pdfForAlphaEqualOne(double x)
{
    x = (x - mu) / sigma - M_2_PI * beta * logSigma; /// Standardize

    double xAdj = std::exp(-M_PI * x * pdfCoef);
    integrandPtr = std::bind(&StableRand::integrandForAlphaEqualOne, this, std::placeholders::_1, xAdj);
    return std::fabs(pdfCoef) * RandMath::integral(integrandPtr, -M_PI_2, M_PI_2) / sigma;
}

double StableRand::integrandForAlphaEqualOne(double theta, double x_adj) const
{
    double cosTheta = qFastCos(theta);
    /// if theta ~ +-pi / 2
    if (std::fabs(cosTheta) < MIN_POSITIVE)
        return 0;
    double thetaAdj = (M_PI_2 + beta * theta) / cosTheta;
    double u = M_2_PI * thetaAdj;
    u *= std::exp(thetaAdj * qFastSin(theta) / beta);
    u *= x_adj;
    return u * std::exp(-u);
}

double StableRand::integrandForCommonAlpha(double theta, double x_adj, double xi_adj) const
{
    double thetaAdj = alpha * (theta + xi_adj);
    double sinThetaAdj = qFastSin(thetaAdj);
    /// if theta ~ 0
    if (std::fabs(sinThetaAdj) < MIN_POSITIVE)
        return 0;
    double cosTheta = qFastCos(theta);
    double y = cosTheta / sinThetaAdj;
    /// if theta ~ pi / 2
    if (std::fabs(y) < MIN_POSITIVE)
        return 0;
    y = std::pow(y, alpha_alpham1);
    y *= qFastCos(thetaAdj - theta);
    y /= cosTheta;
    double u = integrandCoef * x_adj * y;
    return u * std::exp(-u);
}

double StableRand::F(double x) const
{
    return cdfPtr(x);
}
