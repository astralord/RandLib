#include "StableRand.h"

StableRand::StableRand(double exponent, double skewness, double scale, double location) :
    U(-M_PI_2, M_PI_2),
    E(1)
{
    setAlphaAndBeta(exponent, skewness);
    setSigma(scale);
    setMu(location);
}

void StableRand::setAlphaAndBeta(double exponent, double skewness)
{
    alpha = std::min(exponent, 2.0);
    alpha = std::max(alpha, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    if (alpha != 1 && qFabs(beta) > MIN_POSITIVE)
    {
        B = beta * std::tan(M_PI_2 * alpha);
        S = std::pow(1 + B * B, 0.5 * alphaInv);
        B = std::atan(B);
    }
    else
    {
        B = 0;
        S = 1;
    }
}

void StableRand::setSigma(double scale)
{
    sigma = std::max(scale, MIN_POSITIVE);
    logSigma = std::log(sigma);
}

void StableRand::setMu(double location)
{
    mu = location;
}

double StableRand::pdf(double x)
{
    return x;
}

double StableRand::cdf(double x)
{
    return x;
}

double StableRand::value()
{
    // alpha = 2 -> X ~ Normal(mu, 2 sigma^2)
    // alpha = 1 && beta = 0 -> X ~ Cauchy(mu, sigma)
    // alpha = 0.5 && beta = +-1 -> X ~ +-Levy(mu, sigma)
    // Verify that it is one of three known analytically expressed distributions

    double V = U.value();
    double W = E.value();
    double rv;

    if (alpha != 1) /// the most common case
    {
        double alphaV = alpha * V;
        rv = S * qFastSin(alphaV + B);
        W /= qFastCos(V - alphaV - B);
        rv *= W;
        rv *= std::pow(W * qFastCos(V), -alphaInv);
    }
    else
    {
        double pi_2_beta_V = M_PI_2 + beta * V;
        rv = pi_2_beta_V * qTan(V);
        rv -= std::log(W * qFastCos(V) / pi_2_beta_V);
        rv += logSigma;
        rv *= M_2_PI * beta;
    }

    return mu + sigma * rv;
}
