#include "LimitingDistribution.h"

LimitingDistribution::LimitingDistribution(double exponent, double skewness, double scale, double location)
{
    setParameters(exponent, skewness);
    setLocation(location);
    setScale(scale);
}

void LimitingDistribution::setParameters(double exponent, double skewness)
{
    alpha = std::min(exponent, 2.0);
    if (alpha <= 0)
        alpha = 2.0;
    alphaInv = 1.0 / alpha;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    /// Should be cautious, known distributions in priority
    if (RandMath::areClose(alpha, 2))
        alpha = 2;
    else if (RandMath::areClose(alpha, 1))
        alpha = 1;
    else if (RandMath::areClose(alpha, 0.5))
        alpha = 0.5;

    if (alpha != 1 /*&& alpha != 2*/ ) /// Common case
    {
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = std::pow(1 + B * B, .5 * alphaInv);
        B = std::atan(B);
    }
}

void LimitingDistribution::setLocation(double location)
{
    mu = location;
}

void LimitingDistribution::setScale(double scale)
{
    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    if (alpha == 1)
        logSigma = std::log(sigma);
}

double LimitingDistribution::Mean() const
{
    if (alpha > 1)
        return mu;
    return (alpha == 1.0) ? NAN : INFINITY;
}

double LimitingDistribution::Variance() const
{
    return (alpha == 2) ? 2 * sigma * sigma : INFINITY;
}

double LimitingDistribution::Median() const
{
    return (beta == 0) ? mu : ContinuousDistribution::Median();
}

double LimitingDistribution::Mode() const
{
    return (beta == 0) ? mu : ContinuousDistribution::Mode();
}

std::complex<double> LimitingDistribution::psi(double t) const
{
    double x = (alpha == 1) ? beta * M_2_PI * log(std::fabs(t)) : zeta;
    if (t > 0)
        x = -x;
    double re = std::pow(std::fabs(sigma * t), alpha);
    return std::complex<double>(re, re * x + mu * t);
}
