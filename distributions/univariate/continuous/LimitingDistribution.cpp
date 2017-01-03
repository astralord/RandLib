#include "LimitingDistribution.h"

LimitingDistribution::LimitingDistribution(double exponent, double skewness, double scale, double location)
{
    SetParameters(exponent, skewness);
    SetLocation(location);
    SetScale(scale);
}

void LimitingDistribution::SetParameters(double exponent, double skewness)
{
    alpha = std::min(std::max(exponent, 0.1), 2.0);
    alphaInv = 1.0 / alpha;

    beta = std::min(skewness, 1.0);
    beta = std::max(beta, -1.0);

    lambda = std::pow(sigma, alpha);

    if (alpha != 1 && alpha != 2) /// Common case
    {   // TODO: don't do also for Levy
        B = beta * std::tan(M_PI_2 * alpha);
        zeta = -B;
        S = 0.5 * alphaInv * std::log1p(zeta * zeta);
        B = std::atan(B);
    }
}

void LimitingDistribution::SetLocation(double location)
{
    mu = location;
}

void LimitingDistribution::SetScale(double scale)
{
    sigma = scale > 0 ? scale : 1.0;
    if (alpha == 1)
        logSigma = std::log(sigma);
    lambda = std::pow(sigma, alpha);
}

double LimitingDistribution::Mean() const
{
    if (alpha > 1)
        return mu;
    if (beta != 1)
        return NAN;
    return (beta > 0) ? INFINITY : -INFINITY;
}

std::complex<double> LimitingDistribution::psi(double t) const
{
    double fabsT = std::fabs(t);
    double x = (alpha == 1) ? beta * M_2_PI * log(fabsT) : zeta;
    if (t < 0)
        x = -x;
    double re = std::pow(sigma * fabsT, alpha);
    return std::complex<double>(re, re * x - mu * t);
}
