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
}

void LimitingDistribution::SetLocation(double location)
{
    mu = location;
}

void LimitingDistribution::SetScale(double scale)
{
    sigma = scale > 0 ? scale : 1.0;
    if (alpha == 1)
        logsigmaPi_2 = std::log(M_PI_2 * sigma);
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
    double x = (alpha == 1) ? M_2_PI * log(fabsT) : std::tan(M_PI_2 * alpha); // second term can be replaced on zeta / beta
    x *= beta;
    if (t < 0)
        x = -x;
    double re = std::pow(sigma * fabsT, alpha);
    return std::complex<double>(re, re * x - mu * t);
}
