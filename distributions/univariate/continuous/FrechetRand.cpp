#include "FrechetRand.h"
#include "ExponentialRand.h"

FrechetRand::FrechetRand(double shape, double scale, double location)
{
    setParameters(shape, scale, location);
}

std::string FrechetRand::name() const
{
    return "Frechet(" + toStringWithPrecision(getShape()) + ", "
                      + toStringWithPrecision(getScale()) + ", "
                      + toStringWithPrecision(getLocation()) + ")";
}

void FrechetRand::setParameters(double shape, double scale, double location)
{
    alpha = shape;
    if (alpha <= 0)
        alpha = 1.0;
    alphaInv = 1.0 / alpha;

    s = scale;
    if (s <= 0)
        s = 1.0;

    m = location;
}

double FrechetRand::f(double x) const
{
    if (x <= m)
        return 0.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    double y = std::exp(-xPow);
    y /= xAdj;
    y *= xPow;
    y *= alpha;
    return y / s;
}

double FrechetRand::F(double x) const
{
    if (x <= m)
        return 0.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    return std::exp(-xPow);
}

double FrechetRand::variate() const
{
    return m + s / std::pow(ExponentialRand::standardVariate(), alphaInv);
}

double FrechetRand::Mean() const
{
    if (alpha <= 1.0)
        return INFINITY;
    return m + s * std::tgamma(1.0 - alphaInv);
}

double FrechetRand::Variance() const
{
    if (alpha <= 2.0)
        return INFINITY;
    double var = std::tgamma(1.0 - alphaInv);
    var *= var;
    var = std::tgamma(1.0 - 2 * alphaInv) - var;
    return s * s * var;
}

double FrechetRand::quantileImpl(double p) const
{
    double x = -std::log(p);
    x = s / std::pow(x, alphaInv);
    return x + m;
}

double FrechetRand::quantileImpl1m(double p) const
{
    double x = -std::log1p(-p);
    x = s / std::pow(x, alphaInv);
    return x + m;
}

double FrechetRand::Median() const
{
    return m + s / std::pow(M_LN2, alphaInv);
}

double FrechetRand::Mode() const
{
    double y = alpha / (1.0 + alpha);
    y = std::pow(y, alphaInv);
    return m + s * y;
}

double FrechetRand::Skewness() const
{
    if (alpha <= 3.0)
        return INFINITY;
    double x = std::tgamma(1.0 - alphaInv);
    double y = std::tgamma(1.0 - 2.0 * alphaInv);
    double z = std::tgamma(1.0 - 3.0 * alphaInv);
    double numerator = 2 * x * x - 3 * y;
    numerator *= x;
    numerator += z;
    double denominator = y - x * x;
    denominator = std::pow(denominator, 1.5);
    return numerator / denominator;
}

double FrechetRand::ExcessKurtosis() const
{
    if (alpha <= 4.0)
        return INFINITY;
    double x = std::tgamma(1.0 - alphaInv);
    double y = std::tgamma(1.0 - 2.0 * alphaInv);
    double z = std::tgamma(1.0 - 3.0 * alphaInv);
    double w = std::tgamma(1.0 - 4.0 * alphaInv);
    double numerator = w - 4 * z * x + 3 * y * y;
    double denominator = y - x * x;
    denominator *= denominator;
    return numerator / denominator;
}

double FrechetRand::Entropy() const
{
    return 1.0 + M_EULER * (1.0 + alphaInv) + std::log(s / alpha);
}
