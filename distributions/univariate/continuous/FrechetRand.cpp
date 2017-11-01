#include "FrechetRand.h"
#include "ExponentialRand.h"

FrechetRand::FrechetRand(double shape, double scale, double location)
{
    SetParameters(shape, scale, location);
}

String FrechetRand::Name() const
{
    return "Frechet(" + toStringWithPrecision(GetShape()) + ", "
                      + toStringWithPrecision(GetScale()) + ", "
                      + toStringWithPrecision(GetLocation()) + ")";
}

void FrechetRand::SetParameters(double shape, double scale, double location)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Frechet distribution: shape should be positive");
    if (scale <= 0.0)
        throw std::invalid_argument("Frechet distribution: scale should be positive");
    alpha = shape;
    alphaInv = 1.0 / alpha;
    s = scale;
    pdfCoef = std::log(alpha / s);
    m = location;
}

double FrechetRand::f(const double & x) const
{
    return (x <= m) ? 0.0 : std::exp(logf(x));
}

double FrechetRand::logf(const double & x) const
{
    if (x <= m)
        return -INFINITY;
    double xAdj = (x - m) / s;
    double logxAdj = std::log(xAdj);
    double a = alpha * logxAdj;
    double expA = std::exp(-a);
    return pdfCoef - a - expA - logxAdj;
}

double FrechetRand::F(const double & x) const
{
    if (x <= m)
        return 0.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    return std::exp(-xPow);
}

double FrechetRand::S(const double & x) const
{
    if (x <= m)
        return 1.0;
    double xAdj = (x - m) / s;
    double xPow = std::pow(xAdj, -alpha);
    return -std::expm1(-xPow);
}

double FrechetRand::Variate() const
{
    return m + s / std::pow(ExponentialRand::StandardVariate(), alphaInv);
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
    double y = -std::log(p);
    y = s / std::pow(y, alphaInv);
    return y + m;
}

double FrechetRand::quantileImpl1m(double p) const
{
    double y = -std::log1p(-p);
    y = s / std::pow(y, alphaInv);
    return y + m;
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
    return numerator / denominator - 6.0;
}

double FrechetRand::Entropy() const
{
    return 1.0 + M_EULER * (1.0 + alphaInv) + std::log(s / alpha);
}
