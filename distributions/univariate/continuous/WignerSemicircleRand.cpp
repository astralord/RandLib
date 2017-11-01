#include "WignerSemicircleRand.h"

WignerSemicircleRand::WignerSemicircleRand(double radius)
{
    SetRadius(radius);
}

String WignerSemicircleRand::Name() const
{
    return "Wigner Semicircle(" + toStringWithPrecision(GetRadius()) + ")";
}

void WignerSemicircleRand::SetRadius(double radius)
{
    if (radius <= 0.0)
        throw std::invalid_argument("Wigner-Semicircle distribution: radius should be positive");
    R = radius;
    RSq = R * R;
    logRSq = std::log(RSq);
}

double WignerSemicircleRand::f(const double & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
        return 0.0;
    double y = RSq - xSq;
    y = std::sqrt(y);
    y *= M_1_PI / RSq;
    return 2 * y;
}

double WignerSemicircleRand::logf(const double & x) const
{
    double xSq = x * x;
    if (xSq >= RSq)
        return -INFINITY;
    return M_LN2 + 0.5 * std::log(RSq - xSq) - M_LNPI - logRSq;
}

double WignerSemicircleRand::F(const double & x) const
{
    if (x <= -R)
        return 0.0;
    if (x >= R)
        return 1.0;
    double y = RSq - x * x;
    y = x * std::sqrt(y) / RSq;
    double z = std::asin(x / R);
    return 0.5 + (y + z) / M_PI;
}

double WignerSemicircleRand::Variate() const
{
    double x = X.Variate();
    x += x - 1;
    return R * x;
}

double WignerSemicircleRand::Mean() const
{
    return 0.0;
}

double WignerSemicircleRand::Variance() const
{
    return 0.25 * RSq;
}

double WignerSemicircleRand::Median() const
{
    return 0.0;
}

double WignerSemicircleRand::Mode() const
{
    return 0.0;
}

double WignerSemicircleRand::Skewness() const
{
    return 0.0;
}

double WignerSemicircleRand::ExcessKurtosis() const
{
    return -1.0;
}

double WignerSemicircleRand::Entropy() const
{
    return M_LNPI + 0.5 * logRSq - 0.5;
}

