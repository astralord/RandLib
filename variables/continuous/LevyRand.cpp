#include "LevyRand.h"

LevyRand::LevyRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LevyRand::name()
{
    return "Levy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LevyRand::setLocation(double location)
{
    mu = location;
}

void LevyRand::setScale(double scale)
{
    c = scale;
    if (c <= 0)
        c = 1.0;
    sqrtc_2pi = std::sqrt(c);
    X.setSigma(1.0 / sqrtc_2pi); /// X ~ N(0, c ^ -0.5)
    sqrtc_2pi *= M_1_SQRT2PI;
}

double LevyRand::f(double x) const
{
    if (x <= mu)
        return 0;
    double xInv = 1.0 / (x - mu);
    double y = -0.5 * c * xInv;
    y = std::exp(y);
    y *= xInv;
    y *= std::sqrt(xInv);
    return sqrtc_2pi * y;
}

double LevyRand::F(double x) const
{
    if (x <= mu)
        return 0;
    double y = x - mu;
    y = 0.5 * c / y;
    y = std::sqrt(y);
    return std::erfc(y);
}

double LevyRand::variate() const
{
    double rv = X.variate();
    rv *= rv;
    rv = 1.0 / rv;
    return mu + rv;
}

double LevyRand::Mean() const
{
    return INFINITY;
}

double LevyRand::Variance() const
{
    return INFINITY;
}

std::complex<double> LevyRand::CF(double t) const
{
    std::complex<double> y(0.0, -2 * c * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

double LevyRand::Quantile(double p) const
{
    double x = X.Quantile(0.5 * p);
    return mu + x * x;
}

double LevyRand::Mode() const
{
    return c / 3.0 + mu;
}

double LevyRand::Skewness() const
{
    return NAN;
}

double LevyRand::ExcessKurtosis() const
{
    return NAN;
}
