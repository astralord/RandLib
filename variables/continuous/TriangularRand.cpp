#include "TriangularRand.h"
#include "UniformRand.h"

TriangularRand::TriangularRand(double lowerLimit, double mode, double upperLimit)
{
    setParameters(lowerLimit, mode, upperLimit);
}

std::string TriangularRand::name()
{
    return "Triangular("
            + toStringWithPrecision(getLowerLimit()) + ", "
            + toStringWithPrecision(getMode()) + ", "
            + toStringWithPrecision(getUpperLimit()) + ")";
}

void TriangularRand::setParameters(double lowerLimit, double mode, double upperLimit)
{
    a = lowerLimit;
    if (mode > a)
        c = mode;
    else
        c = a + 1.0;
    if (upperLimit > c)
        b = upperLimit;
    else
        b = c + 1.0;
    setConstantsForGenerator();
}

void TriangularRand::setConstantsForGenerator()
{
    constForGenerator = (c - a) / (b - a);
    coefGenerator1 = (b - a) * (c - a);
    coefGenerator2 = (b - a) * (b - c);
}

double TriangularRand::f(double x) const
{
    if (x <= a)
        return 0;
    if (x < c)
        return 2.0 * (x - a) / ((b - a) * (c - a));
    if (x == c)
        return 2.0 / (b - a);
    if (x < b)
        return 2.0 * (b - x) / ((b - a) * (b - c));
    return 0;
}

double TriangularRand::F(double x) const
{
    if (x <= a)
        return 0.0;
    if (x <= c)
        return (x - a) * (x - a) / ((b - a) * (c - a));
    if (x < b)
        return 1.0 - (b - x) * (b - x) / ((b - a) * (b - c));
    return 1.0;
}

double TriangularRand::variate() const
{
    double u = UniformRand::standardVariate();
    if (u < constForGenerator)
        return a + std::sqrt(u * coefGenerator1);
    return b - std::sqrt((1 - u) * coefGenerator2);
}

double TriangularRand::Mean() const
{
    return (a + b + c) / 3.0;
}

double TriangularRand::Variance() const
{
    return (a * (a - b) + b * (b - c) + c * (c - a)) / 18.0;
}

std::complex<double> TriangularRand::CF(double t) const
{
    std::complex<double> x(0, a * t);
    std::complex<double> y(0, c * t);
    std::complex<double> z(0, b * t);
    x = (b - c) * std::exp(x);
    y = (b - a) * std::exp(y);
    z = (c - a) * std::exp(z);
    std::complex<double> numerator = x - y + z;
    double denominator = (b - a);
    denominator *= (c - a);
    denominator *= (b - c);
    denominator *= t * t;
    std::complex<double> frac = -numerator / denominator;
    return frac + frac;
}

double TriangularRand::Median() const
{
    if (c + c > a + b)
        return a + std::sqrt(0.5 * (b - a) * (c - a));
    return b - std::sqrt(0.5 * (b - a) * (b - c));
}

double TriangularRand::Mode() const
{
    return c;
}

double TriangularRand::Skewness() const
{
    double numerator = M_SQRT2;
    numerator *= (a + b - c - c);
    numerator *= (a + a - b - c);
    numerator *= (a - b - b + c);
    double denominator = a * (a - b);
    denominator += b * (b - c);
    denominator += c * (c - a);
    denominator *= std::sqrt(denominator);
    return 0.2 * numerator / denominator;
}

double TriangularRand::ExcessKurtosis() const
{
    return -0.6;
}
