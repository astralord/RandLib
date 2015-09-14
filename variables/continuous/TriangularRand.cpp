#include "TriangularRand.h"

TriangularRand::TriangularRand(double lowerLimit, double mode, double upperLimit)
{
    //TODO: add sanity check and merge setters!
    setLowerLimit(lowerLimit);
    setMode(mode);
    setUpperLimit(upperLimit);
}

std::string TriangularRand::name()
{
    return "Triangular("
            + toStringWithPrecision(getLowerLimit()) + ", "
            + toStringWithPrecision(getMode()) + ", "
            + toStringWithPrecision(getUpperLimit()) + ")";
}

void TriangularRand::setConstantsForGenerator()
{
    constForGenerator = (c - a) / (b - a);
    coefGenerator1 = (b - a) * (c - a);
    coefGenerator2 = (b - a) * (b - c);
}

void TriangularRand::setLowerLimit(double lowerLimit)
{
    a = lowerLimit;
    setConstantsForGenerator();
}

void TriangularRand::setMode(double mode)
{
    c = mode;
    setConstantsForGenerator();
}

void TriangularRand::setUpperLimit(double upperLimit)
{
    b = upperLimit;
    setConstantsForGenerator();
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
