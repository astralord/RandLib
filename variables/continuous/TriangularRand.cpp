#include "TriangularRand.h"

TriangularRand::TriangularRand(double lowerLimit, double mode, double upperLimit)
{
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

double TriangularRand::Var() const
{
    return (a * (a - b) + b * (b - c) + c * (c - a)) / 18;
}
