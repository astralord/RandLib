#include "UniformDiscreteRand.h"


UniformDiscreteRand::UniformDiscreteRand(int minValue, int maxValue)
{
    setBoundaries(minValue, maxValue);
}

std::string UniformDiscreteRand::name()
{
    return "Uniform(" + toStringWithPrecision(getMinValue()) + ", " + toStringWithPrecision(getMaxValue()) + ")";
}

void UniformDiscreteRand::setBoundaries(int minValue, int maxValue)
{
    a = minValue;
    b = maxValue;

    if (b < a)
        SWAP(a, b);

    n = b - a + 1;
}

double UniformDiscreteRand::P(int k) const
{
    if (k < a || k > b)
        return 0;
    return 1.0 / n;
}

double UniformDiscreteRand::F(double x) const
{
    return (std::floor(x) - a + 1) / n;
}

double UniformDiscreteRand::variate() const
{
    return a + RandGenerator::variate() % n;
}

double UniformDiscreteRand::E() const
{
    return .5 * (b + a);
}

double UniformDiscreteRand::Var() const
{
    return (n * n - 1) / 12;
}

double UniformDiscreteRand::quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    double x = p * n + a - 1;
    double xFloor = std::floor(x);
    if (RandMath::areEqual(x, xFloor))
        return xFloor;
    return NAN;
}

double UniformDiscreteRand::Median() const
{
    return .5 * (b + a);
}

double UniformDiscreteRand::Mode() const
{
    /// any from {a, ..., b}
    return variate();
}

double UniformDiscreteRand::Skewness() const
{
    return 0.0;
}

double UniformDiscreteRand::ExcessKurtosis() const
{
    double n2 = n * n;
    return -1.2 * (n2 + 1) / (n2 - 1);
}
