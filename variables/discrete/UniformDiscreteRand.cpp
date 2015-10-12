#include "UniformDiscreteRand.h"


UniformDiscreteRand::UniformDiscreteRand(int minValue, int maxValue)
{
    setBoundaries(minValue, maxValue);
}

std::string UniformDiscreteRand::name()
{
    return "Uniform Discrete(" + toStringWithPrecision(getMinValue()) + ", " + toStringWithPrecision(getMaxValue()) + ")";
}

void UniformDiscreteRand::setBoundaries(int minValue, int maxValue)
{
    a = minValue;
    b = maxValue;

    if (b < a)
        SWAP_INTEGER(a, b);

    n = b - a + 1;
    nInv = 1.0 / n;
}

double UniformDiscreteRand::P(int k) const
{
    if (k < a || k > b)
        return 0;
    return nInv;
}

double UniformDiscreteRand::F(double x) const
{
    if (x < a)
        return 0.0;
    if (x > b)
        return 1.0;
    return (std::floor(x) - a + 1) * nInv;
}

double UniformDiscreteRand::variate() const
{
    return a + RandGenerator::variate() % n;
}

double UniformDiscreteRand::Mean() const
{
    return .5 * (b + a);
}

double UniformDiscreteRand::Variance() const
{
    return (n * n - 1) / 12;
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
    double nSq = n * n;
    return -1.2 * (nSq + 1) / (nSq - 1);
}

std::complex<double> UniformDiscreteRand::CF(double t) const
{
    std::complex<double> y(0.0, a * t);
    y = std::exp(y);
    std::complex<double> z(0.0, b * t + b);
    z = std::exp(z);
    y -= z;
    z = std::complex<double>(0.0, t);
    z = 1.0 - std::exp(z);
    return nInv * y / z;
}
