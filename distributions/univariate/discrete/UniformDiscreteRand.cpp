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
        b = a + 1;

    n = b - a + 1;
    nInv = 1.0 / n;
}

double UniformDiscreteRand::P(int k) const
{
    if (k < a || k > b)
        return 0;
    return nInv;
}

double UniformDiscreteRand::F(int k) const
{
    if (k < a)
        return 0.0;
    if (k > b)
        return 1.0;
    return (k - a + 1) * nInv;
}

int UniformDiscreteRand::variate() const
{
    double x = (RandGenerator::variate() % n);
    return a + x;
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

int UniformDiscreteRand::Mode() const
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
    double res = 2.0 / (n * n - 1);
    ++res;
    return -1.2 * res;
}

std::complex<double> UniformDiscreteRand::CF(double t) const
{
    if (t == 0)
        return std::complex<double>(1, 0);
    double at = a * t;
    double bp1t = (b + 1) * t;
    double reNum = std::cos(at) - std::cos(bp1t);
    double imNum = std::sin(at) - std::sin(bp1t);
    std::complex<double> numen(reNum, imNum);
    std::complex<double>denom(1 - std::cos(t), -std::sin(t));
    return nInv * numen / denom;
}
