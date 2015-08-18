#include "UniformRand.h"

UniformRand::UniformRand(double minValue, double maxValue)
{
    setBoundaries(minValue, maxValue);
}

void UniformRand::setName()
{
    nameStr = "Uniform(" + toStringWithPrecision(getMinValue()) + ", " + toStringWithPrecision(getMaxValue()) + ")";
}

void UniformRand::setBoundaries(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    /// Sanity check
    if (b < a)
        SWAP(a, b);
    if (b - a < MIN_POSITIVE)
        b = a + MIN_POSITIVE;

    c = 1.0 / (b - a);

    setName();
}

double UniformRand::variate() const
{
    return a + standardVariate() * (b - a);
}

double UniformRand::variate(double minValue, double maxValue)
{
    return minValue + standardVariate() * (maxValue - minValue);
}

double UniformRand::standardVariate()
{
    double x;
    unsigned long long a = BasicRandGenerator::variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bits
    *((unsigned long long *)&x) = a; /// Make a double from bits
    return x - 1.0;
}

double UniformRand::f(double x) const
{
    return (x >= a && x <= b) ? c : 0;
}

double UniformRand::F(double x) const
{
    if (x < a)
        return 0;
    if (x > b)
        return 1;
    return c * (x - a);
}

bool UniformRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;
    double maxVar = sample.at(0), minVar = maxVar;
    for (double var : sample) {
        maxVar = std::max(var, maxVar);
        minVar = std::min(var, minVar);
    }

    setBoundaries(minVar, maxVar);
    return true;
}
