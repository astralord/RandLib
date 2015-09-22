#include "UniformRand.h"

UniformRand::UniformRand(double minValue, double maxValue)
{
    setBoundaries(minValue, maxValue);
}

std::string UniformRand::name()
{
    return "Uniform(" + toStringWithPrecision(getMinValue()) + ", " + toStringWithPrecision(getMaxValue()) + ")";
}

void UniformRand::setBoundaries(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    /// Sanity check
    if (b < a)
        SWAP(a, b);
    if (RandMath::areEqual(b, a))
        b = a + 1.0;

    c = 1.0 / (b - a);
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
    unsigned long long a = RandGenerator::variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bits
    *((unsigned long long *)&x) = a; /// Make a double from bits
    return x - 1.0;
}

double UniformRand::Mean() const
{
    return .5 * (b + a);
}

double UniformRand::Variance() const
{
    return (b - a) * (b - a) / 12;
}

std::complex<double> UniformRand::CF(double t) const
{
    std::complex<double> x(0, -t * a), y(0, -t * b);
    std::complex<double> numerator = std::exp(x) - std::exp(y);
    std::complex<double> denominator(0, t * (b - a));
    return numerator / denominator;
}

double UniformRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    return a + (b - a) * p;
}

double UniformRand::Median() const
{
    return .5 * (b + a);
}

double UniformRand::Mode() const
{
    /// this can be any value in [a, b]
    return variate();
}

double UniformRand::Skewness() const
{
    return 0.0;
}

double UniformRand::ExcessKurtosis() const
{
    return -1.2;
}

double UniformRand::Entropy() const
{
    return (b == a) ? -INFINITY : std::log(b - a);
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
