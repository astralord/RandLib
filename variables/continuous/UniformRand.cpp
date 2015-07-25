#include "UniformRand.h"

UniformRand::UniformRand(double minValue, double maxValue)
{
    setBoundaries(minValue, maxValue);
}

void UniformRand::setBoundaries(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    /// Sanity check
    if (b < a)
        swapBoundaries();
    if (b - a < MIN_POSITIVE)
        b = a + MIN_POSITIVE;

    c = 1.0 / (b - a);
    delta = (b - a) * B.maxInv();
}

void UniformRand::swapBoundaries()
{
    a += b;
    b -= a;
    a += b;
    b = -b;
}

double UniformRand::variate()
{
    double rv = B.getRand() * delta; /// ~ U(0, 1)
    return rv + a; /// ~ U(a, b)
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
