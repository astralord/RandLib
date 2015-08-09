#include "GeometricRand.h"

GeometricRand::GeometricRand(double probability)
{
    setProbability(probability);
}

void GeometricRand::setProbability(double probability)
{
    p = std::min(std::max(probability, MIN_POSITIVE), 1.0);

    /// we use two different generators for two different cases
    /// if p < 0.2 then the tail is too heavy
    /// (probability to be in main body is less then 0.977482)
    /// thus we always use low bound of variate exponential distribution
    /// otherwise we choose table method
    if (p < 0.2)
    {
        W.setRate(-std::log(1 - p));
    }
    else
    {
        table[0] = p;
        double prod = p;
        for (unsigned i = 1; i < tableSize; ++i)
        {
            prod *= (1 - p);
            table[i] = table[i - 1] + prod;
        }
    }
}

double GeometricRand::P(int k) const
{
    return p * std::pow(1 - p, k);
}

double GeometricRand::F(double x) const
{
    return 1 - std::pow(1 - p, std::floor(x) + 1);
}

double GeometricRand::variate() const
{
    if (p < 0.2)
        return variateForSmallP();
    return variateForLargeP();
}

double GeometricRand::variateForSmallP() const
{
    return std::floor(W.variate());
}

double GeometricRand::variateForLargeP() const
{
    double U = UniformRand::standardVariate();
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateForLargeP();
    /// handle the main body
    unsigned x = 0;
    while (U > table[x])
        ++x;
    return x;
}

