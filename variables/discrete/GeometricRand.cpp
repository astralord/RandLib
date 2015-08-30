#include "GeometricRand.h"

GeometricRand::GeometricRand(double probability)
{
    setProbability(probability);
}

std::string GeometricRand::name()
{
    return "Geometric(" + toStringWithPrecision(getProbability()) + ")";
}

void GeometricRand::setProbability(double probability)
{
    p = std::min(std::max(probability, MIN_POSITIVE), 1.0);
    q = 1.0 - p;

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
        for (int i = 1; i < tableSize; ++i)
        {
            prod *= q;
            table[i] = table[i - 1] + prod;
        }
    }
}

double GeometricRand::P(int k) const
{
    return p * std::pow(q, k);
}

double GeometricRand::F(double x) const
{
    return 1 - std::pow(q, std::floor(x) + 1);
}

double GeometricRand::variate() const
{
    if (p < 0.2)
        return variateForSmallP();
    return variateForLargeP();
}

void GeometricRand::sample(QVector<double> &outputData)
{
    if (p < 0.2) {
        for (double &var : outputData)
            var = variateForSmallP();
    }
    else {
        for (double &var : outputData)
            var = variateForLargeP();
    }
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
    int x = 0;
    while (U > table[x])
        ++x;
    return x;
}

