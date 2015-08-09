#include "UniformDiscreteRand.h"


UniformDiscreteRand::UniformDiscreteRand(int minValue, int maxValue)
{
    setBoundaries(minValue, maxValue);
}

void UniformDiscreteRand::setBoundaries(int minValue, int maxValue)
{
    a = minValue;
    b = maxValue;

    if (b < a)
        SWAP(a, b);

    n = b - a + 1;
    delta = n * BasicRandGenerator::maxInv();
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
    return a + std::floor(BasicRandGenerator::getRand() * delta);
}
