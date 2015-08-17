#include "BernoulliRand.h"

BernoulliRand::BernoulliRand(double probability)
{
    setProbability(probability);
}

void BernoulliRand::setName()
{
    nameStr = "Bernoulli(" + toStringWithPrecision(getProbability()) + ")";
}

void BernoulliRand::setProbability(double probability)
{
    p = std::min(std::max(probability, 0.0), 1.0);
    generatorEdge = (1 - p) * BasicRandGenerator::maxValue();

    setName();
}

double BernoulliRand::P(int k) const
{
    return (k == 0) ? (1 - p) : ((k == 1) ? p : 0);
}

double BernoulliRand::F(double x) const
{
    return (x < 0) ? 0 : ((x < 1) ? 1 - p : 1);
}

double BernoulliRand::variate() const
{
    return (BasicRandGenerator::variate() < generatorEdge) ? 0 : 1;
}
