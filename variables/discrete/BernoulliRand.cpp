#include "BernoulliRand.h"

BernoulliRand::BernoulliRand(double probability)
{
    setP(probability);
}

void BernoulliRand::setP(double probability)
{
    p = std::min(std::max(probability, 0.0), 1.0);
    generatorEdge = (1 - p) * BasicRandGenerator::max();
}

double BernoulliRand::P(int k) const
{
    return (k == 0) ? 1 - p : ((k == 1) ? p : 0);
}

double BernoulliRand::F(double x) const
{
    return (x < 0) ? 0 : ((x < 1) ? 1 - p : 1);
}

double BernoulliRand::variate()
{
    return (BasicRandGenerator::getRand() < generatorEdge) ? 0 : 1;
}
