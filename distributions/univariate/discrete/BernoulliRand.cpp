#include "BernoulliRand.h"
#include "../continuous/UniformRand.h"
#include "../BasicRandGenerator.h"

BernoulliRand::BernoulliRand(double probability) : BinomialRand(1, probability)
{
    boundary = q * RandGenerator::maxValue();
}

std::string BernoulliRand::name()
{
    return "Bernoulli(" + toStringWithPrecision(getProbability()) + ")";
}

void BernoulliRand::setProbability(double probability)
{
    setParameters(1, probability);
    boundary = q * RandGenerator::maxValue();
}

double BernoulliRand::P(int k) const
{
    return (k == 0) ? q : ((k == 1) ? p : 0);
}

double BernoulliRand::F(int k) const
{
    return (k < 0) ? 0 : ((k < 1) ? q : 1);
}

int BernoulliRand::variate() const
{
    return RandGenerator::variate() > boundary;
}

int BernoulliRand::variate(double p)
{
    return UniformRand::standardVariate() <= p;
}

int BernoulliRand::standardVariate()
{
    static const size_t maxDecimals = RandGenerator::maxDecimals();
    static size_t decimals = 1;
    static unsigned long long X = 0;
    if (decimals == 1)
    {
        /// refresh
        decimals = maxDecimals;
        X = RandGenerator::variate();
    }
    else
    {
        --decimals;
        X >>= 1;
    }
    return X & 1;
}

double BernoulliRand::Entropy()
{
    return -(p * std::log(p) + q * std::log(q));
}
