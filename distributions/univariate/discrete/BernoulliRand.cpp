#include "BernoulliRand.h"
#include "../continuous/UniformRand.h"
#include "../BasicRandGenerator.h"

BernoulliRand::BernoulliRand(double probability) : BinomialRand(1, probability)
{
    boundary = q * RandGenerator::maxValue();
}

std::string BernoulliRand::Name() const
{
    return "Bernoulli(" + toStringWithPrecision(GetProbability()) + ")";
}

void BernoulliRand::SetProbability(double probability)
{
    SetParameters(1, probability);
    boundary = q * RandGenerator::maxValue();
}

double BernoulliRand::P(const int & k) const
{
    return (k == 0) ? q : ((k == 1) ? p : 0);
}

double BernoulliRand::logP(const int & k) const
{
    return (k == 0) ? log1mProb : ((k == 1) ? logProb : 0);
}

double BernoulliRand::F(const int & k) const
{
    return (k < 0) ? 0.0 : ((k < 1) ? q : 1);
}

double BernoulliRand::S(const int & k) const
{
    return (k < 0) ? 1.0 : ((k < 1) ? p : 0.0);
}

int BernoulliRand::Variate() const
{
    return RandGenerator::variate() > boundary;
}

int BernoulliRand::Variate(double p)
{
    return UniformRand::StandardVariate() <= p;
}

int BernoulliRand::StandardVariate()
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
    return -(p * logProb + q * log1mProb);
}
