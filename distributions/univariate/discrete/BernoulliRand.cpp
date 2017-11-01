#include "BernoulliRand.h"
#include "../continuous/UniformRand.h"
#include "../BasicRandGenerator.h"

BernoulliRand::BernoulliRand(double probability) : BinomialDistribution(1, probability)
{
    boundary = q * RandGenerator::MaxValue();
}

String BernoulliRand::Name() const
{
    return "Bernoulli(" + toStringWithPrecision(GetProbability()) + ")";
}

void BernoulliRand::SetProbability(double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Bernoulli distribution: probability parameter should in interval [0, 1]");
    SetParameters(1, probability);
    boundary = q * RandGenerator::MaxValue();
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
    return RandGenerator::Variate() > boundary;
}

int BernoulliRand::Variate(double probability)
{
    return (probability < 0.0 || probability > 1.0) ? -1 : UniformRand::StandardVariate() <= probability;
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
        X = RandGenerator::Variate();
    }
    else
    {
        --decimals;
        X >>= 1;
    }
    return X & 1;
}

void BernoulliRand::Sample(std::vector<int> &outputData) const
{
    if (p == 0.5) {
        for (int & var : outputData)
            var = StandardVariate();
    }
    else {
        for (int & var : outputData)
            var = this->Variate();
    }
}

double BernoulliRand::Entropy()
{
    return -(p * logProb + q * log1mProb);
}
