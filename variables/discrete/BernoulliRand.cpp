#include "BernoulliRand.h"
#include "../continuous/UniformRand.h"

double BernoulliRand::U = UniformRand::standardVariate();

BernoulliRand::BernoulliRand(double probability)
{
    setProbability(probability);
}

std::string BernoulliRand::name()
{
    return "Bernoulli(" + toStringWithPrecision(getProbability()) + ")";
}

void BernoulliRand::setProbability(double probability)
{
    p = std::min(std::max(probability, 0.0), 1.0);
    q = 1.0 - p;

    boundary = q * RandGenerator::maxValue();
}

double BernoulliRand::P(int k) const
{
    return (k == 0) ? q : ((k == 1) ? p : 0);
}

double BernoulliRand::F(double x) const
{
    return (x < 0) ? 0 : ((x < 1) ? q : 1);
}

double BernoulliRand::variate() const
{
    return RandGenerator::variate() > boundary;
}

double BernoulliRand::variate(double p)
{
    return UniformRand::standardVariate() > 1.0 - p;
}

double BernoulliRand::standardVariate()
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

double BernoulliRand::Mean() const
{
    return p;
}

double BernoulliRand::Variance() const
{
    return p * (1 - p);
}

std::complex<double> BernoulliRand::CF(double t) const
{
    return std::complex<double>(q + p * std::cos(t), std::sin(t));
}

double BernoulliRand::Median() const
{
    return (p < 0.5) ? 0 : ((p > 0.5) ? 1 : 0.5);
}

double BernoulliRand::Mode() const
{
    /// if q == p -> this can be any of {0, 1}
    return (p < 0.5) ? 0 : ((p > 0.5) ? 1 : variate());
}

double BernoulliRand::Skewness() const
{
    return (q - p) / std::sqrt(p * q);
}

double BernoulliRand::ExcessKurtosis() const
{
    return 1.0 / (p * q) - 6;
}

double BernoulliRand::Entropy()
{
    return -(p * std::log(p) + q * std::log(q));
}
