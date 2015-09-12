#include "BernoulliRand.h"

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
    generatorEdge = q * RandGenerator::maxValue();
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
    return (RandGenerator::variate() < generatorEdge) ? 0 : 1;
}

double BernoulliRand::variate(double p)
{
    return (RandGenerator::variate() < (1.0 - p) * RandGenerator::maxValue()) ? 0 : 1;
}

std::complex<double> BernoulliRand::CF(double t) const
{
    return std::complex<double>(q + p * std::cos(t), std::sin(t));
}

double BernoulliRand::Quantile(double p) const
{
    if (p == 0)
        return -INFINITY;
    if (p == q)
        return 0;
    if (p == 1)
        return 1;
    return NAN;
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
