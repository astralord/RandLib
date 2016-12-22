#include "GeometricRand.h"
#include "../continuous/UniformRand.h"

GeometricRand::GeometricRand(double probability) : PascalRand(1, probability)
{
}

std::string GeometricRand::Name() const
{
    return "Geometric(" + toStringWithPrecision(GetProbability()) + ")";
}

void GeometricRand::SetProbability(double probability)
{
    SetParameters(1, probability);
}

double GeometricRand::P(int k) const
{
    return (k < 0) ? 0 : p * std::exp(k * logQ);
}

double GeometricRand::F(int k) const
{
    return (k < 0) ? 0 : 1 - std::exp((k + 1) * logQ);
}

int GeometricRand::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == EXPONENTIAL)
        return variateGeometricThroughExponential();
    if (genId == TABLE)
        return variateGeometricByTable();
    return -1; /// unexpected return
}

int GeometricRand::Variate(double probability)
{
    /// here we use 0.05 instead of 0.08 because log(q) wasn't hashed
    if (probability < 0.05)
        return std::floor(ExponentialRand::Variate(-std::log1p(-probability)));

    double U = UniformRand::StandardVariate();
    int x = 0;
    double prod = probability, sum = prod, q = 1 - probability;
    while (U > sum) {
        prod *= q;
        sum += prod;
        ++x;
    }
    return x;
}

void GeometricRand::Sample(std::vector<int> &outputData) const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == EXPONENTIAL) {
        for (int &var : outputData)
            var = variateGeometricThroughExponential();
    }
    else if (genId == TABLE) {
        for (int &var : outputData)
            var = variateGeometricByTable();
    }
}

double GeometricRand::Median() const
{
    return std::floor(-M_LN2 / logQ);
}

double GeometricRand::Entropy() const
{
    double a = -q * logQ;
    double b = -p * logP;
    return (a + b) / (M_LN2 * p);
}

bool GeometricRand::FitMLE(const std::vector<int> &sample)
{
    if (!checkValidity(sample))
        return false;
    SetProbability(1.0 / (sampleMean(sample) + 1));
    return true;
}

bool GeometricRand::FitMM(const std::vector<int> &sample)
{
    return FitMLE(sample);
}

bool GeometricRand::FitBayes(const std::vector<int> &sample, BetaRand &priorDistribution)
{
    int n = sample.size();
    if (!checkValidity(sample))
        return false;
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    priorDistribution.SetParameters(alpha + n, beta + sampleSum(sample));
    SetProbability(priorDistribution.Mean());
    return true;
}
