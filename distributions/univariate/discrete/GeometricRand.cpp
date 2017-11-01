#include "GeometricRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

String GeometricRand::Name() const
{
    return "Geometric(" + toStringWithPrecision(GetProbability()) + ")";
}

void GeometricRand::SetProbability(double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Geometric distribution: probability parameter should in interval [0, 1]");
    SetParameters(1, probability);
}

double GeometricRand::P(const int & k) const
{
    return (k < 0) ? 0 : p * std::exp(k * log1mProb);
}

double GeometricRand::logP(const int & k) const
{
    return (k < 0) ? -INFINITY : logProb + k * log1mProb;
}

double GeometricRand::F(const int & k) const
{
    return (k < 0) ? 0 : -std::expm1((k + 1) * log1mProb);
}

double GeometricRand::S(const int & k) const
{
    return (k < 0) ? 1 : std::exp((k + 1) * log1mProb);
}

int GeometricRand::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == EXPONENTIAL)
        return variateGeometricThroughExponential();
    if (genId == TABLE)
        return variateGeometricByTable();
    /// unexpected return
    return -1;
}

int GeometricRand::Variate(double probability)
{
    if (probability > 1.0 || probability < 0.0)
        return -1;

    /// here we use 0.05 instead of 0.08 because log(q) wasn't hashed
    if (probability < 0.05) {
        double rate = -std::log1p(-probability);
        double X = ExponentialRand::StandardVariate() / rate;
        return std::floor(X);
    }

    double U = UniformRand::StandardVariate();
    int x = 0;
    double prod = probability, sum = prod, qprob = 1.0 - probability;
    while (U > sum) {
        prod *= qprob;
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

int GeometricRand::Median() const
{
    double median = -M_LN2 / log1mProb;
    double flooredMedian = std::floor(median);
    return (RandMath::areClose(median, flooredMedian, 1e-8)) ? flooredMedian - 1 : flooredMedian;
}

double GeometricRand::Entropy() const
{
    double a = -q * log1mProb;
    double b = -p * logProb;
    return (a + b) / (M_LN2 * p);
}
