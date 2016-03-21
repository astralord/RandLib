#include "GeometricRand.h"
#include "../continuous/UniformRand.h"

GeometricRand::GeometricRand(double probability) : PascalRand(1, probability)
{
}

std::string GeometricRand::name()
{
    return "Geometric(" + toStringWithPrecision(getProbability()) + ")";
}

void GeometricRand::setProbability(double probability)
{
    setParameters(1, probability);
}

double GeometricRand::P(int k) const
{
    return (k < 0) ? 0 : p * std::pow(q, k);
}

double GeometricRand::F(double x) const
{
    return (x < 0) ? 0 : 1 - std::pow(q, std::floor(x) + 1);
}

double GeometricRand::variate() const
{
    return (p < 0.2) ? variateGeometricThroughExponential() : variateGeometricByTable();
}

double GeometricRand::variate(double probability)
{
    /// here we use 0.1 instead of 0.2 because log(1-p) wasn't hashed
    if (probability < 0.1)
        return std::floor(ExponentialRand::variate(-std::log(1 - probability)));

    double U = UniformRand::standardVariate();
    int x = 0;
    double prod = probability, sum = prod, q = 1 - probability;
    while (U > sum) {
        prod *= q;
        sum += prod;
        ++x;
    }
    return x;
}

void GeometricRand::sample(QVector<double> &outputData) const
{
    if (p < 0.2) {
        for (double &var : outputData)
            var = variateGeometricThroughExponential();
    }
    else {
        for (double &var : outputData)
            var = variateGeometricByTable();
    }
}

double GeometricRand::Median() const
{
    return std::floor(-M_LN2 / std::log(q));
}

double GeometricRand::Entropy() const
{
    double a = -q * std::log(q);
    double b = -p * std::log(p);
    return (a + b) / (M_LN2 * p);
}

bool GeometricRand::checkValidity(const QVector<double> &sample)
{
    for (int var : sample) {
        if (var < 0)
            return false;
    }
    return true;
}

bool GeometricRand::fitProbabilityMLE(const QVector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setProbability(1.0 / (RandMath::sampleMean(sample) + 1));
    return true;
}

bool GeometricRand::fitProbabilityMM(const QVector<double> &sample)
{
    return fitProbabilityMLE(sample);
}

bool GeometricRand::fitProbabilityBayes(const QVector<double> &sample, BetaRand &priorDistribution)
{
    int n = sample.size();
    if (!checkValidity(sample))
        return false;
    double alpha = priorDistribution.getAlpha();
    double beta = priorDistribution.getBeta();
    priorDistribution.setShapes(alpha + n, beta + RandMath::sum(sample));
    setProbability(priorDistribution.Mean());
    return true;
}
