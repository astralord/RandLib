#include "randomvariable.h"
#include <time.h>

RandomVariable::RandomVariable()
{
    /// maybe we should do it in more elegant way
    seedRand(time(0));
}

void RandomVariable::seedRand(unsigned long seed)
{
    randValue ^= seed;
}

unsigned long RandomVariable::SHR3()
{
    randValue ^= (randValue << 17);
    randValue ^= (randValue >> 13);
    randValue ^= (randValue << 5);
    return randValue;
}

void RandomVariable::_cdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = cdf(x[i]);
}

void RandomVariable::_cdf(const QVector<double> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = cdf(x[i]);
}

void RandomVariable::sample(QVector<double> &outputData)
{
    for (int i = 0; i != outputData.size(); ++i)
        outputData[i] = value();
}

void RandomVariable::sample(std::vector<double> &outputData) {
    for (size_t i = 0; i != outputData.size(); ++i)
        outputData[i] = value();
}

void ContinuousRand::_pdf(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
}

void ContinuousRand::_pdf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
}

double ContinuousRand::likelihood(const std::vector<double> &sample) const
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= pdf(sample[i]);
    return res;
}

double ContinuousRand::loglikelihood(const std::vector<double> &sample) const
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += std::log(pdf(sample[i]));
    return res;
}


void DiscreteIntRand::probs(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteIntRand::probs(const QVector<int> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}


void DiscreteDoubleRand::probs(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteDoubleRand::probs(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}
