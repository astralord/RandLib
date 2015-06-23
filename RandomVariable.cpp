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

void RandomVariable::sample(QVector<double> &outputData)
{
    for (int i = 0; i != outputData.size(); ++i)
        outputData[i] = value();
}

void RandomVariable::sample(std::vector<double> &outputData) {
    for (size_t i = 0; i != outputData.size(); ++i)
        outputData[i] = value();
}


/// CONTINUOUS

void ContinuousRand::pdf(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

void ContinuousRand::pdf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = f(x[i]);
}

void ContinuousRand::cdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

void ContinuousRand::cdf(const QVector<double> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

double ContinuousRand::likelihood(const std::vector<double> &sample) const
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= f(sample[i]);
    return res;
}

double ContinuousRand::loglikelihood(const std::vector<double> &sample) const
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += std::log(f(sample[i]));
    return res;
}


/// DISCRETE INTEGER

void DiscreteIntRand::pmf(const std::vector<int> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteIntRand::pmf(const QVector<int> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteIntRand::cdf(const std::vector<int> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

void DiscreteIntRand::cdf(const QVector<int> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

double DiscreteIntRand::likelihood(const std::vector<int> &sample) const
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= P(sample[i]);
    return res;
}

double DiscreteIntRand::loglikelihood(const std::vector<int> &sample) const
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += std::log(P(sample[i]));
    return res;
}


/// DISCRETE DOUBLE

void DiscreteDoubleRand::pmf(const std::vector<double> &x, std::vector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteDoubleRand::pmf(const QVector<double> &x, QVector<double> &y) const
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = P(x[i]);
}

void DiscreteDoubleRand::cdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}

void DiscreteDoubleRand::cdf(const QVector<double> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = F(x[i]);
}


double DiscreteDoubleRand::likelihood(const std::vector<double> &sample) const
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= P(sample[i]);
    return res;
}

double DiscreteDoubleRand::loglikelihood(const std::vector<double> &sample) const
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += std::log(P(sample[i]));
    return res;
}
