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
    //z ^= randValue;
    //w ^= z;
}

void RandomVariable::_pdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
}

void RandomVariable::_cdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = cdf(x[i]);
}

void RandomVariable::_pdf(const QVector<double> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
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

double RandomVariable::likelihood(const std::vector<double> &sample)
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= pdf(sample[i]);
    return res;
}

double RandomVariable::loglikelihood(const std::vector<double> &sample)
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += qLn(pdf(sample[i]));
    return res;
}

unsigned long RandomVariable::SHR3()
{
    randValue ^= (randValue << 17);
    randValue ^= (randValue >> 13);
    randValue ^= (randValue << 5);
    return randValue;
}

unsigned long long RandomVariable::MWC()
{
    /*z = 36969 * (z & 65535) + (z >> 16);
    w = 18000 * (w & 65535) + (w >> 16);
    return w + (z << 16);*/
    return 0;
}
