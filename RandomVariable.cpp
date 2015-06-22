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

double RandomVariable::factorial(int n)
{
    double res = 1.0;
    for (int i = 2; i <= n; ++i)
        res *= i;
    return res;
}

double RandomVariable::doubleFactorial(int n)
{
    double res = 1.0;
    for (int i = n % 2 + 2; i <= n; i += 2)
        res *= i;
    return res;
}

double RandomVariable::lowerIncGamma(double a, double x)
{
    double sum = 0;
    double term = 1.0 / a;
    int n = 1;
    while (term > MIN_POSITIVE)
    {
        sum = sum + term;
        term *= (x / (a + n));
        ++n;
    }
    return std::pow(x, a) * std::exp(-x) * sum;
}

double RandomVariable::upperIncGamma(double a, double x)
{
    // TODO: find useful approximation
    return std::tgamma(a) - lowerIncGamma(a, x);
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

void ContinuousRand::_pdf(const std::vector<double> &x, std::vector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
}

void ContinuousRand::_pdf(const QVector<double> &x, QVector<double> &y)
{
    size_t size = std::min(x.size(), y.size());
    for (size_t i = 0; i != size; ++i)
        y[i] = pdf(x[i]);
}

double ContinuousRand::likelihood(const std::vector<double> &sample)
{
    double res = 1.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res *= pdf(sample[i]);
    return res;
}

double ContinuousRand::loglikelihood(const std::vector<double> &sample)
{
    double res = 0.0;
    for (size_t i = 0; i != sample.size(); ++i)
        res += std::log(pdf(sample[i]));
    return res;
}
