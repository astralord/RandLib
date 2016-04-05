#include "ChiRand.h"

ChiRand::ChiRand(int degree)
{
    setDegree(degree);
}

void ChiRand::setDegree(int degree)
{
    X.setDegree(degree);
}

double ChiRand::f(double x) const
{
    return (x > 0) ? 2 * x * X.f(x * x) : 0.0;
}

double ChiRand::F(double x) const
{
    return (x > 0) ? X.F(x * x) : 0.0;
}

double ChiRand::variate() const
{
    return std::sqrt(X.variate());
}

void ChiRand::sample(std::vector<double> &outputData) const
{
    X.sample(outputData);
    for (double & var : outputData)
        var = std::sqrt(var);
}

double ChiRand::varianceImpl(double mean) const
{
    return v - mean * mean;
}

double ChiRand::skewnessImpl(double mean, double sigma) const
{
    double variance = sigma * sigma;
    double y = mean * (1 - 2 * variance);
    return y / (sigma * variance);
}

double ChiRand::Mean() const
{
    double halfV = 0.5 * v;
    double x = std::lgamma(halfV + 0.5);
    double y = std::lgamma(halfV);
    return M_SQRT2 * std::exp(x - y);
}

double ChiRand::Variance() const
{
    return varianceImpl(Mean());
}

double ChiRand::Mode() const
{
    return (v > 1) ? std::sqrt(v - 1) : 0;
}

double ChiRand::Skewness() const
{
    double mean = Mean();
    return skewnessImpl(mean, std::sqrt(varianceImpl(mean)));
}

double ChiRand::ExcessKurtosis() const
{
    double mean = Mean();
    double variance = varianceImpl(mean);
    double sigma = std::sqrt(variance);
    double skewness = skewnessImpl(mean, sigma);
    double y = 1 - mean * sigma * skewness;
    y /= variance;
    --y;
    return y + y;
}
