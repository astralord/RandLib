#include "UniformDiscreteRand.h"

UniformDiscreteRand::UniformDiscreteRand(int minValue, int maxValue)
{
    SetBoundaries(minValue, maxValue);
}

std::string UniformDiscreteRand::Name() const
{
    return "Uniform Discrete(" + toStringWithPrecision(MinValue()) + ", " + toStringWithPrecision(MaxValue()) + ")";
}

void UniformDiscreteRand::SetBoundaries(int minValue, int maxValue)
{
    a = minValue;
    b = maxValue;

    if (b < a)
        b = a + 1;

    n = b - a + 1;
    nInv = 1.0 / n;
    logN = std::log(n);
}

double UniformDiscreteRand::P(int k) const
{
    return (k < a || k > b) ? 0.0 : nInv;
}

double UniformDiscreteRand::logP(int k) const
{
    return (k < a || k > b) ? -INFINITY : -logN;
}

double UniformDiscreteRand::F(int k) const
{
    if (k < a)
        return 0.0;
    if (k > b)
        return 1.0;
    return (k - a + 1) * nInv;
}

int UniformDiscreteRand::Variate() const
{
    double x = (RandGenerator::variate() % n);
    return a + x;
}

double UniformDiscreteRand::Mean() const
{
    return 0.5 * (b + a);
}

double UniformDiscreteRand::Variance() const
{
    return (n * n - 1.0) / 12;
}

double UniformDiscreteRand::Median() const
{
    return 0.5 * (b + a);
}

int UniformDiscreteRand::Mode() const
{
    return 0.5;
}

double UniformDiscreteRand::Skewness() const
{
    return 0.0;
}

double UniformDiscreteRand::ExcessKurtosis() const
{
    double res = 2.0 / (n * n - 1);
    ++res;
    return -1.2 * res;
}

std::complex<double> UniformDiscreteRand::CFImpl(double t) const
{
    double at = a * t;
    double bp1t = (b + 1) * t;
    double reNum = std::cos(at) - std::cos(bp1t);
    double imNum = std::sin(at) - std::sin(bp1t);
    std::complex<double> numen(reNum, imNum);
    std::complex<double>denom(1.0 - std::cos(t), -std::sin(t));
    return nInv * numen / denom;
}

double UniformDiscreteRand::Entropy() const
{
    return logN;
}

double UniformDiscreteRand::Likelihood(const std::vector<int> &sample) const
{
    for (const int & var : sample) {
        if (var < a || var > b)
            return 0.0;
    }
    return std::pow(n, -sample.size());
}

double UniformDiscreteRand::LogLikelihood(const std::vector<int> &sample) const
{
    for (const int & var : sample) {
        if (var < a || var > b)
            return -INFINITY;
    }
    return -sample.size() * logN;
}
