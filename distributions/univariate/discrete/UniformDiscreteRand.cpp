#include "UniformDiscreteRand.h"

UniformDiscreteRand::UniformDiscreteRand(int minValue, int maxValue)
{
    SetBoundaries(minValue, maxValue);
}

String UniformDiscreteRand::Name() const
{
    return "Uniform Discrete(" + toStringWithPrecision(MinValue()) + ", " + toStringWithPrecision(MaxValue()) + ")";
}

void UniformDiscreteRand::SetBoundaries(int minValue, int maxValue)
{
    if (minValue >= maxValue)
        throw std::invalid_argument("Uniform discrete distribution: minimal value should be less than maximum value");

    a = minValue;
    b = maxValue;

    n = b - a + 1;
    nInv = 1.0 / n;
    logN = std::log(n);

    unsigned long long MAX_RAND = RandGenerator::MaxValue();
    MAX_RAND_UNBIASED = MAX_RAND - MAX_RAND % n - 1;
}

double UniformDiscreteRand::P(const int & k) const
{
    return (k < a || k > b) ? 0.0 : nInv;
}

double UniformDiscreteRand::logP(const int & k) const
{
    return (k < a || k > b) ? -INFINITY : -logN;
}

double UniformDiscreteRand::F(const int & k) const
{
    if (k < a)
        return 0.0;
    if (k > b)
        return 1.0;
    return (k - a + 1) * nInv;
}

int UniformDiscreteRand::Variate() const
{
    unsigned long intVar;
    do {
        intVar = RandGenerator::Variate();
    } while (intVar > MAX_RAND_UNBIASED);
    return a + (intVar % n);
}

double UniformDiscreteRand::Mean() const
{
    return 0.5 * (b + a);
}

double UniformDiscreteRand::Variance() const
{
    double nm1 = n - 1;
    double np1 = n + 1;
    return nm1 * np1 / 12;
}

int UniformDiscreteRand::Median() const
{
    return (b + a) >> 1;
}

int UniformDiscreteRand::Mode() const
{
    /// this can be any value in [a, b]
    return 0.5 * (a + b);
}

double UniformDiscreteRand::Skewness() const
{
    return 0.0;
}

double UniformDiscreteRand::ExcessKurtosis() const
{
    double kurt = n;
    kurt *= n;
    --kurt;
    kurt = 2.0 / kurt;
    ++kurt;
    return -1.2 * kurt;
}

std::complex<double> UniformDiscreteRand::CFImpl(double t) const
{
    double at = a * t;
    double bp1t = (b + 1) * t;
    double reNum = std::cos(at) - std::cos(bp1t);
    double imNum = std::sin(at) - std::sin(bp1t);
    std::complex<double> numerator(reNum, imNum);
    std::complex<double> denominator(1.0 - std::cos(t), -std::sin(t));
    return nInv * numerator / denominator;
}

double UniformDiscreteRand::Entropy() const
{
    return logN;
}

double UniformDiscreteRand::LikelihoodFunction(const std::vector<int> &sample) const
{
    for (const int & var : sample) {
        if (var < a || var > b)
            return 0.0;
    }
    return std::pow(n, -sample.size());
}

double UniformDiscreteRand::LogLikelihoodFunction(const std::vector<int> &sample) const
{
    for (const int & var : sample) {
        if (var < a || var > b)
            return -INFINITY;
    }
    return -sample.size() * logN;
}
