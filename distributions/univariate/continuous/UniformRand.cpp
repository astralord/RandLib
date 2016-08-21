#include "UniformRand.h"
#include "../BasicRandGenerator.h"

UniformRand::UniformRand(double minValue, double maxValue)
{
    SetSupport(minValue, maxValue);
}

std::string UniformRand::Name() const
{
    return "Uniform(" + toStringWithPrecision(GetMinValue()) + ", " + toStringWithPrecision(GetMaxValue()) + ")";
}

void UniformRand::SetSupport(double minValue, double maxValue)
{
    BetaRand::SetParameters(1, 1, minValue, maxValue);
    bmaInv = 1.0 / (b - a);
}

double UniformRand::f(double x) const
{
    return (x < a || x > b) ? 0 : bmaInv;
}

double UniformRand::F(double x) const
{
    if (x < a)
        return 0;
    return (x > b) ? 1 : bmaInv * (x - a);
}

double UniformRand::Variate() const
{
    return a + StandardVariate() * bma;
}

double UniformRand::Variate(double minValue, double maxValue)
{
    return minValue + StandardVariate() * (maxValue - minValue);
}

double UniformRand::StandardVariate()
{
#ifdef JKISS32RAND
    return RandGenerator::variate() / 4294967296.0;
#else
    double x;
    unsigned long long a = RandGenerator::variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bits
    *((unsigned long long *)&x) = a; /// Make a double from bits
    return x - 1.0;
#endif
}

double UniformRand::Mean() const
{
    return .5 * (b + a);
}

double UniformRand::Variance() const
{
    return bma * bma / 12;
}

std::complex<double> UniformRand::CF(double t) const
{
    if (t == 0)
        return 1;
    std::complex<double> x(0, t * a), y(0, t * b);
    std::complex<double> numerator = std::exp(y) - std::exp(x);
    std::complex<double> denominator(0, t * bma);
    return numerator / denominator;
}

double UniformRand::quantileImpl(double p) const
{
    return a + bma * p;
}

double UniformRand::quantileImpl1m(double p) const
{
    return a + bma - bma * p;
}

double UniformRand::Median() const
{
    return .5 * (b + a);
}

double UniformRand::Mode() const
{
    /// this can be any value in [a, b]
    return Variate();
}

double UniformRand::Skewness() const
{
    return 0.0;
}

double UniformRand::ExcessKurtosis() const
{
    return -1.2;
}

double UniformRand::Entropy() const
{
    return (b == a) ? -INFINITY : std::log(bma);
}

bool UniformRand::FitMinimumMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double minVar = sample.at(0);
    for (double var : sample) {
        if (var > b)
            return false;
        minVar = std::min(var, minVar);
    }
    SetSupport(minVar, b);
    return true;
}

bool UniformRand::FitMaximumMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double maxVar = sample.at(0);
    for (double var : sample) {
        if (var < a)
            return false;
        maxVar = std::max(var, maxVar);
    }
    SetSupport(a, maxVar);
    return true;
}

bool UniformRand::FitMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double maxVar = sample.at(0), minVar = maxVar;
    for (double var : sample) {
        maxVar = std::max(var, maxVar);
        minVar = std::min(var, minVar);
    }
    SetSupport(minVar, maxVar);
    return true;
}

bool UniformRand::FitMinimumMM(const std::vector<double> &sample)
{
    double m = sampleMean(sample);
    SetSupport(m + m - b, b);
    return true;
}

bool UniformRand::FitMaximumMM(const std::vector<double> &sample)
{
    double m = sampleMean(sample);
    SetSupport(a, m + m - a);
    return true;
}

bool UniformRand::FitMM(const std::vector<double> &sample)
{
    double mean = sampleMean(sample);
    double var = sampleVariance(sample, mean);
    double s = std::sqrt(3 * var);
    SetSupport(mean - s, mean + s);
    return true;
}

bool UniformRand::FitMinimumUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double minVar = sample.at(0);
    for (double var : sample) {
        if (var < a)
            return false;
        minVar = std::min(var, minVar);
    }
    
    /// E[min] = b - n / (n + 1) * (b - a)
    minVar = (minVar * (n + 1) - b) / n;
    SetSupport(minVar, b);
    return true;
}

bool UniformRand::FitMaximumUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double maxVar = sample.at(0);
    for (double var : sample) {
        if (var < a)
            return false;
        maxVar = std::max(var, maxVar);
    }
    
    /// E[max] = (b - a) * n / (n + 1) + a
    maxVar = (maxVar * (n + 1) - a) / n;
    SetSupport(a, maxVar);
    return true;
}

bool UniformRand::FitUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double maxVar = sample.at(0), minVar = maxVar;
    for (double var : sample) {
        maxVar = std::max(var, maxVar);
        minVar = std::min(var, minVar);
    }

    /// E[min] = b - n / (n + 1) * (b - a)
    /// E[max] = (b - a) * n / (n + 1) + a

    minVar = (minVar * n - maxVar) / (n - 1);
    maxVar = (maxVar * (n + 1) - minVar) / n;

    SetSupport(minVar, maxVar);
    return true;
}
