#include "UniformRand.h"
#include "../BasicRandGenerator.h"

UniformRand::UniformRand(double minValue, double maxValue)
{
    setSupport(minValue, maxValue);
}

std::string UniformRand::name()
{
    return "Uniform(" + toStringWithPrecision(getMinValue()) + ", " + toStringWithPrecision(getMaxValue()) + ")";
}

void UniformRand::setSupport(double minValue, double maxValue)
{
    BetaRand::setSupport(minValue, maxValue);
    bmaInv = 1.0 / (b - a);
}

double UniformRand::f(double x) const
{
    return (x >= a && x <= b) ? bmaInv : 0;
}

double UniformRand::F(double x) const
{
    if (x < a)
        return 0;
    return (x > b) ? 1 : bmaInv * (x - a);
}

double UniformRand::variate() const
{
    return a + standardVariate() * bma;
}

double UniformRand::variate(double minValue, double maxValue)
{
    return minValue + standardVariate() * (maxValue - minValue);
}

double UniformRand::standardVariate()
{
    double x;
    unsigned long long a = RandGenerator::variate();
    a = (a >> 12) | 0x3FF0000000000000ULL; /// Take upper 52 bits
    *((unsigned long long *)&x) = a; /// Make a double from bits
    return x - 1.0;
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
    std::complex<double> x(0, -t * a), y(0, -t * b);
    std::complex<double> numerator = std::exp(x) - std::exp(y);
    std::complex<double> denominator(0, t * bma);
    return numerator / denominator;
}

double UniformRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    return a + bma * p;
}

double UniformRand::Median() const
{
    return .5 * (b + a);
}

double UniformRand::Mode() const
{
    /// this can be any value in [a, b]
    return variate();
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

bool UniformRand::fitMinimumMLE(const std::vector<double> &sample)
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
    setSupport(minVar, b);
    return true;
}

bool UniformRand::fitMaximumMLE(const std::vector<double> &sample)
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
    setSupport(a, maxVar);
    return true;
}

bool UniformRand::fitSupportMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;
    double maxVar = sample.at(0), minVar = maxVar;
    for (double var : sample) {
        maxVar = std::max(var, maxVar);
        minVar = std::min(var, minVar);
    }
    setSupport(minVar, maxVar);
    return true;
}

bool UniformRand::fitMinimumMM(const std::vector<double> &sample)
{
    double m = RandMath::sampleMean(sample);
    setSupport(m + m - b, b);
    return true;
}

bool UniformRand::fitMaximumMM(const std::vector<double> &sample)
{
    double m = RandMath::sampleMean(sample);
    setSupport(a, m + m - a);
    return true;
}

bool UniformRand::fitSupportMM(const std::vector<double> &sample)
{
    double mean = RandMath::sampleMean(sample);
    double var = RandMath::sampleVariance(sample, mean);
    double s = std::sqrt(3 * var);
    setSupport(mean - s, mean + s);
    return true;
}

bool UniformRand::fitMinimumUMVU(const std::vector<double> &sample)
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
    setSupport(minVar, b);
    return true;
}

bool UniformRand::fitMaximumUMVU(const std::vector<double> &sample)
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
    setSupport(a, maxVar);
    return true;
}

bool UniformRand::fitSupportUMVU(const std::vector<double> &sample)
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

    setSupport(minVar, maxVar);
    return true;
}
