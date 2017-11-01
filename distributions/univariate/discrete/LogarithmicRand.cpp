#include "LogarithmicRand.h"
#include "../continuous/UniformRand.h"

LogarithmicRand::LogarithmicRand(double probability)
{
    SetProbability(probability);
}

String LogarithmicRand::Name() const
{
    return "Logarithmic(" + toStringWithPrecision(GetProbability()) + ")";
}

void LogarithmicRand::SetProbability(double probability)
{
    if (probability <= 0.0 || probability >= 1.0)
        throw std::invalid_argument("Logarithmic distribution: probability parameter should in interval (0, 1)");
    p = probability;
    logProb = std::log(p);
    log1mProb = std::log1p(-p);
}

double LogarithmicRand::P(const int & k) const
{
    return (k < 1) ? 0.0 : -std::pow(p, k) / (k * log1mProb);
}

double LogarithmicRand::logP(const int & k) const
{
    return (k < 1) ? -INFINITY : k * logProb - std::log(-k * log1mProb);
}

double LogarithmicRand::betaFun(int a) const
{
    double denom = a + 1;
    double sum = p * p / (a + 2) + p / (a + 1) + 1.0 / a;
    double add = 1;
    int i = 3;
    do {
        add = std::exp(i * logProb) / (++denom);
        sum += add;
        ++i;
    } while (add > MIN_POSITIVE * sum);
    return std::exp(a * logProb) * sum;
}

double LogarithmicRand::F(const int & k) const
{
    return (k < 1) ? 0.0 : 1 + betaFun(k + 1) / log1mProb;
}

double LogarithmicRand::S(const int & k) const
{
    return (k < 1) ? 1.0 : -betaFun(k + 1) / log1mProb;
}

int LogarithmicRand::Variate() const
{
    /// Kemp's second accelerated generator
    /// p. 548, "Non-Uniform Random Variate Generation" by Luc Devroye
    double V = UniformRand::StandardVariate();
    if (V >= p)
        return 1.0;
    double U = UniformRand::StandardVariate();
    double y = -std::expm1(U * log1mProb);
    if (V > y)
        return 1.0;
    if (V > y * y)
        return 2.0;
    return std::floor(1.0 + std::log(V) / std::log(y));
}

double LogarithmicRand::Mean() const
{
    return -p / (1.0 - p) / log1mProb;
}

double LogarithmicRand::Variance() const
{
    double var = p / log1mProb + 1;
    var /= log1mProb;
    var *= p;
    double q = 1.0 - p;
    return -var / (q * q);
}

std::complex<double> LogarithmicRand::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t), std::sin(t));
    y = std::log(1.0 - p * y);
    return y / log1mProb;
}

int LogarithmicRand::Mode() const
{
    return 1.0;
}
