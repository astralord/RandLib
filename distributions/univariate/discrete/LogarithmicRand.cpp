#include "LogarithmicRand.h"
#include "../continuous/UniformRand.h"

LogarithmicRand::LogarithmicRand(double probability)
{
    setProbability(probability);
}

std::string LogarithmicRand::name() const
{
    return "Logarithmic(" + toStringWithPrecision(getProbability()) + ")";
}

void LogarithmicRand::setProbability(double probability)
{
    p = probability;
    if (p <= 0 || p >= 1)
        p = 0.5;
    q = 1.0 - p;
    logQInv = 1.0 / std::log(q);
}

double LogarithmicRand::P(int k) const
{
    return (k < 1) ? 0 : -logQInv * std::pow(p, k) / k;
}

double LogarithmicRand::F(int k) const
{
    return (k < 1) ? 0 : 1 + logQInv * RandMath::incompleteBetaFun(p, k + 1, 0);
}

int LogarithmicRand::variate() const
{
    /// Kemp's second accelerated generator
    /// p. 548, "Non-Uniform Random Variate Generation" by Luc Devroye
    double V = UniformRand::standardVariate();
    if (V >= p)
        return 1.0;
    double U = UniformRand::standardVariate();
    double y = 1.0 - std::exp(U / logQInv);
    if (V > y)
        return 1.0;
    if (V > y * y)
        return 2.0;
    return std::floor(1.0 + std::log(V) / std::log(y));
}

double LogarithmicRand::Mean() const
{
    return -logQInv * p / q;
}

double LogarithmicRand::Variance() const
{
    double var = p * logQInv + 1;
    var *= logQInv;
    var *= p;
    return -var / (q * q);
}

std::complex<double> LogarithmicRand::CF(double t) const
{
    if (t == 0)
        return 1;
    std::complex<double> y(0, t);
    y = std::exp(y);
    y = std::log(1.0 - p * y);
    return y * logQInv;
}

int LogarithmicRand::Mode() const
{
    return 1.0;
}
