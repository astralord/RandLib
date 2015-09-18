#include "LogarithmicRand.h"

LogarithmicRand::LogarithmicRand(double probability)
{
    setProbability(probability);
}

std::string LogarithmicRand::name()
{
    return "LogarithmicRand(" + toStringWithPrecision(getProbability()) + ")";
}

void LogarithmicRand::setProbability(double probability)
{
    p = std::min(probability, 1.0 - 1e-14);
    if (p <= 0)
        p = MIN_POSITIVE;
    q = 1.0 - p;
    logQInv = 1.0 / std::log(1.0 - p);
}

double LogarithmicRand::P(int k) const
{
    return (k < 1) ? 0 : -logQInv * std::pow(p, k) / k;
}

double LogarithmicRand::F(double x) const
{
    double k = std::floor(x);
    return (k < 1) ? 0 : 1 + logQInv * RandMath::incompleteBetaFun(p, k + 1, 0);
}

double LogarithmicRand::variate() const
{
    //TODO
    return -1.0;
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
    std::complex<double> y(0, t);
    y = std::exp(y);
    y = std::log(1.0 - p * y);
    return y * logQInv;
}

double LogarithmicRand::Mode() const
{
    return 1.0;
}
