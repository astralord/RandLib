#include "LogNormalRand.h"

LogNormalRand::LogNormalRand(double location, double scale) :
    X(0, 1)
{
    setLocation(location);
    setScale(scale);
}

void LogNormalRand::setLocation(double location)
{
    X.setMean(location);
    expMu = std::exp(X.M());
}

void LogNormalRand::setScale(double scale)
{
    X.setSigma(scale * scale);
    expSigmaSq = std::exp(X.Var());
}

double LogNormalRand::f(double x) const
{
    return (x > 0) ? X.f(std::log(x)) / x : 0;
}

double LogNormalRand::F(double x) const
{
    return (x > 0) ? X.F(std::log(x)) : 0;
}

double LogNormalRand::value()
{
    return std::exp(X.value());
}

double LogNormalRand::ExcessKurtosis() const
{
    double expSigma2Sq = expSigmaSq * expSigmaSq;
    double res = expSigma2Sq; /// exp(2s^2)
    res += 2 * expSigmaSq + 3; /// exp(2s^2) + 2exp(s^2) + 3
    res *= expSigma2Sq; /// exp(4s^2) + 2exp(3s^2) + 3exp(2s^2)
    return res - 6;
}
