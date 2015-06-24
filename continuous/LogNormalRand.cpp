#include "LogNormalRand.h"

LogNormalRand::LogNormalRand(double location, double scale)
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
    X.setSigma(scale);
    expVar = std::exp(X.Var());
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
    double expVarSq = expVar * expVar;
    double res = expVarSq; /// exp(2s^2)
    res += 2 * expVar + 3; /// exp(2s^2) + 2exp(s^2) + 3
    res *= expVarSq; /// exp(4s^2) + 2exp(3s^2) + 3exp(2s^2)
    return res - 6;
}

bool LogNormalRand::fitToData(const QVector<double> &sample)
{
    double average = 0.0, deviation = 0.0;
    foreach (double var, sample) {
        if (var <= 0)
            return false;
        average += std::log(var);
    }
    average /= sample.size();
    if (std::isnan(average) || std::isinf(average))
        return false;

    foreach (double var, sample) {
        double currDev = (std::log(var) - average);
        deviation += currDev * currDev;
    }
    deviation /= (sample.size() - 1);
    if (std::isnan(deviation) || std::isinf(deviation))
        return false;

    setLocation(average);
    setScale(std::sqrt(deviation));
    return true;
}
