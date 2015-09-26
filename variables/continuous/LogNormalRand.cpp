#include "LogNormalRand.h"

LogNormalRand::LogNormalRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LogNormalRand::name()
{
    return "Log-Normal(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LogNormalRand::setLocation(double location)
{
    X.setMean(location);
    expMu = std::exp(X.Mean());
}

void LogNormalRand::setScale(double scale)
{
    X.setSigma(scale);
    expVar = std::exp(X.Variance());
}

double LogNormalRand::f(double x) const
{
    return (x > 0) ? X.f(std::log(x)) / x : 0;
}

double LogNormalRand::F(double x) const
{
    return (x > 0) ? X.F(std::log(x)) : 0;
}

double LogNormalRand::variate() const
{
    return std::exp(X.variate());
}

double LogNormalRand::Mean() const
{
    return expMu * std::sqrt(expVar);
}

double LogNormalRand::Variance() const
{
    return (expVar - 1) * expMu * expMu * expVar;
}

double LogNormalRand::Quantile(double p) const
{
    if (p == 0.0)
        return 0.0;
    double y = std::exp(X.Quantile(p));
    /// we try to be more accurate here due to scanty approximation of normal quantile
    double root = y;
    double var = Variance();
    if (RandMath::findRoot([this, p] (double x) {
        return F(x) - p;
    }, y - var, y + var, root))
        return root;
    return y;
}

double LogNormalRand::Median() const
{
    return expMu;
}

double LogNormalRand::Mode() const
{
    return expMu / expVar;
}

double LogNormalRand::Skewness() const
{
    return (expVar + 2) * std::sqrt(expVar - 1);
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
    if (sample.size() == 0)
        return false;

    /// Calculate location
    double average = 0.0;
    for (double var : sample) {
        if (var <= 0)
            return false;
        average += std::log(var);
    }
    average /= sample.size();

    /// Calculate scale
    double deviation = 0.0;
    for (double var : sample) {
        double currDev = (std::log(var) - average);
        deviation += currDev * currDev;
    }
    deviation /= std::max(sample.size() - 1, 1);

    setLocation(average);
    setScale(std::sqrt(deviation));
    return true;
}
