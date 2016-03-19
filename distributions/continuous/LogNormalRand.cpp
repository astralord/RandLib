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
    X.setLocation(location);
    expMu = std::exp(X.Mean());
}

void LogNormalRand::setScale(double scale)
{
    X.setScale(scale);
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

bool LogNormalRand::checkValidity(const QVector<double> &sample)
{
    for (double var : sample) {
        if (var <= 0)
            return false;
    }
    return true;
}

bool LogNormalRand::fitLocation_MM(const QVector<double> &sample)
{
    double average = RandMath::sampleMean(sample);
    double var = X.getVar();
    setLocation(std::log(average) - 0.5 * var);
    return true;
}

bool LogNormalRand::fitScale_MM(const QVector<double> &sample)
{
    double average = RandMath::sampleMean(sample);
    double mu = X.getLocation();
    double aux = std::log(average) - mu;
    setScale(std::sqrt(aux + aux));
    return true;
}

bool LogNormalRand::fit_MLE(const QVector<double> &sample)
{
    int N = sample.size();
    if (N == 0)
        return false;

    /// Calculate location
    long double logMean = 0.0L;
    for (double var : sample) {
        logMean += std::log(var);
    }
    logMean /= N;

    /// Calculate scale
    long double deviation = 0.0L;
    for (double var : sample) {
        double currDev = (std::log(var) - logMean);
        deviation += currDev * currDev;
    }
    deviation /= N;

    setLocation(logMean);
    setScale(std::sqrt(deviation));
    return true;
}
