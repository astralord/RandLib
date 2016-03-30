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

bool LogNormalRand::checkValidity(const std::vector<double> &sample)
{
    for (double var : sample) {
        if (var <= 0)
            return false;
    }
    return true;
}

bool LogNormalRand::fitLocationMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = RandMath::sampleMean(sample);
    double var = X.getVariance();
    setLocation(std::log(average) - 0.5 * var);
    return true;
}

bool LogNormalRand::fitScaleMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = RandMath::sampleMean(sample);
    double mu = X.getLocation();
    double aux = std::log(average) - mu;
    setScale(std::sqrt(aux + aux));
    return true;
}

bool LogNormalRand::fitLocationAndScaleMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = RandMath::sampleMean(sample);
    double secondMoment = RandMath::rawMoment(sample, 2);
    double averageSq = average * average;
    setLocation(0.5 * std::log(averageSq * averageSq / secondMoment));
    setScale(std::sqrt(std::log(secondMoment / averageSq)));
    return true;
}

bool LogNormalRand::fitLocationMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;

    long double logMean = 0.0L;
    for (double var : sample) {
        logMean += std::log(var);
    }
    logMean /= n;

    setLocation(logMean);
    return true;
}

bool LogNormalRand::fitScaleMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;

    long double logVariance = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var);
        logVariance += logVar * logVar;
    }
    double mu = X.getLocation();
    logVariance /= n;
    logVariance -= mu * mu;

    setScale(std::sqrt(logVariance));
    return true;
}

bool LogNormalRand::fitLocationAndScaleMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;

    long double logMean = 0.0L;
    long double logVariance = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var);
        logMean += logVar;
        logVariance += logVar * logVar;
    }
    logMean /= n;
    logVariance /= n;
    logVariance -= logMean * logMean;

    setLocation(logMean);
    setScale(std::sqrt(logVariance));
    return true;
}
