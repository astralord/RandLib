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

double LogNormalRand::logAverage(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0)
         return 0.0;
    long double logSum = 0.0L;
    for (double var : sample) {
        logSum += std::log(var);
    }
    return logSum / n;
}

double LogNormalRand::logSecondMoment(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0)
         return 0.0;
    long double logSum = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var);
        logSum += logVar * logVar;
    }
    return logSum / n;
}

bool LogNormalRand::fitLocationMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = RandMath::sampleMean(sample);
    double var = X.Variance();
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
    setLocation(logAverage(sample));
    return true;
}

bool LogNormalRand::fitScaleMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;
    double mu = X.getLocation();
    double logVariance = logSecondMoment(sample) - mu * mu;
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

bool LogNormalRand::fitLocationBayes(const std::vector<double> &sample, NormalRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double mu0 = priorDistribution.getLocation();
    double tau0 = priorDistribution.getPrecision();
    double tau = X.getPrecision();
    double numerator = n * logAverage(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    priorDistribution.setLocation(numerator / denominator);
    priorDistribution.setVariance(1.0 / denominator);
    setLocation(priorDistribution.Mean());
    return true;
}

bool LogNormalRand::fitScaleBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.getShape();
    double beta = priorDistribution.getRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.getLocation();
    double newBeta = beta + 0.5 * n * (logSecondMoment(sample) - mu * mu);
    priorDistribution.setParameters(newAlpha, 1.0 / newBeta);
    setScale(std::sqrt(priorDistribution.Mean()));
    return true;
}

bool LogNormalRand::fitLocationAndScaleBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.getShape();
    double beta = priorDistribution.getRate();
    double mu0 = priorDistribution.getLocation();
    double lambda = priorDistribution.getPrecision();
    double average = logAverage(sample), sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logSecondMoment(sample) - average * average;
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    priorDistribution.setParameters(newMu0, newLambda, newAlpha, newBeta);
    double2d mean = priorDistribution.Mean();
    setLocation(mean.x);
    setScale(std::sqrt(mean.y));
    return true;
}
