#include "LogNormalRand.h"

LogNormalRand::LogNormalRand(double location, double squaredScale)
{
    SetLocation(location);
    SetScale(squaredScale > 0.0 ? std::sqrt(squaredScale) : 1.0);
}

std::string LogNormalRand::Name() const
{
    return "Log-Normal(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void LogNormalRand::SetLocation(double location)
{
    X.SetLocation(location);
    expMu = std::exp(X.Mean());
}

void LogNormalRand::SetScale(double scale)
{
    X.SetScale(scale);
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

double LogNormalRand::StandardVariate()
{
    return std::exp(NormalRand::StandardVariate());
}

double LogNormalRand::Variate(double location, double scale)
{
    return std::exp(NormalRand::Variate(location, scale));
}

double LogNormalRand::Variate() const
{
    return std::exp(X.Variate());
}

double LogNormalRand::Mean() const
{
    return expMu * std::sqrt(expVar);
}

double LogNormalRand::Variance() const
{
    return (expVar - 1) * expMu * expMu * expVar;
}

double LogNormalRand::quantileImpl(double p) const
{
    return std::exp(X.Quantile(p));
}

double LogNormalRand::quantileImpl1m(double p) const
{
    return std::exp(X.Quantile1m(p));
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

bool LogNormalRand::FitLocationMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = sampleMean(sample);
    double var = X.Variance();
    SetLocation(std::log(average) - 0.5 * var);
    return true;
}

bool LogNormalRand::FitScaleMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = sampleMean(sample);
    double mu = X.GetLocation();
    double aux = std::log(average) - mu;
    SetScale(std::sqrt(aux + aux));
    return true;
}

bool LogNormalRand::FitMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double average = sampleMean(sample);
    double secondMoment = rawMoment(sample, 2);
    double averageSq = average * average;
    SetLocation(0.5 * std::log(averageSq * averageSq / secondMoment));
    SetScale(std::sqrt(std::log(secondMoment / averageSq)));
    return true;
}

bool LogNormalRand::FitLocationMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;
    SetLocation(logAverage(sample));
    return true;
}

bool LogNormalRand::FitScaleMLE(const std::vector<double> &sample)
{
    size_t n = sample.size();
    if (n == 0 || !checkValidity(sample))
        return false;
    double mu = X.GetLocation();
    double logVariance = logSecondMoment(sample) - mu * mu;
    SetScale(std::sqrt(logVariance));
    return true;
}

bool LogNormalRand::FitMLE(const std::vector<double> &sample)
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

    SetLocation(logMean);
    SetScale(std::sqrt(logVariance));
    return true;
}

bool LogNormalRand::FitLocationBayes(const std::vector<double> &sample, NormalRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = X.GetPrecision();
    double numerator = n * logAverage(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    priorDistribution.SetLocation(numerator / denominator);
    priorDistribution.SetVariance(1.0 / denominator);
    SetLocation(priorDistribution.Mean());
    return true;
}

bool LogNormalRand::FitScaleBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.GetLocation();
    double newBeta = beta + 0.5 * n * (logSecondMoment(sample) - mu * mu);
    priorDistribution.SetParameters(newAlpha, 1.0 / newBeta);
    SetScale(std::sqrt(priorDistribution.Mean()));
    return true;
}

bool LogNormalRand::FitBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    double average = logAverage(sample), sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logSecondMoment(sample) - average * average;
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    priorDistribution.SetParameters(newMu0, newLambda, newAlpha, newBeta);
    DoublePair mean = priorDistribution.Mean();
    SetLocation(mean.first);
    SetScale(std::sqrt(mean.second));
    return true;
}
