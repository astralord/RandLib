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
    expMu = std::exp(X.GetLocation());
}

void LogNormalRand::SetScale(double scale)
{
    X.SetScale(scale);
    expHalfSigmaSq = std::exp(0.5 * X.Variance());
}

double LogNormalRand::f(double x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

double LogNormalRand::logf(double x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = X.logf(logX);
    return y - logX;
}

double LogNormalRand::F(double x) const
{
    return (x > 0.0) ? X.F(std::log(x)) : 0.0;
}

double LogNormalRand::S(double x) const
{
    return (x > 0.0) ? X.S(std::log(x)) : 1.0;
}

double LogNormalRand::StandardVariate()
{
    return std::exp(NormalRand::StandardVariate());
}

double LogNormalRand::Variate(double location, double scale)
{
    return std::exp(location + scale * NormalRand::StandardVariate());
}

double LogNormalRand::Variate() const
{
    return std::exp(X.Variate());
}

double LogNormalRand::Mean() const
{
    return expMu * expHalfSigmaSq;
}

double LogNormalRand::Variance() const
{
    double y = expMu * expHalfSigmaSq;
    return y * y * std::expm1(X.Variance());
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
    return expMu / (expHalfSigmaSq * expHalfSigmaSq);
}

double LogNormalRand::Skewness() const
{
    double y = std::expm1(X.Variance());
    return (expHalfSigmaSq * expHalfSigmaSq + 2) * std::sqrt(y);
}

double LogNormalRand::ExcessKurtosis() const
{
    double a = std::pow(expHalfSigmaSq, 8);
    double b = 2 * std::pow(expHalfSigmaSq, 6);
    double c = 3 * std::pow(expHalfSigmaSq, 4);
    return a + b + c - 6;
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
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        return false;
    double average = sampleMean(sample);
    double var = X.Variance();
    SetLocation(std::log(average) - 0.5 * var);
    return true;
}

bool LogNormalRand::FitScaleMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        return false;
    double average = sampleMean(sample);
    double mu = X.GetLocation();
    double aux = std::log(average) - mu;
    SetScale(std::sqrt(aux + aux));
    return true;
}

bool LogNormalRand::FitMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
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
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        return false;
    SetLocation(logAverage(sample));
    return true;
}

bool LogNormalRand::FitScaleMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        return false;
    double mu = X.GetLocation();
    double logVariance = logSecondMoment(sample) - mu * mu;
    SetScale(std::sqrt(logVariance));
    return true;
}

bool LogNormalRand::FitMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        return false;
    size_t n = sample.size();
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
