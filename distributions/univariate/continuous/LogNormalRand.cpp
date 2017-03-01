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
    long double logSum = 0.0L;
    for (double var : sample) {
        logSum += std::log(var);
    }
    return logSum / n;
}

double LogNormalRand::logVariance(const std::vector<double> &sample, double mu)
{
    size_t n = sample.size();
    long double logSum = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var) - mu;
        logSum += logVar * logVar;
    }
    return logSum / n;
}

void LogNormalRand::FitLocationMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double average = sampleMean(sample);
    double var = X.Variance();
    SetLocation(std::log(average) - 0.5 * var);
}

void LogNormalRand::FitScaleMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double average = sampleMean(sample);
    double mu = X.GetLocation();
    double aux = std::log(average) - mu;
    SetScale(std::sqrt(aux + aux));
}

void LogNormalRand::FitMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double average = sampleMean(sample);
    double secondMoment = rawMoment(sample, 2);
    double averageSq = average * average;
    SetLocation(0.5 * std::log(averageSq * averageSq / secondMoment));
    SetScale(std::sqrt(std::log(secondMoment / averageSq)));
}

void LogNormalRand::FitLocationMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    SetLocation(logAverage(sample));
}

void LogNormalRand::FitScaleMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double mu = X.GetLocation();
    SetScale(std::sqrt(logVariance(sample, mu)));
}

void LogNormalRand::FitMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
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
}

NormalRand LogNormalRand::FitLocationBayes(const std::vector<double> &sample, const NormalRand &priorDistribution)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = X.GetPrecision();
    double numerator = n * logAverage(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    NormalRand posteriorDistribution(numerator / denominator, 1.0 / denominator);
    SetLocation(posteriorDistribution.Mean());
    return posteriorDistribution;
}

InverseGammaRand LogNormalRand::FitScaleBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.GetLocation();
    double newBeta = beta + 0.5 * n * logVariance(sample, mu);
    InverseGammaRand posteriorDistribution(newAlpha, newBeta);
    SetScale(std::sqrt(posteriorDistribution.Mean()));
    return posteriorDistribution;
}

NormalInverseGammaRand LogNormalRand::FitBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    double average = logAverage(sample), sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logVariance(sample, average);
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    NormalInverseGammaRand posteriorDistribution(newMu0, newLambda, newAlpha, newBeta);
    DoublePair mean = posteriorDistribution.Mean();
    SetLocation(mean.first);
    SetScale(std::sqrt(mean.second));
    return posteriorDistribution;
}
