#include "LogNormalRand.h"

LogNormalRand::LogNormalRand(double location, double squaredScale)
{
    SetLocation(location);
    SetScale(squaredScale > 0.0 ? std::sqrt(squaredScale) : 1.0);
}

String LogNormalRand::Name() const
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
    if (scale <= 0.0)
        throw std::invalid_argument("Log-Normal distribution: scale should be positive");
    X.SetScale(scale);
    expHalfSigmaSq = std::exp(0.5 * X.Variance());
}

double LogNormalRand::f(const double & x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

double LogNormalRand::logf(const double & x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = X.logf(logX);
    return y - logX;
}

double LogNormalRand::F(const double & x) const
{
    return (x > 0.0) ? X.F(std::log(x)) : 0.0;
}

double LogNormalRand::S(const double & x) const
{
    return (x > 0.0) ? X.S(std::log(x)) : 1.0;
}

double LogNormalRand::Variate() const
{
    return std::exp(X.Variate());
}

double LogNormalRand::StandardVariate(RandGenerator &randGenerator)
{
    return std::exp(NormalRand::StandardVariate(randGenerator));
}

void LogNormalRand::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
}

long double LogNormalRand::Mean() const
{
    return expMu * expHalfSigmaSq;
}

long double LogNormalRand::Variance() const
{
    double y = expMu * expHalfSigmaSq;
    return y * y * std::expm1l(X.Variance());
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

long double LogNormalRand::Skewness() const
{
    double y = std::expm1l(X.Variance());
    return (expHalfSigmaSq * expHalfSigmaSq + 2) * std::sqrt(y);
}

long double LogNormalRand::ExcessKurtosis() const
{
    double temp = expHalfSigmaSq * expHalfSigmaSq;
    double c = temp * temp;
    double b = c * temp;
    double a = c * c;
    return a + 2 * b + 3 * c - 6;
}

void LogNormalRand::FitLocation(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    SetLocation(GetSampleLogMean(sample));
}

void LogNormalRand::FitScale(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double mu = X.GetLocation();
    SetScale(std::sqrt(GetSampleLogVariance(sample, mu)));
}

void LogNormalRand::Fit(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    long double logMean = 0.0L;
    long double logSqDev = 0.0L;
    for (double var : sample) {
        double logVar = std::log(var);
        logMean += logVar;
        logSqDev += logVar * logVar;
    }
    logMean /= n;
    logSqDev /= n;
    logSqDev -= logMean * logMean;

    SetLocation(logMean);
    SetScale(std::sqrt(logSqDev));
}

NormalRand LogNormalRand::FitLocationBayes(const std::vector<double> &sample, const NormalRand &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = X.GetPrecision();
    double numerator = n * GetSampleLogMean(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    NormalRand posteriorDistribution(numerator / denominator, 1.0 / denominator);
    SetLocation(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

InverseGammaRand LogNormalRand::FitScaleBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.GetLocation();
    double newBeta = beta + 0.5 * n * GetSampleLogVariance(sample, mu);
    InverseGammaRand posteriorDistribution(newAlpha, newBeta);
    SetScale(std::sqrt(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean()));
    return posteriorDistribution;
}

NormalInverseGammaRand LogNormalRand::FitBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    DoublePair logStats = GetSampleLogMeanAndVariance(sample);
    double average = logStats.first, sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logStats.second;
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    NormalInverseGammaRand posteriorDistribution(newMu0, newLambda, newAlpha, newBeta);
    DoublePair newParams = MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean();
    SetLocation(newParams.first);
    SetScale(std::sqrt(newParams.second));
    return posteriorDistribution;
}
