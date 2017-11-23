#include "LogNormalRand.h"

template < typename RealType >
LogNormalRand<RealType>::LogNormalRand(double location, double squaredScale)
{
    SetLocation(location);
    SetScale(squaredScale > 0.0 ? std::sqrt(squaredScale) : 1.0);
}

template < typename RealType >
String LogNormalRand<RealType>::Name() const
{
    return "Log-Normal(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void LogNormalRand<RealType>::SetLocation(double location)
{
    X.SetLocation(location);
    expMu = std::exp(X.GetLocation());
}

template < typename RealType >
void LogNormalRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Log-Normal distribution: scale should be positive");
    X.SetScale(scale);
    expHalfSigmaSq = std::exp(0.5 * X.Variance());
}

template < typename RealType >
double LogNormalRand<RealType>::f(const RealType & x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

template < typename RealType >
double LogNormalRand<RealType>::logf(const RealType & x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = X.logf(logX);
    return y - logX;
}

template < typename RealType >
double LogNormalRand<RealType>::F(const RealType & x) const
{
    return (x > 0.0) ? X.F(std::log(x)) : 0.0;
}

template < typename RealType >
double LogNormalRand<RealType>::S(const RealType & x) const
{
    return (x > 0.0) ? X.S(std::log(x)) : 1.0;
}

template < typename RealType >
RealType LogNormalRand<RealType>::Variate() const
{
    return std::exp(X.Variate());
}

template < typename RealType >
RealType LogNormalRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    return std::exp(NormalRand<RealType>::StandardVariate(randGenerator));
}

template < typename RealType >
void LogNormalRand<RealType>::Reseed(unsigned long seed) const
{
    X.Reseed(seed);
}

template < typename RealType >
long double LogNormalRand<RealType>::Mean() const
{
    return expMu * expHalfSigmaSq;
}

template < typename RealType >
long double LogNormalRand<RealType>::Variance() const
{
    double y = expMu * expHalfSigmaSq;
    return y * y * std::expm1l(X.Variance());
}

template < typename RealType >
RealType LogNormalRand<RealType>::quantileImpl(double p) const
{
    return std::exp(X.Quantile(p));
}

template < typename RealType >
RealType LogNormalRand<RealType>::quantileImpl1m(double p) const
{
    return std::exp(X.Quantile1m(p));
}

template < typename RealType >
RealType LogNormalRand<RealType>::Median() const
{
    return expMu;
}

template < typename RealType >
RealType LogNormalRand<RealType>::Mode() const
{
    return expMu / (expHalfSigmaSq * expHalfSigmaSq);
}

template < typename RealType >
long double LogNormalRand<RealType>::Skewness() const
{
    double y = std::expm1l(X.Variance());
    return (expHalfSigmaSq * expHalfSigmaSq + 2) * std::sqrt(y);
}

template < typename RealType >
long double LogNormalRand<RealType>::ExcessKurtosis() const
{
    double temp = expHalfSigmaSq * expHalfSigmaSq;
    double c = temp * temp;
    double b = c * temp;
    double a = c * c;
    return a + 2 * b + 3 * c - 6;
}

template < typename RealType >
void LogNormalRand<RealType>::FitLocation(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    this->SetLocation(this->GetSampleLogMean(sample));
}

template < typename RealType >
void LogNormalRand<RealType>::FitScale(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    double mu = X.GetLocation();
    this->SetScale(std::sqrt(this->GetSampleLogVariance(sample, mu)));
}

template < typename RealType >
void LogNormalRand<RealType>::Fit(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
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

template < typename RealType >
NormalRand<RealType> LogNormalRand<RealType>::FitLocationBayes(const std::vector<RealType> &sample, const NormalRand<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = X.GetPrecision();
    double numerator = n * this->GetSampleLogMean(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    NormalRand<RealType> posteriorDistribution(numerator / denominator, 1.0 / denominator);
    SetLocation(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template < typename RealType >
InverseGammaRand<RealType> LogNormalRand<RealType>::FitScaleBayes(const std::vector<RealType> &sample, const InverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double mu = X.GetLocation();
    double newBeta = beta + 0.5 * n * this->GetSampleLogVariance(sample, mu);
    InverseGammaRand<RealType> posteriorDistribution(newAlpha, newBeta);
    SetScale(std::sqrt(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean()));
    return posteriorDistribution;
}

template < typename RealType >
NormalInverseGammaRand<RealType> LogNormalRand<RealType>::FitBayes(const std::vector<RealType> &sample, const NormalInverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    DoublePair logStats = this->GetSampleLogMeanAndVariance(sample);
    double average = logStats.first, sum = n * average;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = logStats.second;
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    NormalInverseGammaRand<RealType> posteriorDistribution(newMu0, newLambda, newAlpha, newBeta);
    DoublePair newParams = MAP ? static_cast<DoublePair>(posteriorDistribution.Mode()) : static_cast<DoublePair>(posteriorDistribution.Mean());
    SetLocation(newParams.first);
    SetScale(std::sqrt(newParams.second));
    return posteriorDistribution;
}


template class LogNormalRand<float>;
template class LogNormalRand<double>;
template class LogNormalRand<long double>;
