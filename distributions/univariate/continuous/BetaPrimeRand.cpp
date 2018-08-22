#include "BetaPrimeRand.h"

template < typename RealType >
BetaPrimeRand<RealType>::BetaPrimeRand(double shape1, double shape2)
{
    SetShapes(shape1, shape2);
}

template < typename RealType >
String BetaPrimeRand<RealType>::Name() const
{
    return "Beta Prime(" + this->toStringWithPrecision(GetAlpha()) + ", " + this->toStringWithPrecision(GetBeta()) + ")";
}

template < typename RealType >
void BetaPrimeRand<RealType>::SetShapes(double shape1, double shape2)
{
    if (shape1 <= 0 || shape2 <= 0)
        throw std::invalid_argument("Beta-prime distribution: shapes should be positive");
    B.SetShapes(shape1, shape2);
    alpha = B.GetAlpha();
    beta = B.GetBeta();
}

template < typename RealType >
double BetaPrimeRand<RealType>::f(const RealType & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (alpha == 1.0)
            return 1.0 / GetBetaFunction();
        return (alpha > 1) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

template < typename RealType >
double BetaPrimeRand<RealType>::logf(const RealType &x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (alpha == 1.0)
            return -GetLogBetaFunction();
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    double y = (alpha - 1) * std::log(x);
    y -= (alpha + beta) * std::log1pl(x);
    return y - GetLogBetaFunction();
}

template < typename RealType >
double BetaPrimeRand<RealType>::F(const RealType & x) const
{
    return (x > 0) ? B.F(x / (1.0 + x)) : 0;
}

template < typename RealType >
double BetaPrimeRand<RealType>::S(const RealType &x) const
{
    return (x > 0) ? B.S(x / (1.0 + x)) : 1;
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::fromBetaVariate(const RealType & betaVar) const
{
    if (betaVar > 1e-5)
        return betaVar / (1.0 - betaVar);
    RealType logVar = std::log(betaVar), log1mVar = std::log1p(-betaVar);
    return std::exp(logVar - log1mVar);
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::Variate() const
{
    double x = B.Variate();
    return fromBetaVariate(x);
}

template < typename RealType >
void BetaPrimeRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    B.Sample(outputData);
    for (RealType &var : outputData)
        var = fromBetaVariate(var);
}

template < typename RealType >
void BetaPrimeRand<RealType>::Reseed(unsigned long seed) const
{
    B.Reseed(seed);
}

template < typename RealType >
long double BetaPrimeRand<RealType>::Mean() const
{
    return (beta > 1) ? alpha / (beta - 1) : INFINITY;
}

template < typename RealType >
long double BetaPrimeRand<RealType>::Variance() const
{
    if (beta <= 2)
        return INFINITY;
    double betam1 = beta - 1;
    double numerator = alpha * (alpha + betam1);
    double denominator = (betam1 - 1) * betam1 * betam1;
    return numerator / denominator;
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::Median() const
{
    return (alpha == beta) ? 1.0 : quantileImpl(0.5);
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::Mode() const
{
    return (alpha < 1) ? 0 : (alpha - 1) / (beta + 1);
}

template < typename RealType >
long double BetaPrimeRand<RealType>::Skewness() const
{
    if (beta <= 3)
        return INFINITY;
    long double aux = alpha + beta - 1;
    long double skewness = (beta - 2) / (alpha * aux);
    skewness = std::sqrt(skewness);
    aux += alpha;
    aux += aux;
    return aux * skewness / (beta - 3);
}

template < typename RealType >
long double BetaPrimeRand<RealType>::ExcessKurtosis() const
{
    if (beta <= 4)
        return INFINITY;
    long double betam1 = beta - 1;
    long double numerator = betam1 * betam1 * (beta - 2) / (alpha * (alpha + betam1));
    numerator += 5 * beta - 11;
    long double denominator = (beta - 3) * (beta - 4);
    return 6 * numerator / denominator;
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::quantileImpl(double p) const
{
    double x = B.Quantile(p);
    return x / (1.0 - x);
}

template < typename RealType >
RealType BetaPrimeRand<RealType>::quantileImpl1m(double p) const
{
    double x = B.Quantile1m(p);
    return x / (1.0 - x);
}

template < typename RealType >
std::complex<double> BetaPrimeRand<RealType>::CFImpl(double t) const
{
    /// if no singularity - simple numeric integration
    if (alpha >= 1)
        return UnivariateDistribution<RealType>::CFImpl(t);

    double re = this->ExpectedValue([this, t] (double x)
    {
        if (x == 0.0)
            return 0.0;
        return std::cos(t * x) - 1.0;
    }, 0.0, INFINITY) + 1.0;

    double im = this->ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, 0.0, INFINITY);
    return std::complex<double>(re, im);
}

template < typename RealType >
void BetaPrimeRand<RealType>::FitAlpha(const std::vector<RealType> &sample)
{
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    long double lnG = this->GetSampleLogMean(sample) - B.GetSampleLog1pMeanNorm(sample);
    long double mean = 0.5;
    if (beta != 1.0) {
        mean = 0.0;
        for (const double & var : sample)
            mean += var / (1.0 + var);
        mean /= sample.size();
    }
    B.FitAlpha(lnG, mean);
    SetShapes(B.GetAlpha(), beta);
}

template < typename RealType >
void BetaPrimeRand<RealType>::FitBeta(const std::vector<RealType> &sample)
{
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    long double lnG1m = -B.GetSampleLog1pMeanNorm(sample);
    long double mean = 0.5;
    if (alpha != 1.0) {
        mean = 0.0;
        for (const double & var : sample)
            mean += var / (1.0 + var);
        mean /= sample.size();
    }
    B.FitBeta(lnG1m, mean);
    SetShapes(alpha, B.GetBeta());
}

template < typename RealType >
void BetaPrimeRand<RealType>::Fit(const std::vector<RealType> &sample)
{
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    long double lnG1m = -B.GetSampleLog1pMeanNorm(sample);
    long double lnG = this->GetSampleLogMean(sample) + lnG1m;
    long double m = 0.0, v = 0.0;
    int n = sample.size();
    for (int i = 0; i < n; ++i) {
        double x = sample[i] / (1.0 + sample[i]);
        double diff = x - m;
        m += diff / (i + 1);
        v += diff * (x - m);
    }
    B.FitShapes(lnG, lnG1m, m, v / n);
    SetShapes(B.GetAlpha(), B.GetBeta());
}

template class BetaPrimeRand<float>;
template class BetaPrimeRand<double>;
template class BetaPrimeRand<long double>;
