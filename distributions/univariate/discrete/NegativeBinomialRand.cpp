#include "NegativeBinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"
#include "PoissonRand.h"

template< typename T >
NegativeBinomialDistribution<T>::NegativeBinomialDistribution(T number, double probability)
{
    SetParameters(number, probability);
}

template< typename T >
void NegativeBinomialDistribution<T>::SetParameters(T number, double probability)
{
    r = (number > 0) ? number : 1;
    p = (probability > 1.0 || probability < 0.0) ? 0.5 : probability;
    q = 1.0 - p;
    logProb = std::log(p);
    log1mProb = std::log1p(-p);
    GammaRV.SetParameters(r, p / q);
    qDivP = GammaRV.GetScale();
    pdfCoef = r * logProb;
    pdfCoef -= GammaRV.GetLogGammaShape();

    if (GetIdOfUsedGenerator() == TABLE) {
        /// table method
        table[0] = p;
        double prod = p;
        for (int i = 1; i < tableSize; ++i)
        {
            prod *= q;
            table[i] = table[i - 1] + prod;
        }
    }
}

template< typename T >
double NegativeBinomialDistribution<T>::P(const int & k) const
{
    return (k < 0) ? 0.0 : std::exp(logP(k));
}

template< typename T >
double NegativeBinomialDistribution<T>::logP(const int & k) const
{
    if (k < 0)
        return -INFINITY;
    double y = std::lgamma(r + k);
    y -= RandMath::lfact(k);
    y += k * log1mProb;
    y += pdfCoef;
    return y;
}

template< typename T >
double NegativeBinomialDistribution<T>::F(const int & k) const
{
    if (k < 0)
        return 0.0;
    int kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(r, kp1);
    return RandMath::ibeta(p, r, kp1, logBetaFun, logProb, log1mProb);
}

template< typename T >
double NegativeBinomialDistribution<T>::S(const int & k) const
{
    if (k < 0)
        return 0.0;
    int kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(kp1, r);
    return RandMath::ibeta(q, kp1, r, logBetaFun, log1mProb, logProb);
}

template< >
NegativeBinomialRand<int>::GENERATOR_ID NegativeBinomialDistribution<int>::GetIdOfUsedGenerator() const
{
    /// if r is small, we use two different generators for two different cases:
    /// if p < 0.08 then the tail is too heavy
    /// (probability to be in main body is less than 0.75),
    /// then we return highest integer less than variate from exponential distribution
    /// otherwise we choose table method
    if (r < 10)
        return (p < 0.08) ? EXPONENTIAL : TABLE;
    return GAMMA_POISSON;
}

template< >
NegativeBinomialRand<double>::GENERATOR_ID NegativeBinomialDistribution<double>::GetIdOfUsedGenerator() const
{
    return GAMMA_POISSON;
}

template< typename T >
int NegativeBinomialDistribution<T>::variateThroughGammaPoisson() const
{
    return PoissonRand::Variate(GammaRV.Variate());
}

template< >
int NegativeBinomialDistribution<int>::variateGeometricByTable() const
{
    double U = UniformRand::StandardVariate();
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateGeometricByTable();
    /// handle the main body
    int x = 0;
    while (U > table[x])
        ++x;
    return x;
}

template< >
int NegativeBinomialDistribution<int>::variateGeometricThroughExponential() const
{
    return std::floor(-ExponentialRand::StandardVariate() / log1mProb);
}

template< >
int NegativeBinomialDistribution<int>::variateByTable() const
{
    double var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricByTable();
    }
    return var;
}

template< >
int NegativeBinomialDistribution<int>::variateThroughExponential() const
{
    double var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricThroughExponential();
    }
    return var;
}

template< >
int NegativeBinomialDistribution<double>::Variate() const
{
    return variateThroughGammaPoisson();
}

template< >
int NegativeBinomialDistribution<int>::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE)
        return variateByTable();
    return (genId == EXPONENTIAL) ? variateThroughExponential() : variateThroughGammaPoisson();
}

template< >
void NegativeBinomialDistribution<double>::Sample(std::vector<int> &outputData) const
{
    for (int &var : outputData)
        var = variateThroughGammaPoisson();
}

template< >
void NegativeBinomialDistribution<int>::Sample(std::vector<int> &outputData) const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE) {
        for (int & var : outputData)
            var = variateByTable();
    }
    else if (genId == EXPONENTIAL) {
        for (int & var : outputData)
            var = variateThroughExponential();
    }
    else {
        for (int &var : outputData)
            var = variateThroughGammaPoisson();
    }
}

template< typename T >
double NegativeBinomialDistribution<T>::Mean() const
{
    return qDivP * r;
}

template< typename T >
double NegativeBinomialDistribution<T>::Variance() const
{
    return qDivP * r / p;
}

template< typename T >
std::complex<double> NegativeBinomialDistribution<T>::CFImpl(double t) const
{
    std::complex<double> denominator(1.0 - q * std::cos(t), -q * std::sin(t));
    return std::pow(p / denominator, r);
}

template< typename T >
int NegativeBinomialDistribution<T>::Mode() const
{
    return (r > 1) ? std::floor((r - 1) * qDivP) : 0;
}

template< typename T >
double NegativeBinomialDistribution<T>::Skewness() const
{
    return (1 + q) / std::sqrt(q * r);
}

template< typename T >
double NegativeBinomialDistribution<T>::ExcessKurtosis() const
{
    double kurtosis = p / qDivP;
    kurtosis += 6;
    return kurtosis / r;
}

template< typename T >
void NegativeBinomialDistribution<T>::FitProbabilityMM(const std::vector<int> &sample)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    double mean = sampleMean(sample);
    SetParameters(r, r / (r + mean));
}

template< typename T >
BetaRand NegativeBinomialDistribution<T>::FitProbabilityBayes(const std::vector<int> &sample, const BetaDistribution &priorDistribution)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    int n = sample.size();
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand posteriorDistribution(alpha + r * n, beta + sampleSum(sample));
    SetParameters(r, posteriorDistribution.Mean());
    return posteriorDistribution;
}


template< typename T >
std::string NegativeBinomialRand<T>::Name() const
{
    return "Negative Binomial(" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template< >
void NegativeBinomialRand<double>::FitNumberMM(const std::vector<int> &sample)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    double mean = sampleMean(sample);
    SetParameters(mean * p / q, p);
}

template< typename T >
constexpr char NegativeBinomialRand<T>::TOO_SMALL_VARIANCE[];

template< >
void NegativeBinomialRand<double>::FitNumberAndProbabilityMM(const std::vector<int> &sample)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    size_t n = sample.size();
    if (n <= 1)
        throw std::invalid_argument(fitError(TOO_FEW_ELEMENTS, "There should be at least 2 elements"));
    double mean = sampleMean(sample);
    double variance = sampleVariance(sample, mean);
    if (variance <= mean)
        throw std::invalid_argument(fitError(NOT_APPLICABLE, TOO_SMALL_VARIANCE));
    SetParameters(mean * mean / (variance - mean), mean / variance);
}

template< >
void NegativeBinomialRand<double>::FitNumberAndProbabilityMLE(const std::vector<int> &sample)
{
    /// Check positivity of sample
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    /// Initial guess by method of moments
    double mean = sampleMean(sample), variance = sampleVariance(sample, mean);
    /// Method can't be applied in the case of too small variance
    if (variance <= mean)
        throw std::invalid_argument(fitError(NOT_APPLICABLE, TOO_SMALL_VARIANCE));
    double guess = mean * mean / (variance - mean);
    size_t n = sample.size();
    if (!RandMath::findRoot([sample, mean, n] (double x)
    {
        double first = 0.0, second = 0.0;
        for (const double & var : sample) {
            first += RandMath::digamma(var + x);
            second += RandMath::trigamma(var + x);
        }
        first -= n * (RandMath::digamma(x) - std::log(x / (x + mean)));
        second -= n * (RandMath::trigamma(x) - mean / (x * (mean + x)));
        return DoublePair(first, second);
    }, guess))
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding algorithm"));
    if (guess <= 0.0)
        throw std::runtime_error(fitError(WRONG_RETURN, "Number should be positive, but returned value is " + toStringWithPrecision(guess)));
    SetParameters(guess, guess / (guess + mean));
}

template class NegativeBinomialDistribution<int>;
template class NegativeBinomialDistribution<double>;
