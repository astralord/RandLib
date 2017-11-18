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
    if (r <= 0.0)
        throw std::invalid_argument("Negative-Binomial distribution: number parameter should be positive");
    if (probability <= 0.0 || probability >= 1.0)
        throw std::invalid_argument("Negative-Binomial distribution: probability parameter should in interval (0, 1)");
    r = (number > 0) ? number : 1;
    p = probability;
    q = 1.0 - p;
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
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
    double y = std::lgammal(r + k);
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

template< typename T >
int NegativeBinomialDistribution<T>::variateThroughGammaPoisson() const
{
    return PoissonRand::Variate(GammaRV.Variate(), localRandGenerator);
}

template< typename T >
int NegativeBinomialDistribution<T>::variateGeometricByTable() const
{
    double U = UniformRand::StandardVariate(localRandGenerator);
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateGeometricByTable();
    /// handle the main body
    int x = 0;
    while (U > table[x])
        ++x;
    return x;
}

template< typename T >
int NegativeBinomialDistribution<T>::variateGeometricThroughExponential() const
{
    double X = std::floor(-ExponentialRand::StandardVariate(localRandGenerator) / log1mProb);
    return X < INT_MAX ? X : INT_MAX - 1;
}

template< typename T >
int NegativeBinomialDistribution<T>::variateByTable() const
{
    double var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricByTable();
    }
    return var;
}

template< typename T >
int NegativeBinomialDistribution<T>::variateThroughExponential() const
{
    double var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricThroughExponential();
    }
    return var;
}

template< typename T >
int NegativeBinomialDistribution<T>::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE)
        return variateByTable();
    return (genId == EXPONENTIAL) ? variateThroughExponential() : variateThroughGammaPoisson();
}

template< typename T >
void NegativeBinomialDistribution<T>::Sample(std::vector<int> &outputData) const
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
void NegativeBinomialDistribution<T>::Reseed(unsigned long seed) const
{
    localRandGenerator.Reseed(seed);
    GammaRV.Reseed(seed + 1);
}

template< typename T >
long double NegativeBinomialDistribution<T>::Mean() const
{
    return qDivP * r;
}

template< typename T >
long double NegativeBinomialDistribution<T>::Variance() const
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
long double NegativeBinomialDistribution<T>::Skewness() const
{
    return (1 + q) / std::sqrt(q * r);
}

template< typename T >
long double NegativeBinomialDistribution<T>::ExcessKurtosis() const
{
    long double kurtosis = p / qDivP;
    kurtosis += 6;
    return kurtosis / r;
}

template< typename T >
BetaRand NegativeBinomialDistribution<T>::FitProbabilityBayes(const std::vector<int> &sample, const BetaDistribution &priorDistribution, bool MAP)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    int n = sample.size();
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand posteriorDistribution(alpha + r * n, beta + GetSampleSum(sample));
    SetParameters(r, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}


template< >
String NegativeBinomialRand<int>::Name() const
{
    return "Pascal" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template< >
String NegativeBinomialRand<double>::Name() const
{
    return "Polya(" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template< typename T >
constexpr char NegativeBinomialRand<T>::TOO_SMALL_VARIANCE[];

template< >
void NegativeBinomialRand<double>::Fit(const std::vector<int> &sample)
{
    /// Check positivity of sample
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    /// Initial guess by method of moments
    DoublePair stats = GetSampleMeanAndVariance(sample);
    double mean = stats.first, variance = stats.second;
    /// Method can't be applied in the case of too small variance
    if (variance <= mean)
        throw std::invalid_argument(fitErrorDescription(NOT_APPLICABLE, TOO_SMALL_VARIANCE));
    double guess = mean * mean / (variance - mean);
    size_t n = sample.size();
    if (!RandMath::findRoot<double>([sample, mean, n] (double x)
    {
        double first = 0.0, second = 0.0;
        for (const double & var : sample) {
            first += RandMath::digamma(var + x);
            second += RandMath::trigamma(var + x);
        }
        first -= n * (RandMath::digammamLog(x) + std::log(x + mean));
        second -= n * (RandMath::trigamma(x) - mean / (x * (mean + x)));
        return DoublePair(first, second);
    }, guess))
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Error in root-finding algorithm"));
    if (guess <= 0.0)
        throw std::runtime_error(fitErrorDescription(WRONG_RETURN, "Number should be positive, but returned value is " + toStringWithPrecision(guess)));
    SetParameters(guess, guess / (guess + mean));
}

template class NegativeBinomialDistribution<int>;
template class NegativeBinomialDistribution<double>;
