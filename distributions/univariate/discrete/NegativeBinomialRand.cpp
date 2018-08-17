#include "NegativeBinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"
#include "PoissonRand.h"

template< typename IntType, typename T>
NegativeBinomialDistribution<IntType, T>::NegativeBinomialDistribution(T number, double probability)
{
    SetParameters(number, probability);
}

template< typename IntType, typename T>
void NegativeBinomialDistribution<IntType, T>::SetParameters(T number, double probability)
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

template< typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::P(const IntType &k) const
{
    return (k < 0) ? 0.0 : std::exp(logP(k));
}

template< typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::logP(const IntType &k) const
{
    if (k < 0)
        return -INFINITY;
    double y = std::lgammal(r + k);
    y -= RandMath::lfact(k);
    y += k * log1mProb;
    y += pdfCoef;
    return y;
}

template< typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::F(const IntType &k) const
{
    if (k < 0)
        return 0.0;
    int kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(r, kp1);
    return RandMath::ibeta(p, r, kp1, logBetaFun, logProb, log1mProb);
}

template< typename IntType, typename T>
double NegativeBinomialDistribution<IntType, T>::S(const IntType &k) const
{
    if (k < 0)
        return 0.0;
    int kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(kp1, r);
    return RandMath::ibeta(q, kp1, r, logBetaFun, log1mProb, logProb);
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::variateThroughGammaPoisson() const
{
    return PoissonRand<IntType>::Variate(GammaRV.Variate(), this->localRandGenerator);
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::variateGeometricByTable() const
{
    double U = UniformRand<double>::StandardVariate(this->localRandGenerator);
    /// handle tail by recursion
    if (U > table[tableSize - 1])
        return tableSize + variateGeometricByTable();
    /// handle the main body
    IntType x = 0;
    while (U > table[x])
        ++x;
    return x;
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::variateGeometricThroughExponential() const
{
    float X = -ExponentialRand<float>::StandardVariate(this->localRandGenerator) / log1mProb;
    return std::floor(X);
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::variateByTable() const
{
    IntType var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricByTable();
    }
    return var;
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::variateThroughExponential() const
{
    IntType var = 0;
    for (int i = 0; i < r; ++i) {
        var += variateGeometricThroughExponential();
    }
    return var;
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE)
        return variateByTable();
    return (genId == EXPONENTIAL) ? variateThroughExponential() : variateThroughGammaPoisson();
}

template< typename IntType, typename T>
void NegativeBinomialDistribution<IntType, T>::Sample(std::vector<IntType> &outputData) const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == TABLE) {
        for (IntType & var : outputData)
            var = variateByTable();
    }
    else if (genId == EXPONENTIAL) {
        for (IntType & var : outputData)
            var = variateThroughExponential();
    }
    else {
        for (IntType &var : outputData)
            var = variateThroughGammaPoisson();
    }
}

template< typename IntType, typename T>
void NegativeBinomialDistribution<IntType, T>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    GammaRV.Reseed(seed + 1);
}

template< typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Mean() const
{
    return qDivP * r;
}

template< typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Variance() const
{
    return qDivP * r / p;
}

template< typename IntType, typename T>
std::complex<double> NegativeBinomialDistribution<IntType, T>::CFImpl(double t) const
{
    std::complex<double> denominator(1.0 - q * std::cos(t), -q * std::sin(t));
    return std::pow(p / denominator, r);
}

template< typename IntType, typename T>
IntType NegativeBinomialDistribution<IntType, T>::Mode() const
{
    return (r > 1) ? std::floor((r - 1) * qDivP) : 0;
}

template< typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::Skewness() const
{
    return (1 + q) / std::sqrt(q * r);
}

template< typename IntType, typename T>
long double NegativeBinomialDistribution<IntType, T>::ExcessKurtosis() const
{
    long double kurtosis = p / qDivP;
    kurtosis += 6;
    return kurtosis / r;
}

template< typename IntType, typename T>
BetaRand<> NegativeBinomialDistribution<IntType, T>::FitProbabilityBayes(const std::vector<IntType> &sample, const BetaDistribution<> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    int n = sample.size();
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand<> posteriorDistribution(alpha + r * n, beta + this->GetSampleSum(sample));
    SetParameters(r, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template class NegativeBinomialDistribution<int, int>;
template class NegativeBinomialDistribution<long int, int>;
template class NegativeBinomialDistribution<long long int, int>;

template class NegativeBinomialDistribution<int, double>;
template class NegativeBinomialDistribution<long int, double>;
template class NegativeBinomialDistribution<long long int, double>;

template< typename IntType, typename T >
String NegativeBinomialRand<IntType, T>::Name() const
{
    if (std::is_integral_v<T>)
        return "Pascal(" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
    return "Polya(" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template< typename IntType, typename T>
constexpr char NegativeBinomialRand<IntType, T>::TOO_SMALL_VARIANCE[];

template< typename IntType, typename T >
void NegativeBinomialRand<IntType, T>::Fit(const std::vector<IntType> &sample)
{
    /// Check positivity of sample
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    /// Initial guess by method of moments
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    double mean = stats.first, variance = stats.second;
    /// Method can't be applied in the case of too small variance
    if (variance <= mean)
        throw std::invalid_argument(this->fitErrorDescription(this->NOT_APPLICABLE, this->TOO_SMALL_VARIANCE));
    double guess = mean * mean / (variance - mean);
    size_t n = sample.size();
    if (!RandMath::findRoot<double>([sample, mean, n] (double x)
    {
        double first = 0.0, second = 0.0;
        for (const IntType & var : sample) {
            first += RandMath::digamma(var + x);
            second += RandMath::trigamma(var + x);
        }
        first -= n * (RandMath::digammamLog(x) + std::log(x + mean));
        second -= n * (RandMath::trigamma(x) - mean / (x * (mean + x)));
        return DoublePair(first, second);
    }, guess))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding algorithm"));
    if (guess <= 0.0)
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, "Number should be positive, but returned value is " + this->toStringWithPrecision(guess)));
    SetParameters(guess, guess / (guess + mean));
}

template class NegativeBinomialRand<int, int>;
template class NegativeBinomialRand<long int, int>;
template class NegativeBinomialRand<long long int, int>;

template class NegativeBinomialRand<int, double>;
template class NegativeBinomialRand<long int, double>;
template class NegativeBinomialRand<long long int, double>;
