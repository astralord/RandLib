#include "BinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/NormalRand.h"
#include "../continuous/ExponentialRand.h"
#include "BernoulliRand.h"

template< typename IntType >
BinomialDistribution<IntType>::BinomialDistribution(IntType number, double probability)
{
    SetParameters(number, probability);
}

template< typename IntType >
void BinomialDistribution<IntType>::SetGeneratorConstants()
{
    minpq = std::min(p, q);
    npFloor = std::floor(n * minpq);
    pFloor = npFloor / n;
    pRes = RandMath::areClose(npFloor, n * minpq) ? 0.0 : minpq - pFloor;

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == BERNOULLI_SUM)
        return;
    else if (genId == WAITING) {
        G.SetProbability(minpq);
        return;
    }

    nqFloor = n - npFloor;
    double qFloor = 1.0 - pFloor;
    if (pRes > 0)
        G.SetProbability(pRes / qFloor);

    /// Set deltas
    double npq = npFloor * qFloor;
    double coef = 128.0 * n / M_PI;
    delta1 = coef * pFloor / (81.0 * qFloor);
    delta1 = npq * std::log(delta1);
    if (delta1 > 1.0)
        delta1 = std::sqrt(delta1);
    else
        delta1 = 1.0;
    delta2 = coef * qFloor / pFloor;
    delta2 = npq * std::log(delta2);
    if (delta2 > 1.0)
        delta2 = std::sqrt(delta2);
    else
        delta2 = 1.0;

    /// Set sigmas and c
    double npqSqrt = std::sqrt(npq);
    sigma1 = npqSqrt * (1.0 + 0.25 * delta1 / npFloor);
    sigma2 = npqSqrt * (1.0 + 0.25 * delta2 / nqFloor);
    c = 2.0 * delta1 / npFloor;

    /// Set a's
    a1 = 0.5 * std::exp(c) * sigma1 * M_SQRT2PI;

    a2 = 0.5 * sigma2 * M_SQRT2PI;
    a2 += a1;

    coefa3 = 0.5 * delta1 / (sigma1 * sigma1);
    a3 = 1.0 / nqFloor - coefa3;
    a3 *= delta1;
    a3 = std::exp(a3);
    a3 /= coefa3;
    a3 += a2;

    coefa4 = 0.5 * delta2 / (sigma2 * sigma2);
    a4 = std::exp(-delta2 * coefa4) / coefa4;
    a4 += a3;

    logPFloor = std::log(pFloor);
    logQFloor = (pFloor == qFloor) ? logPFloor : std::log1pl(-pFloor);

    logPnpInv = logProbFloor(npFloor);
}

template< typename IntType >
void BinomialDistribution<IntType>::SetParameters(IntType number, double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Binomial distribution: probability parameter should in interval [0, 1]");
    if (number <= 0)
        throw std::invalid_argument("Binomial distribution: number should be positive");
    n = number;
    p = probability;
    q = 1.0 - p;
    np = n * p;
    lfactn = RandMath::lfact(n);
    logProb = std::log(p);
    log1mProb = std::log1pl(-p);
    SetGeneratorConstants();
}

template< typename IntType >
double BinomialDistribution<IntType>::SufficientStatistic(IntType x) const
{
    return x;
}

template< typename IntType >
double BinomialDistribution<IntType>::SourceParameters() const
{
    return this->p;
}

template< typename IntType >
double BinomialDistribution<IntType>::SourceToNatural(double sourceParameters) const
{
    return std::log(sourceParameters) - std::log1p(-sourceParameters);
}

template< typename IntType >
double BinomialDistribution<IntType>::NaturalParameters() const
{
    return this->logProb - this->log1mProb;
}

template< typename IntType >
double BinomialDistribution<IntType>::LogNormalizer(double theta) const
{
    double F = n * RandMath::log1pexp(theta);
    F -= this->lfactn;
    return F;
}

template< typename IntType >
double BinomialDistribution<IntType>::LogNormalizerGradient(double theta) const
{
    double expTheta = std::exp(theta);
    return this->n * expTheta / (1.0 + expTheta);
}

template< typename IntType >
double BinomialDistribution<IntType>::CarrierMeasure(IntType x) const
{
    double k = -RandMath::lfact(x);
    k -= RandMath::lfact(n - x);
    return k;
}

template< typename IntType >
double BinomialDistribution<IntType>::logProbFloor(int k) const
{
    double y = lfactn;
    y -= RandMath::lfact(n - k);
    y -= RandMath::lfact(k);
    y += k * logPFloor;
    y += (n - k) * logQFloor;
    return y;
}

template< typename IntType >
double BinomialDistribution<IntType>::P(const IntType & k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

template< typename IntType >
double BinomialDistribution<IntType>::logP(const IntType & k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    double y = lfactn;
    y -= RandMath::lfact(n - k);
    y -= RandMath::lfact(k);
    y += k * logProb;
    y += (n - k) * log1mProb;
    return y;
}

template< typename IntType >
double BinomialDistribution<IntType>::F(const IntType & k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = RandMath::lfact(n - kp1);
    logBetaFun += RandMath::lfact(k);
    logBetaFun -= lfactn;
    return RandMath::ibeta(q, nmk, kp1, logBetaFun, log1mProb, logProb);
}

template< typename IntType >
double BinomialDistribution<IntType>::S(const IntType & k) const
{
    if (k < 0)
        return 1.0;
    if (k >= n)
        return 0.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(kp1, nmk);
    return RandMath::ibeta(p, kp1, nmk, logBetaFun, logProb, log1mProb);
}

template< typename IntType >
IntType BinomialDistribution<IntType>::variateRejection() const
{
    /// a rejection algorithm by Devroye and Naderlsamanl (1980)
    /// p.533. Non-Uniform Random Variate Generation. Luc Devroye
    /// it can be used only when n * p is integer and p < 0.5
    bool reject = true;
    size_t iter = 0;
    float Y, V;
    IntType X;
    do {
        float U = a4 * UniformRand<float>::StandardVariate(this->localRandGenerator);
        if (U <= a1)
        {
            float N = NormalRand<float>::StandardVariate(this->localRandGenerator);
            Y = sigma1 * std::fabs(N);
            reject = (Y >= delta1);
            if (!reject)
            {
                float W = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
                X = std::floor(Y);
                V = -W - 0.5 * N * N + c;
            }
        }
        else if (U <= a2)
        {
            float N = NormalRand<float>::StandardVariate(this->localRandGenerator);
            Y = sigma2 * std::fabs(N);
            reject = (Y >= delta2);
            if (!reject)
            {
                float W = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
                X = std::floor(-Y);
                V = -W - 0.5 * N * N;
            }
        }
        else if (U <= a3)
        {
            float W1 = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
            float W2 = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
            Y = delta1 + W1 / coefa3;
            X = std::floor(Y);
            V = -W2 - coefa3 * Y + delta1 / nqFloor;
            reject = false;
        }
        else
        {
            float W1 = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
            float W2 = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
            Y = delta2 + W1 / coefa4;
            X = std::floor(-Y);
            V = -W2 - coefa4 * Y;
            reject = false;
        }

        if (!reject) {
            X += npFloor;
            if (X >= 0 && X <= n && V <= logProbFloor(X) - logPnpInv)
                return X;
        }
    } while (++iter <= ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Binomial distribution: sampling failed");
}

template< typename IntType >
IntType BinomialDistribution<IntType>::variateWaiting(IntType number) const
{
    /// waiting algorithm, using
    /// sum of geometrically distributed variables
    IntType X = -1, sum = 0;
    do {
        sum += G.Variate() + 1;
        ++X;
    } while (sum <= number);
    return X;
}

template< typename IntType >
IntType BinomialDistribution<IntType>::variateWaiting(IntType number, double probability, RandGenerator &randGenerator)
{
    IntType X = -1;
    double sum = 0;
    do {
        IntType add = GeometricRand<IntType>::Variate(probability, randGenerator) + 1;
        if (add < 0) /// we catched overflow
            return X + 1;
        sum += add;
        ++X;
    } while (sum <= number);
    return X;
}

template< typename IntType >
IntType BinomialDistribution<IntType>::variateBernoulliSum(IntType number, double probability, RandGenerator &randGenerator)
{
    IntType var = 0;
    if (RandMath::areClose(probability, 0.5)) {
        for (int i = 0; i != number; ++i)
            var += BernoulliRand::StandardVariate(randGenerator);
    }
    else {
        for (int i = 0; i != number; ++i)
            var += BernoulliRand::Variate(probability, randGenerator);
    }
    return var;
}

template< typename IntType >
IntType BinomialDistribution<IntType>::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING:
    {
        IntType var = variateWaiting(n);
        return (p <= 0.5) ? var : n - var;
    }
    case REJECTION:
    {
        /// if X ~ Bin(n, p') and Y ~ Bin(n - X, (p - p') / (1 - p'))
        /// then Z = X + Y ~ Bin(n, p)
        IntType Z = variateRejection();
        if (pRes > 0)
            Z += variateWaiting(n - Z);
        return (p > 0.5) ? n - Z : Z;
    }
    case BERNOULLI_SUM:
        return variateBernoulliSum(n, p, this->localRandGenerator);
    default:
        throw std::invalid_argument("Binomial distribution: invalid generator id");
    }
}

template< typename IntType >
IntType BinomialDistribution<IntType>::Variate(IntType number, double probability, RandGenerator &randGenerator)
{
    /// sanity check
    if (number < 0)
        throw std::invalid_argument("Binomial distribution: number should be positive, but it's equal to "
                                    + std::to_string(number));
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Binomial distribution: probability parameter should in interval [0, 1], but it's equal to "
                                    + std::to_string(probability));
    if (probability == 0.0)
        return 0;
    if (probability == 1.0)
        return number;

    if (number < 10)
        return variateBernoulliSum(number, probability, randGenerator);
    if (probability < 0.5)
        return variateWaiting(number, probability, randGenerator);
    return number - variateWaiting(number, 1.0 - probability, randGenerator);
}

template< typename IntType >
void BinomialDistribution<IntType>::Sample(std::vector<IntType> &outputData) const
{
    if (p == 0.0) {
        std::fill(outputData.begin(), outputData.end(), 0);
        return;
    }
    if (RandMath::areClose(p, 1.0)) {
        std::fill(outputData.begin(), outputData.end(), n);
        return;
    }

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING:
    {
        if (p <= 0.5) {
            for (IntType &var : outputData)
               var = variateWaiting(n);
        }
        else {
            for (IntType &var : outputData)
               var = n - variateWaiting(n);
        }
        return;
    }
    case REJECTION:
    {
        for (IntType &var : outputData)
            var = variateRejection();
        if (pRes > 0) {
            for (IntType &var : outputData)
                var += variateWaiting(n - var);
        }
        if (p > 0.5) {
            for (IntType &var : outputData)
               var = n - var;
        }
        return;
    }
    case BERNOULLI_SUM:
    default:
    {
        for (IntType &var : outputData)
           var = variateBernoulliSum(n, p, this->localRandGenerator);
        return;
    }
    }
}

template< typename IntType >
void BinomialDistribution<IntType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    G.Reseed(seed);
}

template< typename IntType >
long double BinomialDistribution<IntType>::Mean() const
{
    return np;
}

template< typename IntType >
long double BinomialDistribution<IntType>::Variance() const
{
    return np * q;
}

template< typename IntType >
IntType BinomialDistribution<IntType>::Median() const
{
    return std::floor(np);
}

template< typename IntType >
IntType BinomialDistribution<IntType>::Mode() const
{
    return std::floor(np + p);
}

template< typename IntType >
long double BinomialDistribution<IntType>::Skewness() const
{
    return (q - p) / std::sqrt(np * q);
}

template< typename IntType >
long double BinomialDistribution<IntType>::ExcessKurtosis() const
{
    long double y = 1.0 / (p * q);
    y -= 6.0;
    return y / n;
}

template< typename IntType >
std::complex<double> BinomialDistribution<IntType>::CFImpl(double t) const
{
    std::complex<double> y(q + p * std::cos(t), p * std::sin(t));
    return std::pow(y, n);
}

template< typename IntType >
void BinomialDistribution<IntType>::FitProbability(const std::vector<IntType> &sample)
{
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    if (!this->allElementsAreNotGreaterThan(n, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(n)));
    SetParameters(n, this->GetSampleMean(sample) / n);
}

template< typename IntType >
BetaRand<> BinomialDistribution<IntType>::FitProbabilityBayes(const std::vector<IntType> &sample, const BetaDistribution<> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    if (!this->allElementsAreNotGreaterThan(n, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(n)));
    int N = sample.size();
    double sum = this->GetSampleSum(sample);
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand posteriorDistribution(sum + alpha, N * n - sum + beta);
    SetParameters(n, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template< typename IntType >
BetaRand<> BinomialDistribution<IntType>::FitProbabilityMinimax(const std::vector<IntType> &sample)
{
    double shape = 0.5 * std::sqrt(n);
    BetaRand B(shape, shape);
    return FitProbabilityBayes(sample, B);
}

template class BinomialDistribution<int>;
template class BinomialDistribution<long int>;
template class BinomialDistribution<long long int>;

template< typename IntType >
String BinomialRand<IntType>::Name() const
{
    return "Binomial(" + this->toStringWithPrecision(this->GetNumber()) + ", " + this->toStringWithPrecision(this->GetProbability()) + ")";
}

template class BinomialRand<int>;
template class BinomialRand<long int>;
template class BinomialRand<long long int>;
