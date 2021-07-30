#include "PoissonRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

template < typename IntType >
PoissonRand<IntType>::PoissonRand(double rate)
{
    SetRate(rate);
}

template < typename IntType >
String PoissonRand<IntType>::Name() const
{
    return "Poisson(" + this->toStringWithPrecision(GetRate()) + ")";
}

template < typename IntType >
void PoissonRand<IntType>::SetGeneratorConstants()
{
    delta = std::round(std::sqrt(2 * mu * std::log(M_1_PI * 128 * mu)));
    delta = std::max(6.0, std::min(mu, delta));

    c1 = std::sqrt(0.5 * M_PI * mu * M_E);

    c2 = 1.0 / (2 * mu + delta);
    c2 += 0.5 * (M_LNPI - M_LN2 + std::log(mu + 0.5 * delta));
    c2 = c1 + std::exp(c2);

    c3 = c2 + 1.0;
    c4 = c3 + 1.0;

    zeta = (4.0 * mu) / delta + 2;
    c = std::exp(-(2.0 + delta) / zeta);
    c *= zeta;
    c += c4;

    sqrtMu = std::sqrt(mu);
    sqrtMupHalfDelta = std::sqrt(mu + 0.5 * delta);
    lfactMu = RandMath::lfact(mu);
}

template < typename IntType >
void PoissonRand<IntType>::SetRate(double rate)
{
    if (rate <= 0.0)
        throw std::invalid_argument("Poisson distribution: rate should be positive, but it's equal to " + std::to_string(rate));
    lambda = rate;

    logLambda = std::log(lambda);
    mu = std::floor(lambda);
    Fmu = F(mu);
    Pmu = this->P(mu);

    if (!generateByInversion())
        SetGeneratorConstants();
}

template < typename IntType >
double PoissonRand<IntType>::SufficientStatistic(IntType x) const
{
    return x;
}

template < typename IntType >
double PoissonRand<IntType>::SourceParameters() const
{
    return lambda;
}

template < typename IntType >
double PoissonRand<IntType>::SourceToNatural(double sourceParameters) const
{
    return std::log(sourceParameters);
}

template < typename IntType >
double PoissonRand<IntType>::NaturalParameters() const
{
    return logLambda;
}

template < typename IntType >
double PoissonRand<IntType>::LogNormalizer(double theta) const
{
    return std::exp(theta);
}

template < typename IntType >
double PoissonRand<IntType>::LogNormalizerGradient(double theta) const
{
    return std::exp(theta);
}

template < typename IntType >
double PoissonRand<IntType>::CarrierMeasure(IntType x) const
{
    return -RandMath::lfact(x);
}

template < typename IntType >
double PoissonRand<IntType>::CrossEntropyAdjusted(double parameters) const
{
    return parameters - lambda * std::log(parameters);
}

template < typename IntType >
double PoissonRand<IntType>::EntropyAdjusted() const
{
    return lambda - lambda * logLambda;
}

template < typename IntType >
double PoissonRand<IntType>::logP(const IntType & k) const
{
    if (k < 0)
        return -INFINITY;
    double y = k * logLambda - lambda;
    return y - RandMath::lfact(k);
}

template < typename IntType >
double PoissonRand<IntType>::F(const IntType & k) const
{
    return (k >= 0.0) ? RandMath::qgamma(k + 1, lambda, logLambda) : 0.0;
}

template < typename IntType >
double PoissonRand<IntType>::S(const IntType & k) const
{
    return (k >= 0.0) ? RandMath::pgamma(k + 1, lambda, logLambda) : 1.0;
}

template < typename IntType >
double PoissonRand<IntType>::acceptanceFunction(IntType X) const
{
    if (X == 0)
        return 0.0;
    double q = X * logLambda;
    q += lfactMu;
    q -= RandMath::lfact(X + mu);
    return q;
}

template < typename IntType >
bool PoissonRand<IntType>::generateByInversion() const
{
    /// the inversion generator is much faster than rejection,
    /// however precision loss for large rate increases drastically
    return lambda < 10;
}

template < typename IntType >
IntType PoissonRand<IntType>::variateRejection() const
{
    size_t iter = 0;
    IntType X = 0;
    do {
        bool reject = false;
        float W = 0.0;
        float U = c * UniformRand<float>::StandardVariate(this->localRandGenerator);
        if (U <= c1) {
            float N = NormalRand<float>::StandardVariate(this->localRandGenerator);
            float Y = -std::fabs(N) * sqrtMu;
            X = std::floor(Y);
            if (X < -mu) {
                reject = true;
            }
            else {
                W = -0.5 * (N * N - 1.0);
            }
        }
        else if (U <= c2) {
            float N = NormalRand<float>::StandardVariate(this->localRandGenerator);
            float Y = 1.0 + std::fabs(N) * sqrtMupHalfDelta;
            X = std::ceil(Y);
            if (X > delta) {
                reject = true;
            }
            else {
                W = Y * (2.0 - Y) / (2.0 * mu + delta);
            }
        }
        else if (U <= c3) {
            return mu;
        }
        else if (U <= c4) {
            X = 1;
        }
        else {
            float V = ExponentialRand<float>::StandardVariate(this->localRandGenerator);
            float Y = delta + V * zeta;
            X = std::ceil(Y);
            W = -(2.0 + Y) / zeta;
        }

        if (!reject && W - ExponentialRand<float>::StandardVariate(this->localRandGenerator) <= acceptanceFunction(X)) {
            return X + mu;
        }

    } while (++iter < ProbabilityDistribution<IntType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Poisson distribution: sampling failed");
}

template < typename IntType >
IntType PoissonRand<IntType>::variateInversion() const
{
    double U = UniformRand<double>::StandardVariate(this->localRandGenerator);
    IntType k = mu;
    double s = Fmu, p = Pmu;
    if (s < U)
    {
        do {
            ++k;
            p *= lambda / k;
            s += p;
        } while (s < U && p > 0);
    }
    else
    {
        s -= p;
        while (k > 0 && s > U) {
            p /= lambda / k;
            --k;
            s -= p;
        }
    }
    return k;
}

template < typename IntType >
IntType PoissonRand<IntType>::Variate() const
{
    return generateByInversion() ? variateInversion() : variateRejection();
}

template < typename IntType >
IntType PoissonRand<IntType>::Variate(double rate, RandGenerator &randGenerator)
{
    /// check validness of parameter
    if (rate <= 0.0)
        throw std::invalid_argument("Poisson distribution: rate should be positive");
    if (rate > 1000) {
        /// approximate with normal distribution
        float X = NormalRand<float>::StandardVariate(randGenerator);
        return std::floor(rate + std::sqrt(rate) * X);
    }
    int k = -1;
    double s = 0;
    do {
        s += ExponentialRand<double>::StandardVariate(randGenerator);
        ++k;
    } while (s < rate);
    return k;
}

template < typename IntType >
void PoissonRand<IntType>::Sample(std::vector<IntType> &outputData) const
{
    if (generateByInversion()) {
        for (IntType & var : outputData)
            var = variateInversion();
    }
    else {
        for (IntType & var : outputData)
            var = variateRejection();
    }
}

template < typename IntType >
long double PoissonRand<IntType>::Mean() const
{
    return lambda;
}

template < typename IntType >
long double PoissonRand<IntType>::Variance() const
{
    return lambda;
}

template < typename IntType >
std::complex<double> PoissonRand<IntType>::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t) - 1.0, std::sin(t));
    return std::exp(lambda * y);
}

template < typename IntType >
IntType PoissonRand<IntType>::Median() const
{
    /// this value is approximate
    return std::max(std::floor(lambda + 1.0 / 3 - 0.02 / lambda), 0.0);
}

template < typename IntType >
IntType PoissonRand<IntType>::Mode() const
{
    if (RandMath::areClose(mu, lambda)) {
        return (Pmu < this->P(mu + 1)) ? mu + 1 : mu;
    }
    return mu;
}

template < typename IntType >
long double PoissonRand<IntType>::Skewness() const
{
    return 1.0 / std::sqrt(lambda);
}

template < typename IntType >
long double PoissonRand<IntType>::ExcessKurtosis() const
{
    return 1.0 / lambda;
}

template < typename IntType >
void PoissonRand<IntType>::Fit(const std::vector<IntType> &sample)
{
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    SetRate(this->GetSampleMean(sample));
}

template < typename IntType >
void PoissonRand<IntType>::Fit(const std::vector<IntType> &sample, DoublePair &confidenceInterval, double significanceLevel)
{
    size_t n = sample.size();

    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_LEVEL, "Alpha is equal to " + this->toStringWithPrecision(significanceLevel)));

    Fit(sample);

    double halfAlpha = 0.5 * significanceLevel;
    ErlangRand<double> ErlangRV(n);
    confidenceInterval.first = ErlangRV.Quantile(halfAlpha);
    ErlangRV.SetShape(n + 1);
    confidenceInterval.second = ErlangRV.Quantile1m(halfAlpha);
}

template < typename IntType >
GammaRand<> PoissonRand<IntType>::FitBayes(const std::vector<IntType> &sample, const GammaDistribution<> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNonNegative(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->NON_NEGATIVITY_VIOLATION));
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    GammaRand<> posteriorDistribution(alpha + this->GetSampleSum(sample), beta + sample.size());
    SetRate(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}


template class PoissonRand<int>;
template class PoissonRand<long int>;
template class PoissonRand<long long int>;
