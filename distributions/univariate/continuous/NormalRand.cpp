#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "../BasicRandGenerator.h"
#include "GammaRand.h"
#include "StudentTRand.h"

template < typename RealType >
NormalRand<RealType>::NormalRand(double location, double variance)
    : StableDistribution<RealType>(2.0, 0.0, 1.0, location)
{
    SetVariance(variance);
}

template < typename RealType >
String NormalRand<RealType>::Name() const
{
    return "Normal(" + this->toStringWithPrecision(this->GetLocation()) + ", " + this->toStringWithPrecision(this->Variance()) + ")";
}

template < typename RealType >
void NormalRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Normal distribution: scale should be positive, but it's equal to " + std::to_string(scale));
    sigma = scale;
    StableDistribution<RealType>::SetScale(sigma * M_SQRT1_2);
}

template < typename RealType >
void NormalRand<RealType>::SetVariance(double variance)
{
    if (variance <= 0.0)
        throw std::invalid_argument("Variance of Normal distribution should be positive, but it's equal to " + std::to_string(variance));
    SetScale(std::sqrt(variance));
}

template < typename RealType >
double NormalRand<RealType>::f(const RealType & x) const
{
    return this->pdfNormal(x);
}

template < typename RealType >
double NormalRand<RealType>::logf(const RealType & x) const
{
    return this->logpdfNormal(x);
}

template < typename RealType >
double NormalRand<RealType>::F(const RealType & x) const
{
    return this->cdfNormal(x);
}

template < typename RealType >
double NormalRand<RealType>::S(const RealType & x) const
{
    return this->cdfNormalCompl(x);
}

template < typename RealType >
RealType NormalRand<RealType>::Variate() const
{
    return this->mu + sigma * StandardVariate(this->localRandGenerator);
}

template < typename RealType >
RealType NormalRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    /// Ziggurat algorithm by George Marsaglia using 256 strips
    size_t iter = 0;
    do {
        unsigned long long B = randGenerator.Variate();
        int stairId = B & 255;
        RealType x = UniformRand<RealType>::StandardVariate(randGenerator) * ziggurat[stairId].second; /// Get horizontal coordinate
        if (x < ziggurat[stairId + 1].second)
            return ((signed)B > 0) ? x : -x;
        if (stairId == 0) /// handle the base layer
        {
            static thread_local RealType z = -1;
            if (z > 0) /// we don't have to generate another exponential variable as we already have one
            {
                x = ExponentialRand<RealType>::StandardVariate(randGenerator) / ziggurat[1].second;
                z -= 0.5 * x * x;
            }
            if (z <= 0) /// if previous generation wasn't successful
            {
                do {
                    x = ExponentialRand<RealType>::StandardVariate(randGenerator) / ziggurat[1].second;
                    z = ExponentialRand<RealType>::StandardVariate(randGenerator) - 0.5 * x * x; /// we storage this value as after acceptance it becomes exponentially distributed
                } while (z <= 0);
            }
            x += ziggurat[1].second;
            return ((signed)B > 0) ? x : -x;
        }
        /// handle the wedges of other stairs
        RealType height = ziggurat[stairId].first - ziggurat[stairId - 1].first;
        if (ziggurat[stairId - 1].first + height * UniformRand<RealType>::StandardVariate(randGenerator) < std::exp(-.5 * x * x))
            return ((signed)B > 0) ? x : -x;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Normal distribution: sampling failed");
}

template < typename RealType >
void NormalRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    for (RealType & var : outputData)
        var = this->Variate();
}

template < typename RealType >
std::complex<double> NormalRand<RealType>::CFImpl(double t) const
{
    return this->cfNormal(t);
}

template < typename RealType >
RealType NormalRand<RealType>::quantileImpl(double p) const
{
    return this->quantileNormal(p);
}

template < typename RealType >
RealType NormalRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileNormal1m(p);
}

template < typename RealType >
long double NormalRand<RealType>::Moment(size_t n) const
{
    if (n == 0)
        return 1;
    return (n & 1) ? std::exp(n * this->GetLogScale() + RandMath::ldfact(n - 1)) : 0.0;
}

template < typename RealType >
double NormalRand<RealType>::KullbackLeiblerDivergence(const NormalRand &Y)
{
    double divLogVar = Y.GetLogScale() - this->GetLogScale();
    double divMean = Y.Mean() - this->Mean();
    double div = (divMean * divMean + this->GetScale()) / Y.GetScale();
    --div;
    div *= 0.5;
    div += divLogVar;
    return div;
}

template < typename RealType >
void NormalRand<RealType>::FitLocation(const std::vector<RealType> &sample)
{
    this->SetLocation(this->GetSampleMean(sample));
}

template < typename RealType >
void NormalRand<RealType>::FitLocation(const std::vector<RealType> &sample, DoublePair &confidenceInterval, double significanceLevel)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_LEVEL, "Input level is equal to " + this->toStringWithPrecision(significanceLevel)));

    FitLocation(sample);

    int n = sample.size();
    NormalRand<RealType> NormalRV(0, 1);
    double halfAlpha = 0.5 * significanceLevel;
    double interval = NormalRV.Quantile1m(halfAlpha) * sigma / std::sqrt(n);
    confidenceInterval.first = this->mu - interval;
    confidenceInterval.second = this->mu + interval;
}

template < typename RealType >
void NormalRand<RealType>::FitVariance(const std::vector<RealType> &sample)
{
    this->SetVariance(this->GetSampleVariance(sample, this->mu));
}

template < typename RealType >
void NormalRand<RealType>::FitVariance(const std::vector<RealType> &sample, DoublePair &confidenceInterval, double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_LEVEL, "Input level is equal to " + this->toStringWithPrecision(significanceLevel)));

    FitVariance(sample);

    size_t n = sample.size();
    double halfAlpha = 0.5 * significanceLevel;
    ChiSquaredRand ChiSqRV(n);
    double numerator = sigma * sigma * (unbiased ? n : (n - 1));
    confidenceInterval.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceInterval.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

template < typename RealType >
void NormalRand<RealType>::FitScale(const std::vector<RealType> &sample, bool unbiased)
{
    if (unbiased == false)
        return FitVariance(sample);
    size_t n = sample.size();
    double halfN = 0.5 * n;
    double s = this->GetSampleVariance(sample, this->mu);
    s *= halfN;
    s = 0.5 * std::log(s);
    s -= std::lgammal(halfN + 0.5);
    s += std::lgammal(halfN);
    SetScale(std::exp(s));
}

template < typename RealType >
void NormalRand<RealType>::Fit(const std::vector<RealType> &sample, bool unbiased)
{
    double adjustment = 1.0;
    if (unbiased == true) {
        size_t n = sample.size();
        if (n <= 1)
            throw std::invalid_argument(this->fitErrorDescription(this->TOO_FEW_ELEMENTS, "There should be at least 2 elements"));
        adjustment = static_cast<double>(n) / (n - 1);
    }
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    this->SetLocation(stats.first);
    this->SetVariance(stats.second * adjustment);
}

template < typename RealType >
void NormalRand<RealType>::Fit(const std::vector<RealType> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_LEVEL, "Input level is equal to " + this->toStringWithPrecision(significanceLevel)));

    Fit(sample, unbiased);

    size_t n = sample.size();
    double sigmaAdj = unbiased ? sigma : (sigma * n) / (n - 1);
    /// calculate confidence interval for mean
    double halfAlpha = 0.5 * significanceLevel;
    StudentTRand<RealType> tRV(n - 1);
    double interval = tRV.Quantile1m(halfAlpha) * sigmaAdj / std::sqrt(n);
    confidenceIntervalForMean.first = this->mu - interval;
    confidenceIntervalForMean.second = this->mu + interval;

    /// calculate confidence interval for variance
    ChiSquaredRand<RealType> ChiSqRV(n - 1);
    double numerator = (n - 1) * sigmaAdj * sigmaAdj;
    confidenceIntervalForVariance.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceIntervalForVariance.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

template < typename RealType >
NormalRand<RealType> NormalRand<RealType>::FitLocationBayes(const std::vector<RealType> &sample, const NormalRand<RealType> &priorDistribution, bool MAP)
{
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = GetPrecision();
    double numerator = this->GetSampleSum(sample) * tau + tau0 * mu0;
    double denominator = sample.size() * tau + tau0;
    NormalRand<RealType> posteriorDistribution(numerator / denominator, 1.0 / denominator);
    this->SetLocation(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template < typename RealType >
InverseGammaRand<RealType> NormalRand<RealType>::FitVarianceBayes(const std::vector<RealType> &sample, const InverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    double halfN = 0.5 * sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double alphaPosterior = alphaPrior + halfN;
    double betaPosterior = betaPrior + halfN * this->GetSampleVariance(sample, this->mu);
    InverseGammaRand<RealType> posteriorDistribution(alphaPosterior, betaPosterior);
    SetVariance(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template < typename RealType >
NormalInverseGammaRand<RealType> NormalRand<RealType>::FitBayes(const std::vector<RealType> &sample, const NormalInverseGammaRand<RealType> &priorDistribution, bool MAP)
{
    size_t n = sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double muPrior = priorDistribution.GetLocation();
    double lambdaPrior = priorDistribution.GetPrecision();
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    double lambdaPosterior = lambdaPrior + n;
    double muPosterior = (lambdaPrior * muPrior + n * stats.first) / lambdaPosterior;
    double halfN = 0.5 * n;
    double alphaPosterior = alphaPrior + halfN;
    double aux = muPrior - stats.first;
    double betaPosterior = betaPrior + halfN * (stats.second + lambdaPrior / lambdaPosterior * aux * aux);
    NormalInverseGammaRand<RealType> posteriorDistribution(muPosterior, lambdaPosterior, alphaPosterior, betaPosterior);
    DoublePair newParams = MAP ? static_cast<DoublePair>(posteriorDistribution.Mode()) : static_cast<DoublePair>(posteriorDistribution.Mean());
    this->SetLocation(newParams.first);
    this->SetVariance(newParams.second);
    return posteriorDistribution;
}

template class NormalRand<float>;
template class NormalRand<double>;
template class NormalRand<long double>;

