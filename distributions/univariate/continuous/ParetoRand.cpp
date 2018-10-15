#include "ParetoRand.h"
#include "UniformRand.h"

template < typename RealType >
ParetoRand<RealType>::ParetoRand(double shape, double scale)
{
    SetShape(shape);
    SetScale(scale);
}

template < typename RealType >
String ParetoRand<RealType>::Name() const
{
    return "Pareto(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void ParetoRand<RealType>::SetShape(double shape)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Pareto distribution: shape should be positive");
    alpha = shape;
    logAlpha = std::log(alpha);
}

template < typename RealType >
void ParetoRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Pareto distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

template < typename RealType >
double ParetoRand<RealType>::f(const RealType &x) const
{
    return (x < sigma) ? 0.0 : std::exp(logf(x));
}

template < typename RealType >
double ParetoRand<RealType>::logf(const RealType & x) const
{
    if (x < sigma)
        return -INFINITY;
    double logX = std::log(x);
    double y = logSigma - logX;
    y *= alpha;
    y -= logX;
    y += logAlpha;
    return y;
}

template < typename RealType >
double ParetoRand<RealType>::F(const RealType & x) const
{
    return (x > sigma) ? -std::expm1l(alpha * std::log(sigma / x)) : 0.0;
}

template < typename RealType >
double ParetoRand<RealType>::S(const RealType & x) const
{
    return (x > sigma) ? std::pow(sigma / x, alpha) : 1.0;
}

template < typename RealType >
RealType ParetoRand<RealType>::variateForAlphaEqualOne(RandGenerator &randGenerator)
{
    return 1.0 / UniformRand<RealType>::StandardVariate(randGenerator);
}

template < typename RealType >
RealType ParetoRand<RealType>::variateForAlphaEqualTwo(RandGenerator &randGenerator)
{
    return 1.0 / std::sqrt(UniformRand<RealType>::StandardVariate(randGenerator));
}

template < typename RealType >
RealType ParetoRand<RealType>::variateForGeneralAlpha(double shape, RandGenerator &randGenerator)
{
    return std::exp(ExponentialRand<RealType>::StandardVariate(randGenerator) / shape);
}

template < typename RealType >
RealType ParetoRand<RealType>::Variate() const
{
    return sigma * StandardVariate(alpha, this->localRandGenerator);
}

template < typename RealType >
RealType ParetoRand<RealType>::StandardVariate(double shape, RandGenerator &randGenerator)
{
    if (RandMath::areClose(shape, 1.0))
        return variateForAlphaEqualOne(randGenerator);
    if (RandMath::areClose(shape, 2.0))
        return variateForAlphaEqualTwo(randGenerator);
    return variateForGeneralAlpha(shape, randGenerator);
}

template < typename RealType >
void ParetoRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    if (RandMath::areClose(alpha, 1.0)) {
        for (RealType &var : outputData)
            var = sigma * variateForAlphaEqualOne(this->localRandGenerator);
    }
    else if (RandMath::areClose(alpha, 2.0)) {
        for (RealType &var : outputData)
            var = sigma * variateForAlphaEqualTwo(this->localRandGenerator);
    }
    else {
        for (RealType &var : outputData)
            var = sigma * variateForGeneralAlpha(alpha, this->localRandGenerator);
    }
}

template < typename RealType >
long double ParetoRand<RealType>::Mean() const
{
    return (alpha > 1) ? alpha * sigma / (alpha - 1) : INFINITY;
}

template < typename RealType >
long double ParetoRand<RealType>::Variance() const
{
    if (alpha > 2)
    {
        long double var = sigma / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

template < typename RealType >
RealType ParetoRand<RealType>::Median() const
{
    return std::exp(logSigma + M_LN2 / alpha);
}

template < typename RealType >
RealType ParetoRand<RealType>::Mode() const
{
    return sigma;
}

template < typename RealType >
long double ParetoRand<RealType>::Skewness() const
{
    if (alpha <= 3)
        return INFINITY;
    double skewness = (alpha - 2.0) / alpha;
    skewness = std::sqrt(skewness);
    skewness *= (1 + alpha) / (alpha - 3);
    return skewness + skewness;
}

template < typename RealType >
long double ParetoRand<RealType>::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    double numerator = alpha + 1;
    numerator *= alpha;
    numerator -= 6;
    numerator *= alpha;
    numerator -= 2;
    double denominator = alpha * (alpha - 3) * (alpha - 4);
    return 6.0 * numerator / denominator;
}

template < typename RealType >
RealType ParetoRand<RealType>::quantileImpl(double p) const
{
    double y = logSigma - std::log1pl(-p) / alpha;
    return std::exp(y);
}

template < typename RealType >
RealType ParetoRand<RealType>::quantileImpl1m(double p) const
{
    double y = logSigma - std::log(p) / alpha;
    return std::exp(y);
}

template < typename RealType >
long double ParetoRand<RealType>::Entropy() const
{
    return logSigma - logAlpha + 1.0 / alpha + 1;
}

template < typename RealType >
void ParetoRand<RealType>::FitShape(const std::vector<RealType> &sample, bool unbiased)
{
    if (!this->allElementsAreNotSmallerThan(sigma, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(sigma)));
    double invShape = this->GetSampleLogMean(sample) - logSigma;
    if (invShape == 0.0)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, "Possibly all the elements of the sample coincide with the lower boundary σ."));
    double shape = 1.0 / invShape;
    if (unbiased)
        shape *= (1.0 - 1.0 / sample.size());
    SetShape(shape);
}

template < typename RealType >
void ParetoRand<RealType>::FitScale(const std::vector<RealType> &sample, bool unbiased)
{
    RealType minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    double scale = unbiased ? (1.0 - 1.0 / (sample.size() * alpha)) * minVar : minVar;
    SetScale(scale);
}

template < typename RealType >
void ParetoRand<RealType>::Fit(const std::vector<RealType> &sample, bool unbiased)
{
    RealType minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));

    double logBiasedSigma = std::log(minVar);
    double invShape = this->GetSampleLogMean(sample) - logBiasedSigma;
    if (invShape == 0.0)
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, "Possibly all the elements of the sample are the same."));
    double shape = 1.0 / invShape;

    double scale = minVar;
    if (unbiased) {
        int n = sample.size();
        shape *= 1.0 - 2.0 / n;
        scale *= 1.0 - invShape / (n - 1);
    }
    SetScale(scale);
    SetShape(shape);
}

template < typename RealType >
GammaRand<RealType> ParetoRand<RealType>::FitShapeBayes(const std::vector<RealType> &sample, const GammaDistribution<RealType> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNotSmallerThan(sigma, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(sigma)));
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newRate = priorDistribution.GetRate() + n * (this->GetSampleLogMean(sample) - logSigma);
    GammaRand<RealType> posteriorDistribution(newShape, newRate);
    SetShape(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template class ParetoRand<float>;
template class ParetoRand<double>;
template class ParetoRand<long double>;
