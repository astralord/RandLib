#include "ExponentialRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

template < typename RealType >
String ExponentialRand<RealType>::Name() const
{
    return "Exponential(" + this->toStringWithPrecision(this->GetRate()) + ")";
}

template < typename RealType >
double ExponentialRand<RealType>::SufficientStatistic(RealType x) const
{
    return x;
}

template < typename RealType >
double ExponentialRand<RealType>::SourceParameters() const
{
    return this->beta;
}

template < typename RealType >
double ExponentialRand<RealType>::SourceToNatural(double rate) const
{
    return -rate;
}

template < typename RealType >
double ExponentialRand<RealType>::LogNormalizer(double thetaP) const
{
    return -std::log(-thetaP);
}

template < typename RealType >
double ExponentialRand<RealType>::LogNormalizerGradient(double thetaP) const
{
    return -1.0 / thetaP;
}

template < typename RealType >
double ExponentialRand<RealType>::CarrierMeasure(RealType) const
{
    return 0.0;
}

template < typename RealType >
double ExponentialRand<RealType>::CrossEntropyAdjusted(double rate) const
{
    return rate * this->theta - std::log(rate);
}

template < typename RealType >
double ExponentialRand<RealType>::EntropyAdjusted() const
{
    return 1.0 - this->logBeta;
}

template < typename RealType >
double ExponentialRand<RealType>::f(const RealType &x) const
{
    return (x < 0.0) ? 0.0 : this->beta * std::exp(-this->beta * x);
}

template < typename RealType >
double ExponentialRand<RealType>::logf(const RealType & x) const
{
    return (x < 0.0) ? -INFINITY : this->logBeta - this->beta * x;
}

template < typename RealType >
double ExponentialRand<RealType>::F(const RealType & x) const
{
    return (x > 0.0) ? -std::expm1l(-this->beta * x) : 0.0;
}

template < typename RealType >
double ExponentialRand<RealType>::S(const RealType & x) const
{
    return (x > 0.0) ? std::exp(-this->beta * x) : 1.0;
}

template < typename RealType >
RealType ExponentialRand<RealType>::Variate() const
{
    return this->theta * StandardVariate(this->localRandGenerator);
}

template < typename RealType >
void ExponentialRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    for (RealType & var : outputData)
        var = this->Variate();
}

template < typename RealType >
RealType ExponentialRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    /// Ziggurat algorithm
    size_t iter = 0;
    do {
        int stairId = randGenerator.Variate() & 255;
        /// Get horizontal coordinate
        RealType x = UniformRand<RealType>::StandardVariate(randGenerator) * ziggurat[stairId].second;
        if (x < ziggurat[stairId + 1].second) /// if we are under the upper stair - accept
            return x;
        if (stairId == 0) /// if we catch the tail
            return ziggurat[1].second + StandardVariate(randGenerator);
        RealType height = ziggurat[stairId].first - ziggurat[stairId - 1].first;
        if (ziggurat[stairId - 1].first + height * UniformRand<RealType>::StandardVariate(randGenerator) < std::exp(-x)) /// if we are under the curve - accept
            return x;
        /// rejection - go back
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    /// fail due to some error
    throw std::runtime_error("Exponential distribution: sampling failed");
}

template < typename RealType >
std::complex<double> ExponentialRand<RealType>::CFImpl(double t) const
{
    return 1.0 / std::complex<double>(1.0, -this->theta * t);
}

template < typename RealType >
long double ExponentialRand<RealType>::Entropy() const
{
    return this->EntropyAdjusted();
}

template < typename RealType >
long double ExponentialRand<RealType>::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return std::exp(RandMath::lfact(n) - n * this->logBeta);
}


template class ExponentialRand<float>;
template class ExponentialRand<double>;
template class ExponentialRand<long double>;
