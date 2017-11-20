#include "LevyRand.h"
#include "NormalRand.h"

template < typename RealType >
LevyRand<RealType>::LevyRand(double location, double scale)
    : StableDistribution<RealType>(0.5, 1, scale, location)
{
}

template < typename RealType >
String LevyRand<RealType>::Name() const
{
    return "Levy(" + this->toStringWithPrecision(this->GetLocation()) + ", " + this->toStringWithPrecision(this->GetScale()) + ")";
}

template < typename RealType >
double LevyRand<RealType>::f(const RealType &x) const
{
    return pdfLevy(x);
}

template < typename RealType >
double LevyRand<RealType>::logf(const RealType & x) const
{
    return logpdfLevy(x);
}

template < typename RealType >
double LevyRand<RealType>::F(const RealType & x) const
{
    return cdfLevy(x);
}

template < typename RealType >
double LevyRand<RealType>::S(const RealType & x) const
{
    return cdfLevyCompl(x);
}

template < typename RealType >
RealType LevyRand<RealType>::Variate() const
{
    RealType rv = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    rv *= rv;
    rv = gamma / rv;
    return this->mu + rv;
}

template < typename RealType >
RealType LevyRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    RealType rv = NormalRand<RealType>::StandardVariate(randGenerator);
    return 1.0 / (rv * rv);
}

template < typename RealType >
RealType LevyRand<RealType>::quantileImpl(double p) const
{
    return this->quantileLevy(p);
}

template < typename RealType >
RealType LevyRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileLevy1m(p);
}

template < typename RealType >
std::complex<double> LevyRand<RealType>::CFImpl(double t) const
{
    return this->cfLevy(t);
}

template < typename RealType >
void LevyRand<RealType>::FitScale(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsAreNotSmallerThan(this->mu, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->mu)));
    long double invSum = 0.0;
    for (RealType var : sample)
        invSum += 1.0 / (var - this->mu);
    invSum = 1.0 / invSum;
    SetScale(sample.size() * invSum);
}
