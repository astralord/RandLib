#include "CauchyRand.h"
#include "UniformRand.h"

template < typename RealType >
CauchyRand<RealType>::CauchyRand(double location, double scale)
    : StableDistribution<RealType>(1, 0, scale, location)
{
}

template < typename RealType >
String CauchyRand<RealType>::Name() const
{
    return "Cauchy(" + this->toStringWithPrecision(this->GetLocation()) + ", " + this->toStringWithPrecision(this->GetScale()) + ")";
}

template < typename RealType >
double CauchyRand<RealType>::f(const RealType & x) const
{
    return this->pdfCauchy(x);
}

template < typename RealType >
double CauchyRand<RealType>::F(const RealType & x) const
{
    return this->cdfCauchy(x);
}

template < typename RealType >
double CauchyRand<RealType>::S(const RealType & x) const
{
    return this->cdfCauchyCompl(x);
}

template < typename RealType >
RealType CauchyRand<RealType>::Variate() const
{
    return this->mu + this->gamma * StandardVariate(this->localRandGenerator);
}

template < typename RealType >
RealType CauchyRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    double x, y;
    do {
        x = 2 * UniformRand<RealType>::StandardVariate(randGenerator) - 1;
        y = 2 * UniformRand<RealType>::StandardVariate(randGenerator) - 1;
    } while (y == 0.0 || x * x + y * y > 1.0);
    return x / y;
}

template < typename RealType >
std::complex<double> CauchyRand<RealType>::CFImpl(double t) const
{
    return this->cfCauchy(t);
}

template < typename RealType >
RealType CauchyRand<RealType>::quantileImpl(double p) const
{
    return this->quantileCauchy(p);
}

template < typename RealType >
RealType CauchyRand<RealType>::quantileImpl1m(double p) const
{
    return this->quantileCauchy1m(p);
}

template < typename RealType >
long double CauchyRand<RealType>::Entropy() const
{
    return 2 * M_LN2 + this->logGamma + M_LNPI;
}

template class CauchyRand<float>;
template class CauchyRand<double>;
template class CauchyRand<long double>;
