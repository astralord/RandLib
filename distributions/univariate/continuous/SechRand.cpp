#include "SechRand.h"
#include "CauchyRand.h"

template < typename RealType >
SechRand<RealType>::SechRand()
{
}

template < typename RealType >
String SechRand<RealType>::Name() const
{
    return "Hyperbolic secant";
}

template < typename RealType >
double SechRand<RealType>::f(const RealType & x) const
{
    return 0.5 / std::cosh(randlib::consts::half_pi * x);
}

template < typename RealType >
double SechRand<RealType>::logf(const RealType & x) const
{
    return randlib::consts::half_pi * x - RandMath::softplus(randlib::consts::pi * x);
}

template < typename RealType >
double SechRand<RealType>::F(const RealType & x) const
{
    double y = std::exp(randlib::consts::half_pi * x);
    return randlib::consts::two_inv_pi * RandMath::atan(y);
}

template < typename RealType >
RealType SechRand<RealType>::Variate() const
{
    RealType y = std::fabs(CauchyRand<RealType>::StandardVariate(this->localRandGenerator));
    return randlib::consts::two_inv_pi * std::log(y);
}

template < typename RealType >
long double SechRand<RealType>::Mean() const
{
    return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::Variance() const
{
    return 1.0;
}

template < typename RealType >
std::complex<double> SechRand<RealType>::CFImpl(double t) const
{
    return 1.0 / std::cosh(t);
}

template < typename RealType >
RealType SechRand<RealType>::quantileImpl(double p) const
{
    RealType x = randlib::consts::half_pi * p;
    x = std::tan(x);
    x = std::log(x);
    return randlib::consts::two_inv_pi * x;
}

template < typename RealType >
RealType SechRand<RealType>::quantileImpl1m(double p) const
{
    RealType x = randlib::consts::half_pi * p;
    x = std::tan(x);
    x = -std::log(x);
    return randlib::consts::two_inv_pi * x;
}

template < typename RealType >
RealType SechRand<RealType>::Median() const
{
    return 0.0;
}

template < typename RealType >
RealType SechRand<RealType>::Mode() const
{
    return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::Skewness() const
{
    return 0.0;
}

template < typename RealType >
long double SechRand<RealType>::ExcessKurtosis() const
{
    return 2.0;
}

template < typename RealType >
long double SechRand<RealType>::Entropy() const
{
    return 2.0 * randlib::consts::two_inv_pi * randlib::consts::catalan;
}


template class SechRand<float>;
template class SechRand<double>;
template class SechRand<long double>;
