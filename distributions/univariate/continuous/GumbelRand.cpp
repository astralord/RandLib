#include "GumbelRand.h"
#include "ExponentialRand.h"

template < typename RealType >
GumbelRand<RealType>::GumbelRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
String GumbelRand<RealType>::Name() const
{
    return "Gumbel(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void GumbelRand<RealType>::SetLocation(double location)
{
    mu = location;
}

template < typename RealType >
void GumbelRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Gumbel distribution: scale should be positive");
    beta = scale;
    logBeta = std::log(beta);
}

template < typename RealType >
double GumbelRand<RealType>::f(const RealType & x) const
{
    return std::exp(logf(x));
}

template < typename RealType >
double GumbelRand<RealType>::logf(const RealType & x) const
{
    double z = (mu - x) / beta;
    double y = std::exp(z);
    return z - y - logBeta;
}

template < typename RealType >
double GumbelRand<RealType>::F(const RealType & x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
    return std::exp(-y);
}

template < typename RealType >
double GumbelRand<RealType>::S(const RealType &x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
    return -std::expm1l(-y);
}

template < typename RealType >
RealType GumbelRand<RealType>::Variate() const
{
    return mu + beta * GumbelRand<RealType>::StandardVariate(this->localRandGenerator);
}

template < typename RealType >
RealType GumbelRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    RealType w = ExponentialRand<RealType>::StandardVariate(randGenerator);
    return -std::log(w);
}

template < typename RealType >
long double GumbelRand<RealType>::Mean() const
{
    return mu + beta * M_EULER;
}

template < typename RealType >
long double GumbelRand<RealType>::Variance() const
{
    double v = M_PI * beta;
    return v * v / 6;
}

template < typename RealType >
RealType GumbelRand<RealType>::quantileImpl(double p) const
{
    return mu - beta * std::log(-std::log(p));
}

template < typename RealType >
RealType GumbelRand<RealType>::quantileImpl1m(double p) const
{
    return mu - beta * std::log(-std::log1pl(-p));
}

template < typename RealType >
RealType GumbelRand<RealType>::Median() const
{
    static constexpr double M_LN_LN2 = std::log(M_LN2);
    return mu - beta * M_LN_LN2;
}

template < typename RealType >
RealType GumbelRand<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double GumbelRand<RealType>::Skewness() const
{
    static constexpr long double skew = 12 * M_SQRT2 * M_SQRT3 * M_APERY / (M_PI_SQ * M_PI);
    return skew;
}

template < typename RealType >
long double GumbelRand<RealType>::ExcessKurtosis() const
{
    return 2.4l;
}

template < typename RealType >
long double GumbelRand<RealType>::Entropy() const
{
    return logBeta + M_EULER + 1.0;
}

template class GumbelRand<float>;
template class GumbelRand<double>;
template class GumbelRand<long double>;
