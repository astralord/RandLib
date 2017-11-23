#include "LogisticRand.h"

template < typename RealType >
LogisticRand<RealType>::LogisticRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
String LogisticRand<RealType>::Name() const
{
    return "Logistic(" + this->toStringWithPrecision(GetLocation()) + ", " + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void LogisticRand<RealType>::SetLocation(double location)
{
    mu = location;
}

template < typename RealType >
void LogisticRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Logistic distribution: scale of should be positive");
    s = scale;
    logS = std::log(s);
}

template < typename RealType >
double LogisticRand<RealType>::f(const RealType & x) const
{
    double numerator = std::exp((mu - x) / s);
    double denominator = (1 + numerator);
    denominator *= denominator;
    denominator *= s;
    return numerator / denominator;
}

template < typename RealType >
double LogisticRand<RealType>::logf(const RealType & x) const
{
    double x0 = (mu - x) / s;
    double y = RandMath::log1pexp(x0);
    y *= 2;
    y += logS;
    return x0 - y;
}

template < typename RealType >
double LogisticRand<RealType>::F(const RealType & x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

template < typename RealType >
double LogisticRand<RealType>::S(const RealType & x) const
{
    double expX = std::exp((mu - x) / s);
    return expX / (1 + expX);
}

template < typename RealType >
RealType LogisticRand<RealType>::Variate() const
{
    /// there can be used rejection method from Laplace or Cauchy (Luc Devroye, p. 471) or ziggurat
    return mu + s * std::log(1.0 / UniformRand<RealType>::StandardVariate(this->localRandGenerator) - 1);
}

template < typename RealType >
long double LogisticRand<RealType>::Mean() const
{
    return mu;
}

template < typename RealType >
long double LogisticRand<RealType>::Variance() const
{
    double sPi = s * M_PI;
    return sPi * sPi / 3;
}

template < typename RealType >
std::complex<double> LogisticRand<RealType>::CFImpl(double t) const
{
    double pist = M_PI * s * t;
    std::complex<double> y(0.0, t * mu);
    y = std::exp(y);
    y *= pist;
    y /= std::sinh(pist);
    return y;
}

template < typename RealType >
double LogisticRand<RealType>::Entropy() const
{
    return 2 + logS;
}

template < typename RealType >
RealType LogisticRand<RealType>::quantileImpl(double p) const
{
    return mu - s * (std::log1pl(-p) - std::log(p));
}

template < typename RealType >
RealType LogisticRand<RealType>::quantileImpl1m(double p) const
{
    return mu - s * (std::log(p) - std::log1pl(-p));
}

template < typename RealType >
RealType LogisticRand<RealType>::Median() const
{
    return mu;
}

template < typename RealType >
RealType LogisticRand<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double LogisticRand<RealType>::Skewness() const
{
    return 0;
}

template < typename RealType >
long double LogisticRand<RealType>::ExcessKurtosis() const
{
    return 1.2;
}

template < typename RealType >
void LogisticRand<RealType>::FitLocation(const std::vector<RealType> &sample)
{
    double nHalf = 0.5 * sample.size();
    RealType root = 0;
    if (!RandMath::findRoot<RealType>([this, sample, nHalf](RealType m)
    {
        double f1 = 0, f2 = 0;
        for (const double & x : sample)
        {
            double aux = std::exp((m - x) / s);
            double denom = 1.0 + aux;
            f1 += 1.0 / denom;
            denom *= denom;
            f2 -= aux / denom;
        }
        f1 -= nHalf;
        return DoublePair(f1, f2);
    }, root))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetLocation(root);
}


template class LogisticRand<float>;
template class LogisticRand<double>;
template class LogisticRand<long double>;
