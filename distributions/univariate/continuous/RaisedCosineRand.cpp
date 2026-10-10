#include "RaisedCosineRand.h"
#include "UniformRand.h"

template < typename RealType >
RaisedCosineDistribution<RealType>::RaisedCosineDistribution(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
void RaisedCosineDistribution<RealType>::SetLocation(double location)
{
    if (!std::isfinite(location))
        throw std::invalid_argument("Raised cosine distribution: location should be finite");
    mu = location;
}

template < typename RealType >
void RaisedCosineDistribution<RealType>::SetScale(double scale)
{
    if (!(scale > 0.0))
        throw std::invalid_argument("Raised-Cosine distribution: scale should be positive");
    s = scale;
    s_pi = s * randlib::consts::inv_pi;
    log2S = std::log(2 * s);
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::f(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -static_cast<double>(randlib::consts::pi) || xAdj >= static_cast<double>(randlib::consts::pi))
        return 0.0;
    double y = std::cos(xAdj) + 1.0;
    return 0.5 * y / s;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::logf(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -static_cast<double>(randlib::consts::pi) || xAdj >= static_cast<double>(randlib::consts::pi))
        return -INFINITY;
    double y = std::cos(xAdj);
    y = std::log1pl(y);
    return y - log2S;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::F(const RealType &x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 0.0;
    if (xAdj >= 1)
        return 1.0;
    double y = std::sin(xAdj * randlib::consts::pi);
    y *= randlib::consts::inv_pi;
    y += xAdj + 1;
    return 0.5 * y;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::S(const RealType &x) const
{
    double xAdj = (x - mu) / s;
    if (xAdj <= -1)
        return 1.0;
    if (xAdj >= 1)
        return 0.0;
    double y = std::sin(xAdj * randlib::consts::pi);
    y *= randlib::consts::inv_pi;
    y += xAdj;
    return 0.5 - 0.5 * y;
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    RealType X = randlib::consts::pi * UniformRand<RealType>::StandardVariate(randGenerator) - randlib::consts::half_pi;
    RealType XSq = X * X;
    RealType U = 2 * UniformRand<RealType>::StandardVariate(randGenerator);
    int a = 0, b = -1;
    RealType W = 0.0, V = 1.0;
    size_t iter = 0;
    do {
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W += V;
        if (U >= W)
            return X;
        a += 2;
        b += 2;
        V *= XSq / (a * b);
        W -= V;
        if (U <= W)
        {
            if (X == 0.0)
                return 0.0;
            return (X > 0.0) ? randlib::consts::pi - X : -randlib::consts::pi - X;
        }
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Raised-Cosine distribution: sampling failed");
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::Variate() const
{
    return mu + s_pi * StandardVariate(this->localRandGenerator);
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Mean() const
{
    return mu;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Variance() const
{
    static constexpr double coef = 1.0 / 3 - 2.0 / randlib::consts::pi_sq;
    return s * s * coef;
}

template < typename RealType >
std::complex<double> RaisedCosineDistribution<RealType>::CFImpl(double t) const
{
    double st = s * t;
    double numerator = randlib::consts::pi_sq * std::sin(st);
    double denominator = st * (randlib::consts::pi_sq - st * st);
    std::complex<double> y(0.0, mu * t);
    y = std::exp(y);
    return numerator / denominator * y;
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::Median() const
{
    return mu;
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::Skewness() const
{
    return 0.0l;
}

template < typename RealType >
long double RaisedCosineDistribution<RealType>::ExcessKurtosis() const
{
    static constexpr long double numerator = 1.2 * (90.0 - randlib::consts::pi_sq * randlib::consts::pi_sq);
    static constexpr long double denominator = randlib::consts::pi_sq - 6.0;
    static constexpr long double y = numerator / (denominator * denominator);
    return y;
}

template class RaisedCosineDistribution<float>;
template class RaisedCosineDistribution<double>;
template class RaisedCosineDistribution<long double>;

template < typename RealType >
String RaisedCosineRand<RealType>::Name() const
{
    return "Raised cosine(" + this->toStringWithPrecision(this->GetLocation()) + ", " + this->toStringWithPrecision(this->GetScale()) + ")";
}

template class RaisedCosineRand<float>;
template class RaisedCosineRand<double>;
template class RaisedCosineRand<long double>;

template < typename RealType >
String RaabGreenRand<RealType>::Name() const
{
    return "Raab Green";
}

template class RaabGreenRand<float>;
template class RaabGreenRand<double>;
template class RaabGreenRand<long double>;
