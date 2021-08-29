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
    mu = location;
}

template < typename RealType >
void RaisedCosineDistribution<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Raised-Cosine distribution: scale should be positive");
    s = scale;
    s_pi = s * M_1_PI;
    log2S = std::log(2 * s);
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::f(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -M_PI || xAdj >= M_PI)
        return 0.0;
    double y = std::cos(xAdj) + 1.0;
    return 0.5 * y / s;
}

template < typename RealType >
double RaisedCosineDistribution<RealType>::logf(const RealType &x) const
{
    double xAdj = (x - mu) / s_pi;
    if (xAdj <= -M_PI || xAdj >= M_PI)
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
    double y = std::sin(xAdj * M_PI);
    y *= M_1_PI;
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
    double y = std::sin(xAdj * M_PI);
    y *= M_1_PI;
    y += xAdj;
    return 0.5 - 0.5 * y;
}

template < typename RealType >
RealType RaisedCosineDistribution<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    /// p. 160. Non-Uniform Random Variate Generation. Luc Devroye
    RealType X = M_PI * UniformRand<RealType>::StandardVariate(randGenerator) - M_PI_2;
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
            return (X > 0.0) ? M_PI - X : -M_PI - X;
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
    static constexpr double coef = 1.0 / 3 - 2.0 / M_PI_SQ;
    return s * s * coef;
}

template < typename RealType >
std::complex<double> RaisedCosineDistribution<RealType>::CFImpl(double t) const
{
    double st = s * t;
    double numerator = M_PI_SQ * std::sin(st);
    double denominator = st * (M_PI_SQ - st * st);
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
    static constexpr long double numerator = 1.2 * (90.0 - M_PI_SQ * M_PI_SQ);
    static constexpr long double denominator = M_PI_SQ - 6.0;
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
