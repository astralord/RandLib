#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/NormalRand.h"

template< typename RealType >
NormalInverseGammaRand<RealType>::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    SetParameters(location, precision, shape, rate);
}

template< typename RealType >
String NormalInverseGammaRand<RealType>::Name() const
{
    return "Normal-Inverse-Gamma(" + this->toStringWithPrecision(GetLocation()) + ", "
                                   + this->toStringWithPrecision(GetPrecision()) + ", "
                                   + this->toStringWithPrecision(GetShape()) + ", "
                                   + this->toStringWithPrecision(GetRate()) + ")";
}

template< typename RealType >
void NormalInverseGammaRand<RealType>::SetParameters(double location, double precision, double shape, double rate)
{
    if (precision <= 0.0)
        throw std::invalid_argument("Precision of Normal-Inverse-Gamma distribution should be positive");
    if (shape <= 0.0)
        throw std::invalid_argument("Shape of Normal-Inverse-Gamma distribution should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Rate of Normal-Inverse-Gamma distribution should be positive");

    mu = location;
    lambda = precision;

    this->Y.SetParameters(shape, rate);
    alpha = this->Y.GetShape();
    beta = this->Y.GetRate();
    this->X.SetDegree(2 * alpha);
    this->X.SetLocation(mu);
    this->X.SetScale(std::sqrt(alpha * lambda / beta));

    pdfCoef = 0.5 * std::log(0.5 * lambda / M_PI);
    pdfCoef += alpha * this->Y.GetLogRate() - this->Y.GetLogGammaShape();
}

template< typename RealType >
double NormalInverseGammaRand<RealType>::f(const Pair<RealType> &point) const
{
    return (point.second > 0.0) ? std::exp(logf(point)) : 0.0;
}

template< typename RealType >
double NormalInverseGammaRand<RealType>::logf(const Pair<RealType> &point) const
{
    double sigmaSq = point.second;
    if (sigmaSq <= 0)
        return -INFINITY;
    double x = point.first;
    double y = (alpha + 1.5) * std::log(sigmaSq);
    double degree = x - mu;
    degree *= degree;
    degree *= lambda;
    degree += 2 * beta;
    degree *= 0.5 / sigmaSq;
    y += degree;
    return pdfCoef - y;
}

template< typename RealType >
double NormalInverseGammaRand<RealType>::F(const Pair<RealType> &point) const
{
    double sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double x = point.first;
    double y = 0.5 * lambda;
    double xmmu = x - mu;
    y *= xmmu * xmmu / sigmaSq;
    y = std::erfc(-std::sqrt(y));
    double z = beta / sigmaSq;
    double temp = alpha * std::log(z) - z;
    y *= std::exp(temp - this->Y.GetLogGammaShape());
    y *= 0.5 / sigmaSq;
    return y;
}

template< typename RealType >
Pair<RealType> NormalInverseGammaRand<RealType>::Variate() const
{
    Pair<RealType> var;
    var.second = this->Y.Variate();
    double coef = std::sqrt(var.second) / lambda;
    var.first = mu + coef * NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    return var;
}

template< typename RealType >
long double NormalInverseGammaRand<RealType>::Correlation() const
{
    return 0.0;
}

template< typename RealType >
Pair<RealType> NormalInverseGammaRand<RealType>::Mode() const
{
    return std::make_pair(mu, 2 * beta / (2 * alpha + 3));
}

template class NormalInverseGammaRand<float>;
template class NormalInverseGammaRand<double>;
template class NormalInverseGammaRand<long double>;
