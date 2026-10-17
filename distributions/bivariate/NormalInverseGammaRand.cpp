#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/NormalRand.h"
#include "../../math/NumericMath.h"

#include <algorithm>

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
    if (!std::isfinite(location))
        throw std::invalid_argument("Normal-Inverse-Gamma distribution: location should be finite");
    if (!(precision > 0.0))
        throw std::invalid_argument("Precision of Normal-Inverse-Gamma distribution should be positive, but it's equal to "
                                    + std::to_string(precision));
    if (!(shape > 0.0))
        throw std::invalid_argument("Shape of Normal-Inverse-Gamma distribution should be positive, but it's equal to "
                                    + std::to_string(shape));
    if (!(rate > 0.0))
        throw std::invalid_argument("Rate of Normal-Inverse-Gamma distribution should be positive, but it's equal to "
                                    + std::to_string(rate));

    mu = location;
    lambda = precision;

    this->Y.SetParameters(shape, rate);
    alpha = this->Y.GetShape();
    beta = this->Y.GetRate();
    this->X.SetDegree(2 * alpha);
    this->X.SetLocation(mu);
    this->X.SetScale(std::sqrt(beta / (alpha * lambda)));

    pdfCoef = 0.5 * std::log(0.5 * lambda / randlib::consts::pi);
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
    RealType sigmaSq = point.second;
    if (sigmaSq <= 0)
        return -INFINITY;
    RealType x = point.first;
    RealType y = (alpha + 1.5) * std::log(sigmaSq);
    RealType degree = x - mu;
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
    RealType sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    RealType x = point.first;
    if (x == -INFINITY)
        return 0.0;
    double marginalY = this->Y.F(sigmaSq);
    if (x == INFINITY)
        return marginalY;
    if (sigmaSq == INFINITY)
        return this->X.F(x);
    if (x == mu)
        return 0.5 * marginalY;

    double lower = beta / sigmaSq;
    double normalCoef = (x - mu) * std::sqrt(lambda / beta);
    double logGammaAlpha = this->Y.GetLogGammaShape();
    double joint = RandMath::integral([=](double t) {
        double oneMinusT = 1.0 - t;
        double u = lower + t / oneMinusT;
        double logWeight = (alpha - 1.0) * std::log(u) - u - logGammaAlpha;
        logWeight -= 2.0 * std::log(oneMinusT);
        double conditional = 0.5 * std::erfc(-normalCoef * std::sqrt(0.5 * u));
        return std::exp(logWeight) * conditional;
    }, 0.0, 1.0);
    return std::clamp(joint, 0.0, std::min(this->X.F(x), marginalY));
}

template< typename RealType >
Pair<RealType> NormalInverseGammaRand<RealType>::Variate() const
{
    Pair<RealType> var;
    var.second = this->Y.Variate();
    double coef = std::sqrt(var.second / lambda);
    var.first = mu + coef * NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    return var;
}

template< typename RealType >
long double NormalInverseGammaRand<RealType>::Correlation() const
{
    return (alpha > 2.0) ? 0.0 : NAN;
}

template< typename RealType >
Pair<RealType> NormalInverseGammaRand<RealType>::Mode() const
{
    return std::make_pair(mu, 2 * beta / (2 * alpha + 3));
}

template class NormalInverseGammaRand<float>;
template class NormalInverseGammaRand<double>;
template class NormalInverseGammaRand<long double>;
