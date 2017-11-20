#include "BivariateNormalRand.h"


template< typename RealType >
BivariateNormalRand<RealType>::BivariateNormalRand(double location1, double location2, double scale1, double scale2, double correlation)
{
    SetLocations(location1, location2);
    SetCovariance(scale1, scale2, correlation);
}

template< typename RealType >
String BivariateNormalRand<RealType>::Name() const
{
    return "Bivariate Normal( (" + this->toStringWithPrecision(GetFirstLocation()) + ", "
                               + this->toStringWithPrecision(GetSecondLocation()) + "), ("
                               + this->toStringWithPrecision(GetFirstScale()) + ", "
                               + this->toStringWithPrecision(GetSecondScale()) + ", "
                               + this->toStringWithPrecision(GetCorrelation()) + ") )";
}

template< typename RealType >
void BivariateNormalRand<RealType>::SetLocations(double location1, double location2)
{
    mu1 = location1;
    mu2 = location2;
    this->X.SetLocation(mu1);
    this->Y.SetLocation(mu2);
}

template< typename RealType >
void BivariateNormalRand<RealType>::SetCovariance(double scale1, double scale2, double correlation)
{
    if (scale1 <= 0.0 || scale2 <= 0.0)
        throw std::invalid_argument("Scales of Bivariate-Normal distribution should be positive");
    if (std::fabs(correlation) > 1.0)
        throw std::invalid_argument("Correlation of Bivariate-Normal distribution should be in interval [-1, 1]");
    sigma1 = scale1;
    sigma2 = scale2;
    rho = correlation;

    this->X.SetScale(sigma1);
    this->Y.SetScale(sigma2);

    sqrt1mroSq = 0.5 * std::log1pl(-rho * rho);
    pdfCoef = this->X.GetLogScale() + this->Y.GetLogScale() + sqrt1mroSq + M_LN2 + M_LNPI;
    sqrt1mroSq = std::exp(sqrt1mroSq);
}

template< typename RealType >
double BivariateNormalRand<RealType>::f(const Pair<RealType> &point) const
{
    RealType xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    if (rho == 1.0) /// f(x, y) = δ(xAdj - yAdj)
        return (xAdj - yAdj == 0) ? INFINITY : 0.0;
    if (rho == -1.0) /// f(x, y) = δ(xAdj + yAdj)
        return (xAdj + yAdj == 0) ? INFINITY : 0.0;
    return std::exp(logf(point));
}

template< typename RealType >
double BivariateNormalRand<RealType>::logf(const Pair<RealType> &point) const
{
    if (rho == 0.0) /// log(f(x, y)) = log(f(x)) + log(f(y))
        return this->X.logf(point.first) + this->Y.logf(point.second);
    RealType xAdj = (point.first - mu1) / sigma1, yAdj = (point.second - mu2) / sigma2;
    if (rho == 1.0) /// log(f(x, y)) = log(δ(xAdj - yAdj))
        return (xAdj - yAdj == 0) ? INFINITY : -INFINITY;
    if (rho == -1.0) /// log(f(x, y)) = log(δ(xAdj + yAdj))
        return (xAdj + yAdj == 0) ? INFINITY : -INFINITY;
    RealType z = xAdj * xAdj - 2 * rho * xAdj * yAdj + yAdj * yAdj;
    return -(pdfCoef + 0.5 * z / (sqrt1mroSq * sqrt1mroSq));
}

template< typename RealType >
double BivariateNormalRand<RealType>::F(const Pair<RealType> &point) const
{
    if (rho == 0.0) /// F(x, y) = F(x)F(y)
        return this->X.F(point.first) * this->Y.F(point.second);
    /// Unnumbered equation between (3) and (4) in Section 2.2 of Genz (2004),
    /// integrating in terms of theta between asin(ρ) and +/- π/2
    RealType p1, p2 = 0.0;
    RealType x = point.first, y = point.second;
    RealType xAdj = (x - mu1) / sigma1, yAdj = (y - mu2) / sigma2;
    if (rho > 0)
        p1 = (xAdj < yAdj) ? this->X.F(x) : this->Y.F(y);
    else
        p1 = std::max(this->X.F(x) - this->Y.F(-y), 0.0);
    if (std::fabs(rho) < 1.0) {
        RealType lowLimit = std::asin(rho);
        RealType highLimit = RandMath::sign(rho) * M_PI_2;
        p2 = RandMath::integral([this, xAdj, yAdj] (double theta) {
            /// Integrand is exp(-(x^2 + y^2 - 2xysin(θ)) / (2cos(θ)^2))
            RealType cosTheta = std::cos(theta);
            RealType tanTheta = std::tan(theta);
            RealType integrand = xAdj * tanTheta - yAdj / cosTheta;
            integrand *= integrand;
            integrand += xAdj * xAdj;
            integrand = std::exp(-0.5 * integrand);
            return integrand;
        }, lowLimit, highLimit);
    }
    return p1 - 0.5 * p2 / M_PI;
}

template< typename RealType >
Pair<RealType> BivariateNormalRand<RealType>::Variate() const
{
    RealType Z1 = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType Z2 = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    RealType x = mu1 + sigma1 * Z1;
    RealType y = mu2 + sigma2 * (rho * Z1 + sqrt1mroSq * Z2);
    return std::make_pair(x, y);
}

template< typename RealType >
long double BivariateNormalRand<RealType>::Correlation() const
{
    return rho;
}

template< typename RealType >
Pair<RealType> BivariateNormalRand<RealType>::Mode() const
{
    return std::make_pair(mu1, mu2);
}

template class BivariateNormalRand<float>;
template class BivariateNormalRand<double>;
template class BivariateNormalRand<long double>;
