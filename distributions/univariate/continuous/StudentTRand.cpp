#include "StudentTRand.h"
#include "NormalRand.h"
#include "CauchyRand.h"

template < typename RealType >
StudentTRand<RealType>::StudentTRand(double degree, double location, double scale)
{
    SetDegree(degree);
    SetLocation(location);
    SetScale(scale);
}

template < typename RealType >
String StudentTRand<RealType>::Name() const
{
    if (mu == 0.0 && sigma == 1.0)
        return "Student-t(" + this->toStringWithPrecision(GetDegree()) + ")";
    return "Student-t(" + this->toStringWithPrecision(GetDegree()) + ", "
                          + this->toStringWithPrecision(GetLocation()) + ", "
                          + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void StudentTRand<RealType>::SetDegree(double degree)
{
    if (degree <= 0.0)
        throw std::invalid_argument("Student-t distribution: degree parameter should be positive");
    nu = degree;
    Y.SetParameters(0.5 * nu, 1.0);

    nup1Half = 0.5 * (nu + 1);
    pdfCoef = Y.GetLogGammaShapeRatio();
    pdfCoef -= 0.5 * M_LNPI;
    logBetaFun = -pdfCoef;
    pdfCoef -= 0.5 * std::log(nu);
}

template < typename RealType >
void StudentTRand<RealType>::SetLocation(double location)
{
    mu = location;
}

template < typename RealType >
void StudentTRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Student-t distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

template < typename RealType >
double StudentTRand<RealType>::f(const RealType & x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double xSq = x0 * x0;
    if (nu == 1) /// Cauchy distribution
        return M_1_PI / (sigma * (1 + xSq));
    if (nu == 2)
        return 1.0 / (sigma * std::pow(2.0 + xSq, 1.5));
    if (nu == 3) {
        double y = 3 + xSq;
        return 6 * M_SQRT3 * M_1_PI / (sigma * y * y);
    }
    return std::exp(logf(x));
}

template < typename RealType >
double StudentTRand<RealType>::logf(const RealType & x) const
{
    /// adjustment
    double x0 = x - mu;
    x0 /= sigma;
    double xSq = x0 * x0;
    if (nu == 1) /// Cauchy distribution
        return -logSigma - M_LNPI - std::log1pl(xSq);
    if (nu == 2)
        return -logSigma - 1.5 * std::log(2.0 + xSq);
    if (nu == 3) {
        double y = -2 * std::log(3.0 + xSq);
        y += M_LN2 + 1.5 * M_LN3 - M_LNPI - logSigma;
        return y;
    }
    double y = -nup1Half * std::log1pl(xSq / nu);
    return pdfCoef + y - logSigma;
}

template < typename RealType >
double StudentTRand<RealType>::F(const RealType & x) const
{
    double x0 = x - mu;
    x0 /= sigma;
    if (x0 == 0.0)
        return 0.5;
    if (nu == 1) {
        /// Cauchy distribution
        return 0.5 + M_1_PI * RandMath::atan(x0);
    }
    double xSq = x0 * x0;
    if (nu == 2) {
        return 0.5 + 0.5 * x0 / std::sqrt(2 + xSq);
    }
    if (nu == 3) {
        double y = M_SQRT3 * x0 / (xSq + 3);
        y += RandMath::atan(x0 / M_SQRT3);
        return 0.5 + M_1_PI * y;
    }
    double t = nu / (xSq + nu);
    double y = 0.5 * RandMath::ibeta(t, 0.5 * nu, 0.5, logBetaFun, std::log(t), std::log1pl(-t));
    return (x0 > 0.0) ? (1.0 - y) : y;
}

template < typename RealType >
double StudentTRand<RealType>::S(const RealType & x) const
{
    double x0 = x - mu;
    x0 /= sigma;
    if (x0 == 0.0)
        return 0.5;
    if (nu == 1) {
        /// Cauchy distribution
        return 0.5 + M_1_PI * RandMath::atan(-x0);
    }
    double xSq = x0 * x0;
    if (nu == 2) {
        return 0.5 - 0.5 * x0 / std::sqrt(2 + xSq);
    }
    if (nu == 3) {
        double y = M_SQRT3 * x0 / (xSq + 3);
        y += RandMath::atan(x0 / M_SQRT3);
        return 0.5 - M_1_PI * y;
    }
    double t = nu / (xSq + nu);
    double y = 0.5 * RandMath::ibeta(t, 0.5 * nu, 0.5, logBetaFun, std::log(t), std::log1pl(-t));
    return (x0 > 0.0) ? y : 1.0 - y;
}

template < typename RealType >
RealType StudentTRand<RealType>::Variate() const
{
    if (nu == 1)
        return mu + sigma * CauchyRand<RealType>::StandardVariate(this->localRandGenerator);
    return mu + sigma * NormalRand<RealType>::StandardVariate(this->localRandGenerator) / Y.Variate();
}

template < typename RealType >
void StudentTRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    if (nu == 1) {
        for (RealType &var : outputData)
            var = mu + sigma * CauchyRand<RealType>::StandardVariate(this->localRandGenerator);
    }
    else {
        Y.Sample(outputData);
        for (RealType &var : outputData)
            var = mu + sigma * NormalRand<RealType>::StandardVariate(this->localRandGenerator) / var;
    }
}

template < typename RealType >
void StudentTRand<RealType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    Y.Reseed(seed + 1);
}

template < typename RealType >
long double StudentTRand<RealType>::Mean() const
{
    return (nu > 1) ? mu : NAN;
}

template < typename RealType >
long double StudentTRand<RealType>::Variance() const
{
    if (nu > 2)
        return sigma * sigma * nu / (nu - 2);
    return (nu > 1) ? INFINITY : NAN;
}

template < typename RealType >
std::complex<double> StudentTRand<RealType>::CFImpl(double t) const
{
    double x = std::sqrt(nu) * t * sigma; // value of sqrt(nu) can be hashed
    double vHalf = 0.5 * nu;
    double y = vHalf * std::log(x);
    y -= Y.GetLogGammaFunction();
    y -= (vHalf - 1) * M_LN2;
    y += RandMath::logBesselK(vHalf, x);
    double costmu = std::cos(t * mu), sintmu = std::sin(t * mu);
    std::complex<double> cf(costmu, sintmu);
    return std::exp(y) * cf;
}

template < typename RealType >
RealType StudentTRand<RealType>::quantileImpl(double p) const
{
    double temp = p - 0.5;
    if (nu == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2.0 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * RandMath::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution<RealType>::quantileImpl(p);
}

template < typename RealType >
RealType StudentTRand<RealType>::quantileImpl1m(double p) const
{
    double temp = 0.5 - p;
    if (nu == 1)
        return std::tan(M_PI * temp) * sigma + mu;
    double pq = p * (1.0 - p);
    if (nu == 2)
        return sigma * 2 * temp * std::sqrt(0.5 / pq) + mu;
    if (nu == 4)
    {
        double alpha = 2 * std::sqrt(pq);
        double beta = std::cos(std::acos(alpha) / 3.0) / alpha - 1;
        return mu + sigma * 2 * RandMath::sign(temp) * std::sqrt(beta);
    }
    return ContinuousDistribution<RealType>::quantileImpl1m(p);
}

template < typename RealType >
RealType StudentTRand<RealType>::Median() const
{
    return mu;
}

template < typename RealType >
RealType StudentTRand<RealType>::Mode() const
{
    return mu;
}

template < typename RealType >
long double StudentTRand<RealType>::Skewness() const
{
    return (nu > 3) ? 0.0 : NAN;
}

template < typename RealType >
long double StudentTRand<RealType>::ExcessKurtosis() const
{
    if (nu > 4)
        return 6.0 / (nu - 4);
    return (nu > 2) ? INFINITY : NAN;
}
