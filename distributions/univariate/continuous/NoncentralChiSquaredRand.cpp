#include "NoncentralChiSquaredRand.h"

template < typename RealType >
NoncentralChiSquaredRand<RealType>::NoncentralChiSquaredRand(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

template < typename RealType >
String NoncentralChiSquaredRand<RealType>::Name() const
{
    return "Noncentral Chi-Squared(" + this->toStringWithPrecision(GetDegree()) + ", "
            + this->toStringWithPrecision(GetNoncentrality()) + ")";
}

template < typename RealType >
void NoncentralChiSquaredRand<RealType>::SetParameters(double degree, double noncentrality)
{
    if (degree <= 0.0)
        throw std::invalid_argument("Noncentral Chi-Squared distribution: degree parameter should be positive");
    if (noncentrality <= 0.0)
        throw std::invalid_argument("Noncentral Chi-Squared distribution: noncentrality parameter should be positive");

    k = degree;
    halfK = 0.5 * k;

    lambda = noncentrality;
    halfLambda = 0.5 * lambda;
    sqrtLambda = std::sqrt(lambda);
    logLambda = std::log(lambda);

    if (k < 1)
        Y.SetRate(halfLambda);
}

template < typename RealType >
double NoncentralChiSquaredRand<RealType>::f(const RealType & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (k == 2)
            return 0.5 * std::exp(-halfLambda);
        return (k > 2) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

template < typename RealType >
double NoncentralChiSquaredRand<RealType>::logf(const RealType & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (k == 2)
            return -halfLambda - M_LN2;
        return (k > 2) ? -INFINITY : INFINITY;
    }
    double halfKm1 = halfK - 1;
    double y = RandMath::logBesselI(halfKm1, std::sqrt(lambda * x));
    double z = halfKm1 * (std::log(x) - logLambda);
    z -= x + lambda;
    return y + 0.5 * z - M_LN2;
}

template < typename RealType >
double NoncentralChiSquaredRand<RealType>::F(const RealType &x) const
{
    if (x <= 0.0)
        return 0.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return RandMath::MarcumP(halfK, halfLambda, halfX, sqrtLambda / M_SQRT2, sqrtHalfX, logLambda - M_LN2, logHalfX);
}

template < typename RealType >
double NoncentralChiSquaredRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return RandMath::MarcumQ(halfK, halfLambda, halfX, sqrtLambda / M_SQRT2, sqrtHalfX, logLambda - M_LN2, logHalfX);
}

template < typename RealType >
RealType NoncentralChiSquaredRand<RealType>::variateForDegreeEqualOne() const
{
    RealType y = sqrtLambda + NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    return y * y;
}

template < typename RealType >
RealType NoncentralChiSquaredRand<RealType>::Variate(double degree, double noncentrality, RandGenerator &randGenerator)
{
    if (degree <= 0.0)
        throw std::invalid_argument("Noncentral Chi-Squared distribution: degree parameter should be positive");
    if (noncentrality <= 0.0)
        throw std::invalid_argument("Noncentral Chi-Squared distribution: noncentrality parameter should be positive");

    if (degree >= 1) {
        RealType rv = (degree == 1) ? 0.0 : 2 * GammaDistribution<RealType>::StandardVariate(0.5 * degree - 0.5, randGenerator);
        RealType y = std::sqrt(noncentrality) + NormalRand<RealType>::StandardVariate(randGenerator);
        return rv + y * y;
    }
    RealType shape = 0.5 * degree + PoissonRand<int>::Variate(0.5 * noncentrality, randGenerator);
    return 2 * GammaDistribution<RealType>::StandardVariate(shape, randGenerator);
}

template < typename RealType >
RealType NoncentralChiSquaredRand<RealType>::Variate() const
{
    if (k < 1)
        return 2 * GammaDistribution<RealType>::StandardVariate(halfK + Y.Variate(), this->localRandGenerator);
    double X = variateForDegreeEqualOne();
    if (k > 1)
        X += 2 * GammaDistribution<RealType>::StandardVariate(halfK - 0.5, this->localRandGenerator);
    return X;
}

template < typename RealType >
void NoncentralChiSquaredRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    if (k >= 1) {
        for (RealType & var : outputData)
            var = variateForDegreeEqualOne();
        double halfKmHalf = halfK - 0.5;
        if (halfKmHalf == 0)
            return;
        for (RealType & var : outputData)
            var += 2 * GammaDistribution<RealType>::StandardVariate(halfKmHalf, this->localRandGenerator);
    }
    else {
        for (RealType & var : outputData)
            var = 2 * GammaDistribution<RealType>::StandardVariate(halfK + Y.Variate(), this->localRandGenerator);
    }
}

template < typename RealType >
void NoncentralChiSquaredRand<RealType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    Y.Reseed(seed + 1);
}

template < typename RealType >
long double NoncentralChiSquaredRand<RealType>::Mean() const
{
    return k + lambda;
}

template < typename RealType >
long double NoncentralChiSquaredRand<RealType>::Variance() const
{
    return 2 * (k + 2 * lambda);
}

template < typename RealType >
RealType NoncentralChiSquaredRand<RealType>::Mode() const
{
    return (k <= 2) ? 0.0 : ContinuousDistribution<RealType>::Mode();
}

template < typename RealType >
long double NoncentralChiSquaredRand<RealType>::Skewness() const
{
    long double y = k + 2 * lambda;
    y = 2.0 / y;
    long double z = y * std::sqrt(y);
    return z * (k + 3 * lambda);
}

template < typename RealType >
long double NoncentralChiSquaredRand<RealType>::ExcessKurtosis() const
{
    long double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}

template < typename RealType >
std::complex<double> NoncentralChiSquaredRand<RealType>::CFImpl(double t) const
{
    std::complex<double> aux(1, -2 * t);
    std::complex<double> y(0, lambda * t);
    y /= aux;
    y -= halfK * std::log(aux);
    return std::exp(y);
}

template class NoncentralChiSquaredRand<float>;
template class NoncentralChiSquaredRand<double>;
template class NoncentralChiSquaredRand<long double>;
