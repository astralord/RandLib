#include "NoncentralChiSquaredRand.h"

NoncentralChiSquaredRand::NoncentralChiSquaredRand(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

String NoncentralChiSquaredRand::Name() const
{
    return "Noncentral Chi-Squared(" + toStringWithPrecision(GetDegree()) + ", "
            + toStringWithPrecision(GetNoncentrality()) + ")";
}

void NoncentralChiSquaredRand::SetParameters(double degree, double noncentrality)
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

double NoncentralChiSquaredRand::f(const double & x) const
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

double NoncentralChiSquaredRand::logf(const double & x) const
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

double NoncentralChiSquaredRand::F(const double & x) const
{
    if (x <= 0.0)
        return 0.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return RandMath::MarcumP(halfK, halfLambda, halfX, sqrtLambda / M_SQRT2, sqrtHalfX, logLambda - M_LN2, logHalfX);
}

double NoncentralChiSquaredRand::S(const double & x) const
{
    if (x <= 0.0)
        return 1.0;
    double halfX = 0.5 * x;
    double sqrtHalfX = std::sqrt(halfX);
    double logHalfX = std::log(halfX);
    return RandMath::MarcumQ(halfK, halfLambda, halfX, sqrtLambda / M_SQRT2, sqrtHalfX, logLambda - M_LN2, logHalfX);
}

double NoncentralChiSquaredRand::variateForDegreeEqualOne() const
{
    double y = sqrtLambda + NormalRand::StandardVariate();
    return y * y;
}

double NoncentralChiSquaredRand::Variate(double degree, double noncentrality)
{
    if (degree <= 0 || noncentrality <= 0)
        return NAN;

    if (degree >= 1) {
        double rv = (degree == 1) ? 0.0 : 2 * GammaDistribution::StandardVariate(0.5 * degree - 0.5);
        double y = std::sqrt(noncentrality) + NormalRand::StandardVariate();
        return rv + y * y;
    }
    double shape = 0.5 * degree + PoissonRand::Variate(0.5 * noncentrality);
    return 2 * GammaDistribution::StandardVariate(shape);
}

double NoncentralChiSquaredRand::Variate() const
{
    if (k < 1)
        return 2 * GammaDistribution::StandardVariate(halfK + Y.Variate());
    double X = variateForDegreeEqualOne();
    if (k > 1)
        X += 2 * GammaDistribution::StandardVariate(halfK - 0.5);
    return X;
}

void NoncentralChiSquaredRand::Sample(std::vector<double> &outputData) const
{
    if (k >= 1) {
        for (double & var : outputData)
            var = variateForDegreeEqualOne();
        if (RandMath::areClose(k, 1))
            return;
        double halfKmHalf = halfK - 0.5;
        for (double & var : outputData)
            var += 2 * GammaDistribution::StandardVariate(halfKmHalf);
    }
    else {
        for (double & var : outputData)
            var = 2 * GammaDistribution::StandardVariate(halfK + Y.Variate());
    }
}

double NoncentralChiSquaredRand::Mean() const
{
    return k + lambda;
}

double NoncentralChiSquaredRand::Variance() const
{
    return 2 * (k + 2 * lambda);
}

double NoncentralChiSquaredRand::Mode() const
{
    return (k < 2) ? 0.0 : ContinuousDistribution::Mode();
}

std::complex<double> NoncentralChiSquaredRand::CFImpl(double t) const
{
    std::complex<double> aux(1, -2 * t);
    std::complex<double> y(0, lambda * t);
    y /= aux;
    y -= halfK * std::log(aux);
    return std::exp(y);
}

double NoncentralChiSquaredRand::Skewness() const
{
    double y = k + 2 * lambda;
    y = 2.0 / y;
    double z = y * std::sqrt(y);
    return z * (k + 3 * lambda);
}

double NoncentralChiSquaredRand::ExcessKurtosis() const
{
    double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}
