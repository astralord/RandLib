#include "NoncentralChiSquaredRand.h"

NoncentralChiSquaredRand::NoncentralChiSquaredRand(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

std::string NoncentralChiSquaredRand::Name() const
{
    return "Noncentral Chi-Squared(" + toStringWithPrecision(GetDegree()) + ", "
            + toStringWithPrecision(GetNoncentrality()) + ")";
}

void NoncentralChiSquaredRand::SetParameters(double degree, double noncentrality)
{
    k = (degree <= 0.0) ? 1.0 : degree;
    halfK = 0.5 * k;

    lambda = (noncentrality <= 0.0) ? 1.0 : noncentrality;
    sqrtLambda = std::sqrt(lambda);
    logLambda = std::log(lambda);

    if (k < 1)
        Y.SetRate(0.5 * lambda);
}

double NoncentralChiSquaredRand::f(double x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (k == 2)
            return 0.5 * std::exp(-0.5 * lambda);
        return (k > 2) ? 0.0 : INFINITY;
    }
    double halfKm1 = halfK - 1;
    double y = RandMath::logModifiedBesselFirstKind(std::sqrt(lambda * x), halfKm1);
    double z = halfKm1 * (std::log(x) - logLambda);
    z -= x + lambda;
    return 0.5 * std::exp(y + 0.5 * z);
}

double NoncentralChiSquaredRand::F(double x) const
{
    return (x > 0) ? RandMath::MarcumP(halfK, 0.5 * lambda, 0.5 * x) : 0.0;
}

double NoncentralChiSquaredRand::S(double x) const
{
    return (x > 0) ? RandMath::MarcumQ(halfK, 0.5 * lambda, 0.5 * x) : 1.0;
}

double NoncentralChiSquaredRand::variateForDegreeEqualOne() const
{
    double y = sqrtLambda + NormalRand::StandardVariate();
    return y * y;
}

double NoncentralChiSquaredRand::Variate(double degree, double noncentrality)
{
    if (degree <= 0 || noncentrality < 0)
        return NAN; /// wrong parameters

    if (degree >= 1) {
        double rv = (degree == 1) ? 0.0 : 2 * GammaRand::StandardVariate(0.5 * degree - 0.5);
        double y = std::sqrt(noncentrality) + NormalRand::StandardVariate();
        return rv + y * y;
    }
    return 2 * GammaRand::StandardVariate(0.5 * degree + PoissonRand::Variate(0.5 * noncentrality));
}

double NoncentralChiSquaredRand::Variate() const
{
    if (k == 1)
        return variateForDegreeEqualOne();
    if (k > 1)
        return 2 * GammaRand::StandardVariate(halfK - 0.5) + variateForDegreeEqualOne();
    return 2 * GammaRand::StandardVariate(halfK + Y.Variate());
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
            var += 2 * GammaRand::StandardVariate(halfKmHalf);
    }
    else {
        for (double & var : outputData)
            var = 2 * GammaRand::StandardVariate(halfK + Y.Variate());
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
