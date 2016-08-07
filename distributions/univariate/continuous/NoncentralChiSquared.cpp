#include "NoncentralChiSquared.h"

NoncentralChiSquared::NoncentralChiSquared(double degree, double noncentrality)
{
    setParameters(degree, noncentrality);
}

std::string NoncentralChiSquared::name() const
{
    return "Noncentral Chi-Squared(" + toStringWithPrecision(getDegree()) + ", "
            + toStringWithPrecision(getNoncentrality()) + ")";
}

void NoncentralChiSquared::setParameters(double degree, double noncentrality)
{
    k = degree;
    if (k <= 0)
        k = 1;
    halfK = 0.5 * k;

    lambda = noncentrality;
    if (lambda <= 0)
        lambda = 1;
    sqrtLambda = std::sqrt(lambda);
    logLambda = std::log(lambda);

    if (k < 1)
        Y.setRate(0.5 * lambda);

    if (k < 2) {
        cdfCoef = lambda + k * M_LN2;
        cdfCoef *= 0.5;
        cdfCoef += std::lgamma(halfK);
    }
}

double NoncentralChiSquared::f(double x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (k == 2)
            return 0.5 * std::exp(-0.5 * lambda);
        return (k > 2) ? 0.0 : INFINITY;
    }
    double halfKm1 = halfK - 1;
    double y = RandMath::modifiedBesselFirstKind(std::sqrt(lambda * x), halfKm1);
    double z = halfKm1 * (std::log(x) - logLambda);
    z -= x + lambda;
    y = 0.5 * z + std::log(y);
    return 0.5 * std::exp(y);
}

double NoncentralChiSquared::F(double x) const
{
    if (x <= 0)
        return 0.0;
    if (k >= 2) {
        return RandMath::integral([this] (double t)
        {
            return f(t);
        }, 0, x);
    }

    /// in this case we have singularity point at 0,
    /// so we get rid of it by subtracting the function
    /// which has the same behaviour at this point
    double y = std::log(x) * halfK;
    y -= cdfCoef;
    y = std::exp(y) / halfK;

    double halfKm1 = halfK - 1.0;
    y += RandMath::integral([this, halfKm1] (double t)
    {
        if (t <= 0)
            return 0.0;
        /// Calculate log of leveling factor
        double logT = std::log(t);
        double exponent = halfKm1 * logT;
        exponent -= cdfCoef;

        /// Calculate log(2f(t))
        double bessel = RandMath::modifiedBesselFirstKind(std::sqrt(lambda * t), halfKm1);
        double z = halfKm1 * (logT - logLambda);
        z -= t + lambda;
        double log2F = 0.5 * z + std::log(bessel);

        /// Return difference f(t) - factor
        return 0.5 * std::exp(log2F) - std::exp(exponent);
    }, 0, x);

    return y;
}

double NoncentralChiSquared::variateForDegreeEqualOne() const
{
    double y = sqrtLambda + NormalRand::standardVariate();
    return y * y;
}

double NoncentralChiSquared::variate(double degree, double noncentrality)
{
    if (degree <= 0 || noncentrality < 0)
        return NAN; /// wrong parameters

    if (degree >= 1) {
        double rv = (degree == 1) ? 0.0 : 2 * GammaRand::standardVariate(0.5 * degree - 0.5);
        double y = std::sqrt(noncentrality) + NormalRand::standardVariate();
        return rv + y * y;
    }
    return 2 * GammaRand::standardVariate(0.5 * degree + PoissonRand::variate(0.5 * noncentrality));
}

double NoncentralChiSquared::variate() const
{
    if (k == 1)
        return variateForDegreeEqualOne();
    if (k > 1)
        return 2 * GammaRand::standardVariate(halfK - 0.5) + variateForDegreeEqualOne();
    return 2 * GammaRand::standardVariate(halfK + Y.variate());
}

void NoncentralChiSquared::sample(std::vector<double> &outputData) const
{
    if (k >= 1) {
        for (double & var : outputData)
            var = variateForDegreeEqualOne();
        if (RandMath::areClose(k, 1))
            return;
        double halfKmHalf = halfK - 0.5;
        for (double & var : outputData)
            var += 2 * GammaRand::standardVariate(halfKmHalf);
    }
    else {
        for (double & var : outputData)
            var = 2 * GammaRand::standardVariate(halfK + Y.variate());
    }
}

double NoncentralChiSquared::Mean() const
{
    return k + lambda;
}

double NoncentralChiSquared::Variance() const
{
    return 2 * (k + 2 * lambda);
}

std::complex<double> NoncentralChiSquared::CF(double t) const
{
    if (t == 0)
        return 1;
    std::complex<double> aux(1, -2 * t);
    std::complex<double> y(0, lambda * t);
    y = std::exp(y / aux);
    return y / std::pow(aux, halfK);
}

double NoncentralChiSquared::Skewness() const
{
    double y = k + 2 * lambda;
    y = std::pow(2.0 / y, 1.5);
    return y * (k + 3 * lambda);
}

double NoncentralChiSquared::ExcessKurtosis() const
{
    double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}
