#include "NoncentralChiSquared.h"

NoncentralChiSquared::NoncentralChiSquared(double degree, double noncentrality)
{
    SetParameters(degree, noncentrality);
}

std::string NoncentralChiSquared::Name() const
{
    return "Noncentral Chi-Squared(" + toStringWithPrecision(GetDegree()) + ", "
            + toStringWithPrecision(GetNoncentrality()) + ")";
}

void NoncentralChiSquared::SetParameters(double degree, double noncentrality)
{
    k = (degree <= 0.0) ? 1.0 : degree;
    halfK = 0.5 * k;

    lambda = (noncentrality <= 0.0) ? 1.0 : noncentrality;
    sqrtLambda = std::sqrt(lambda);
    logLambda = std::log(lambda);

    if (k < 1)
        Y.SetRate(0.5 * lambda);

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

    /// for k == 2 https://pdfs.semanticscholar.org/e410/e73d7a92d3869205e347add025f1bda666ac.pdf
    /// also look http://www.tlc.unipr.it/ferrari/Publications/Journals/CoFe02.pdf

    if (x == lambda && k == 2) {
        double y = RandMath::modifiedBesselFirstKind(x, 0);
        y *= std::exp(-x);
        return 0.5 * (1 - y);
    }

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
    double y = sqrtLambda + NormalRand::StandardVariate();
    return y * y;
}

double NoncentralChiSquared::Variate(double degree, double noncentrality)
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

double NoncentralChiSquared::Variate() const
{
    if (k == 1)
        return variateForDegreeEqualOne();
    if (k > 1)
        return 2 * GammaRand::StandardVariate(halfK - 0.5) + variateForDegreeEqualOne();
    return 2 * GammaRand::StandardVariate(halfK + Y.Variate());
}

void NoncentralChiSquared::Sample(std::vector<double> &outputData) const
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

double NoncentralChiSquared::Mean() const
{
    return k + lambda;
}

double NoncentralChiSquared::Variance() const
{
    return 2 * (k + 2 * lambda);
}

std::complex<double> NoncentralChiSquared::CFImpl(double t) const
{
    std::complex<double> aux(1, -2 * t);
    std::complex<double> y(0, lambda * t);
    y /= aux;
    y -= halfK * std::log(aux);
    return std::exp(y);
}

double NoncentralChiSquared::Skewness() const
{
    double y = k + 2 * lambda;
    y = 2.0 / y;
    double z = y * std::sqrt(y);
    return z * (k + 3 * lambda);
}

double NoncentralChiSquared::ExcessKurtosis() const
{
    double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}
