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

    lambda = std::max(noncentrality, 0.0);
    sqrtLambda = std::sqrt(lambda);

    if (k > 1) {
        X.setParameters(0.5 * (k - 1), 0.5);
    }
    else {
        X.setParameters(0.5 * k, 0.5);
        Y.setRate(0.5 * lambda);
    }

    if (k < 2) {
        cdfCoef = lambda + k * M_LN2;
        cdfCoef *= 0.5;
        cdfCoef += std::lgamma(0.5 * k);
        cdfCoef = std::exp(-cdfCoef);
    }
}

double NoncentralChiSquared::f(double x) const
{
    if (x < 0)
        return 0.0;
    double halfkm1 = 0.5 * k - 1;
    double y = RandMath::modifiedBesselFirstKind(std::sqrt(lambda * x), halfkm1);
    y *= std::pow(x / lambda, 0.5 * halfkm1);
    return 0.5 * y * std::exp(-0.5 * (x + lambda));
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

    double halfK = 0.5 * k;
    double halfKm1 = halfK - 1.0;
    double y = cdfCoef / halfK * std::pow(x, halfK);

    y += RandMath::integral([this, halfKm1] (double t)
    {
        if (t <= 0)
            return 0.0;
        return f(t) - cdfCoef * std::pow(t, halfKm1);
    }, 0, x);

    return y;
}

double NoncentralChiSquared::variateForDegreeEqualOne() const
{
    double y = sqrtLambda + NormalRand::standardVariate();
    return y * y;
}

double NoncentralChiSquared::variate() const
{
    double y = (k == 1) ? X.variate() : 0.0;
    return variateForDegreeEqualOne() + y;
}

void NoncentralChiSquared::sample(std::vector<double> &outputData) const
{
    if (k != 1)
        X.sample(outputData);
    else
        std::fill(outputData.begin(), outputData.end(), 0.0);
    if (k >= 1)
    {
        for (double & var : outputData)
            var += variateForDegreeEqualOne();
    }
    else
    {
        for (double & var : outputData)
            var += GammaRand::variate(Y.variate(), 0.5);
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
    return y / std::pow(aux, 0.5 * k);
}

double NoncentralChiSquared::Skewness() const
{
    double y = k + 2 * lambda;
    y = std::pow(2 / y, 1.5);
    return y * (k + 3 * lambda);
}

double NoncentralChiSquared::ExcessKurtosis() const
{
    double y = k + 2 * lambda;
    return 12 * (k + 4 * lambda) / (y * y);
}
