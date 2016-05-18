#include "NoncentralChiSquared.h"

NoncentralChiSquared::NoncentralChiSquared(int degree, double noncentrality)
{
    setParameters(degree, noncentrality);
}

std::string NoncentralChiSquared::name()
{
    return "Noncentral Chi-Squared(" + toStringWithPrecision(getDegree()) + ", "
            + toStringWithPrecision(getNoncentrality()) + ")";
}

void NoncentralChiSquared::setParameters(int degree, double noncentrality)
{
    k = std::max(degree, 1);
    X.setDegree(k - 1);

    lambda = std::max(noncentrality, 0.0);
    sqrtLambda = std::sqrt(lambda);
}

double NoncentralChiSquared::f(double x) const
{
    // TODO:
    return x + NAN;
}

double NoncentralChiSquared::F(double x) const
{
    // TODO:
    return x + NAN;
}

double NoncentralChiSquared::variate() const
{
    double y = sqrtLambda + NormalRand::standardVariate();
    return y * y + X.variate();
}

void NoncentralChiSquared::sample(std::vector<double> &outputData) const
{
    X.sample(outputData);
    for (double & var : outputData) {
        double y = sqrtLambda + NormalRand::standardVariate();
        var += y * y;
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
