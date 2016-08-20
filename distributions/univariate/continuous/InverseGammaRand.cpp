#include "InverseGammaRand.h"

InverseGammaRand::InverseGammaRand(double shape, double rate)
{
    setParameters(shape, rate);
}

std::string InverseGammaRand::name() const
{
    return "Inverse-Gamma(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

void InverseGammaRand::setParameters(double shape, double rate)
{
    X.setParameters(shape, rate);
    alpha = X.getShape();
    beta = X.getRate();
}

double InverseGammaRand::f(double x) const
{
    if (x <= 0)
        return 0.0;
    return X.f(1.0 / x) / (x * x);
}

double InverseGammaRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    double y = RandMath::logUpperIncGamma(alpha, beta / x);
    return std::exp(y - getLogGammaFunction());
}

double InverseGammaRand::variate() const
{
    return 1.0 / X.variate();
}

void InverseGammaRand::sample(std::vector<double> &outputData) const
{
    X.sample(outputData);
    for (double &var : outputData)
        var = 1.0 / var;
}

double InverseGammaRand::Mean() const
{
    return (alpha > 1) ? beta / (alpha - 1) : INFINITY;
}

double InverseGammaRand::Variance() const
{
    if (alpha <= 2)
        return INFINITY;
    double var = beta / (alpha - 1);
    var *= var;
    return var / (alpha - 2);
}

double InverseGammaRand::QuantileImpl(double p) const
{
    return 1.0 / X.QuantileImpl(1 - p);
}

double InverseGammaRand::Mode() const
{
    return beta / (alpha + 1);
}

double InverseGammaRand::Skewness() const
{
    return (alpha > 3) ? 4 * std::sqrt(alpha - 2) / (alpha - 3) : INFINITY;
}

double InverseGammaRand::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    double numerator = 30 * alpha - 66;
    double denominator = (alpha - 3) * (alpha - 4);
    return numerator / denominator;
}
