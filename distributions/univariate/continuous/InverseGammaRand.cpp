#include "InverseGammaRand.h"

InverseGammaRand::InverseGammaRand(double shape, double scale) : GammaRand()
{
    setParameters(shape, scale);
}

std::string InverseGammaRand::name() const
{
    return "Inverse-Gamma(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

double InverseGammaRand::f(double x) const
{
    if (x <= 0)
        return 0.0;
    double y = -(alpha + 1) * std::log(x);
    y -= beta / x;
    y += pdfCoef;
    return std::exp(y);
}

double InverseGammaRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    double y = RandMath::logUpperIncGamma(alpha, beta / x);
    return std::exp(mLgammaShape + y);
}

double InverseGammaRand::variate() const
{
    return 1.0 / GammaRand::variate();
}

void InverseGammaRand::sample(std::vector<double> &outputData) const
{
    GammaRand::sample(outputData);
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
