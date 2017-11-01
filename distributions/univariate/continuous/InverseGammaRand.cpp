#include "InverseGammaRand.h"

InverseGammaRand::InverseGammaRand(double shape, double rate)
{
    SetParameters(shape, rate);
}

String InverseGammaRand::Name() const
{
    return "Inverse-Gamma(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetRate()) + ")";
}

void InverseGammaRand::SetParameters(double shape, double rate)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Inverse-Gamma distribution: shape should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Inverse-Gamma distribution: rate should be positive");
    X.SetParameters(shape, rate);
    alpha = X.GetShape();
    beta = X.GetRate();
    pdfCoef = -X.GetLogGammaShape() + alpha * X.GetLogRate();
}

double InverseGammaRand::f(const double & x) const
{
    return (x > 0.0) ? std::exp(logf(x)) : 0.0;
}

double InverseGammaRand::logf(const double & x) const
{
    if (x <= 0.0)
        return -INFINITY;
    double logX = std::log(x);
    double y = -(alpha - 1.0) * logX;
    y -= beta / x;
    y += pdfCoef;
    return y - 2 * logX;
}

double InverseGammaRand::F(const double & x) const
{
    return (x > 0.0) ? X.S(1.0 / x) : 0.0;
}

double InverseGammaRand::S(const double & x) const
{
    return (x > 0.0) ? X.F(1.0 / x) : 1.0;
}

double InverseGammaRand::Variate() const
{
    return 1.0 / X.Variate();
}

void InverseGammaRand::Sample(std::vector<double> &outputData) const
{
    X.Sample(outputData);
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

double InverseGammaRand::quantileImpl(double p) const
{
    return 1.0 / X.Quantile1m(p);
}

double InverseGammaRand::quantileImpl1m(double p) const
{
    return 1.0 / X.Quantile(p);
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
    double numerator = 30 * alpha - 66.0;
    double denominator = (alpha - 3) * (alpha - 4);
    return numerator / denominator;
}
