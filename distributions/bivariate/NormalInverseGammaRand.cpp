#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/NormalRand.h"

NormalInverseGammaRand::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    SetParameters(location, precision, shape, rate);
}

String NormalInverseGammaRand::Name() const
{
    return "Normal-Inverse-Gamma(" + toStringWithPrecision(GetLocation()) + ", "
                                   + toStringWithPrecision(GetPrecision()) + ", "
                                   + toStringWithPrecision(GetShape()) + ", "
                                   + toStringWithPrecision(GetRate()) + ")";
}

void NormalInverseGammaRand::SetParameters(double location, double precision, double shape, double rate)
{
    if (precision <= 0.0)
        throw std::invalid_argument("Precision of Normal-Inverse-Gamma distribution should be positive");
    if (shape <= 0.0)
        throw std::invalid_argument("Shape of Normal-Inverse-Gamma distribution should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Rate of Normal-Inverse-Gamma distribution should be positive");

    mu = location;
    lambda = precision;

    Y.SetParameters(shape, rate);
    alpha = Y.GetShape();
    beta = Y.GetRate();
    X.SetDegree(2 * alpha);
    X.SetLocation(mu);
    X.SetScale(std::sqrt(alpha * lambda / beta));

    pdfCoef = 0.5 * std::log(0.5 * lambda / M_PI);
    pdfCoef += alpha * Y.GetLogRate() - Y.GetLogGammaShape();
}

double NormalInverseGammaRand::f(const DoublePair &point) const
{
    return (point.second > 0.0) ? std::exp(logf(point)) : 0.0;
}

double NormalInverseGammaRand::logf(const DoublePair &point) const
{
    double sigmaSq = point.second;
    if (sigmaSq <= 0)
        return -INFINITY;
    double x = point.first;
    double y = alpha + 1 - 1.5 * std::log(sigmaSq);
    double degree = x - mu;
    degree *= degree;
    degree *= lambda;
    degree += 2 * beta;
    degree *= 0.5 / sigmaSq;
    y -= degree;
    return pdfCoef + y;
}

double NormalInverseGammaRand::F(const DoublePair &point) const
{
    double sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double x = point.first;
    double sigma = std::sqrt(sigmaSq);
    double y = 0.5 * lambda;
    y *= (mu - x) / sigma;
    y = std::erfc(y);
    double z = beta / sigmaSq;
    double temp = alpha * std::log(z) - z;
    y *= std::exp(temp - Y.GetLogGammaShape());
    y *= 0.5 / sigmaSq;
    return y;
}

DoublePair NormalInverseGammaRand::Variate() const
{
    DoublePair var;
    var.second = Y.Variate();
    double coef = std::sqrt(var.second) / lambda;
    var.first = mu + coef * NormalRand::StandardVariate();
    return var;
}

double NormalInverseGammaRand::Correlation() const
{
    return 0.0;
}
