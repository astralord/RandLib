#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/NormalRand.h"

NormalInverseGammaRand::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    SetParameters(location, precision, shape, rate);
}

std::string NormalInverseGammaRand::Name() const
{
    return "Normal-Inverse-Gamma(" + toStringWithPrecision(GetLocation()) + ", "
                                   + toStringWithPrecision(GetPrecision()) + ", "
                                   + toStringWithPrecision(GetShape()) + ", "
                                   + toStringWithPrecision(GetRate()) + ")";
}

void NormalInverseGammaRand::SetParameters(double location, double precision, double shape, double rate)
{
    mu = location;
    lambda = (precision > 0.0) ? precision : 1.0;

    //StudentTRand X(2 * Y.GetShape(), mu, std::sqrt(alpha * lambda / beta));
    X.SetDegree(2 * shape);
    X.SetLocation(mu);
    X.SetScale(std::sqrt(alpha * lambda / beta));
    Y.SetParameters(shape, rate);
    alpha = Y.GetShape();
    beta = Y.GetRate();

    pdfCoef = 0.5 * std::log(0.5 * lambda / M_PI);
    pdfCoef += alpha * std::log(beta) - Y.GetLogGammaFunction();
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
    y *= std::exp(temp - Y.GetLogGammaFunction());
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

DoublePair NormalInverseGammaRand::Mean() const
{
    DoublePair mean;
    mean.first = mu;
    mean.second = Y.Mean();
    return mean;
}

double NormalInverseGammaRand::Correlation() const
{
    return 0.0;
}
