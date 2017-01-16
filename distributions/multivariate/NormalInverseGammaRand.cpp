#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/StudentTRand.h"
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
    lambda = precision;
    if (lambda <= 0)
        lambda = 1.0;
    Y.SetParameters(shape, rate);
    alpha = Y.GetShape();
    beta = Y.GetRate();

    pdfCoef = 0.5 * std::log(0.5 * lambda / M_PI);
    pdfCoef += alpha * std::log(beta) - Y.GetLogGammaFunction();
}

double NormalInverseGammaRand::f(DoublePair point) const
{
    double x = point.first, sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double sigma = std::sqrt(sigmaSq);
    double y = (alpha + 1) - std::log(sigmaSq);
    double degree = (x - mu);
    degree *= degree;
    degree *= lambda;
    degree += beta + beta;
    degree /= (sigmaSq + sigmaSq);
    y -= degree;
    y = std::exp(pdfCoef + y);
    y /= sigma;
    return y;
}

double NormalInverseGammaRand::F(DoublePair point) const
{
    double x = point.first, sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double sigma = std::sqrt(sigmaSq);
    double y = 0.5 * lambda;
    y *= (x - mu) / sigma;
    y = std::erfc(-y);
    double z = beta /sigmaSq;
    double temp = alpha * std::log(z) - z;
    y *= std::exp(temp - Y.GetLogGammaFunction());
    y /= (sigmaSq + sigmaSq);
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

DoubleTriplet NormalInverseGammaRand::Covariance() const
{
    double var1;
    if (alpha <= 0.5)
        var1 = NAN;
    else if (alpha <= 1)
        var1 = INFINITY;
    else
        var1 = alpha * alpha * lambda / (beta * (alpha - 1));
    double var2 = Y.Variance();
    return std::make_tuple(var1, 0, var2);
}

double NormalInverseGammaRand::Correlation() const
{
    return 0.0;
}

void NormalInverseGammaRand::GetMarginalDistributions(UnivariateProbabilityDistribution<double> &distribution1, UnivariateProbabilityDistribution<double> &distribution2) const
{
    StudentTRand X(2 * Y.GetShape(), mu, std::sqrt(alpha * lambda / beta));
    distribution1 = X;
    distribution2 = Y;
}
