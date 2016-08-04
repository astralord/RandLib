#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"

NormalInverseGammaRand::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    setParameters(location, precision, shape, rate);
}

std::string NormalInverseGammaRand::name() const
{
    return "Normal-Inverse-Gamma(" + toStringWithPrecision(getLocation()) + ", "
                                   + toStringWithPrecision(getPrecision()) + ", "
                                   + toStringWithPrecision(getShape()) + ", "
                                   + toStringWithPrecision(getRate()) + ")";
}

void NormalInverseGammaRand::setParameters(double location, double precision, double shape, double rate)
{
    mu = location;
    lambda = precision;
    if (lambda <= 0)
        lambda = 1.0;
    Y.setParameters(shape, rate);
    alpha = Y.getShape();
    beta = Y.getRate();

    cdfCoef = -Y.getLogGammaFunction();
    pdfCoef = 0.5 * std::log(0.5 * lambda / M_PI);
    pdfCoef += alpha * std::log(beta) + cdfCoef;
}

double NormalInverseGammaRand::f(DoublePair point) const
{
    double x = point.first, sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double sigma = std::sqrt(sigmaSq);
    double y = (alpha + 1) * std::log(1.0 / sigmaSq);
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
    y = std::erf(y);
    ++y;
    double z = beta /sigmaSq;
    double temp = alpha * std::log(z) - z;
    y *= std::exp(temp + cdfCoef);
    y /= (sigmaSq + sigmaSq);
    return y;
}

DoublePair NormalInverseGammaRand::variate() const
{
    DoublePair var;
    var.second = Y.variate();
    double coef = std::sqrt(var.second) / lambda;
    var.first = mu + coef * NormalRand::standardVariate();
    return var;
}

DoublePair NormalInverseGammaRand::Mean() const
{
    DoublePair mean;
    mean.first = mu;
    mean.second = Y.Mean();
    return mean;
}

void NormalInverseGammaRand::Covariance(SquareMatrix<2> &matrix) const
{
    if (alpha <= 0.5)
        matrix(0, 0) = NAN;
    else if (alpha <= 1)
        matrix(0, 0) = INFINITY;
    else
        matrix(0, 0) = alpha * alpha * lambda / (beta * (alpha - 1));
    matrix(0, 1) = matrix(1, 0) = 0;
    matrix(1, 1) = Y.Variance();
}

double NormalInverseGammaRand::Correlation() const
{
    return 0.0;
}

void NormalInverseGammaRand::getFirstMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    StudentTRand X(2 * Y.getShape(), mu, std::sqrt(alpha * lambda / beta));
    distribution = X;
}

void NormalInverseGammaRand::getSecondMarginalDistribution(UnivariateProbabilityDistribution<double> &distribution) const
{
    distribution = Y;
}

