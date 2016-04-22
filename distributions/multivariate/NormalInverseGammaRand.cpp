#include "NormalInverseGammaRand.h"
#include "../univariate/continuous/StudentTRand.h"
#include "../univariate/continuous/NormalRand.h"

NormalInverseGammaRand::NormalInverseGammaRand(double location, double precision, double shape, double rate)
{
    setParameters(location, precision, shape, rate);
}

std::string NormalInverseGammaRand::name()
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

    pdfCoef = std::sqrt(0.5 * lambda / M_PI);
    double alpha = Y.getShape();
    double beta = Y.getRate();
    cdfCoef = 1.0 / std::tgamma(alpha);
    pdfCoef *= std::pow(beta, alpha) * cdfCoef;
}

double NormalInverseGammaRand::f(DoublePair point) const
{
    double x = point.first, sigmaSq = point.second;
    if (sigmaSq <= 0)
        return 0.0;
    double sigma = std::sqrt(sigmaSq);
    double alpha = Y.getShape();
    double beta = Y.getRate();
    double y = std::pow(1.0 / sigmaSq, alpha + 1);
    y /= sigma;
    double degree = (x - mu);
    degree *= degree;
    degree *= lambda;
    degree += beta + beta;
    degree /= (sigmaSq + sigmaSq);
    y *= std::exp(-degree);
    return pdfCoef * y;
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
    double alpha = Y.getShape();
    double beta = Y.getRate();
    double z = beta /sigmaSq;
    y *= std::pow(z, alpha);
    y *= std::exp(-z);
    y /= (sigmaSq + sigmaSq);
    return cdfCoef * y;
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
    double alpha = Y.getShape();
    double beta = Y.getRate();
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

void NormalInverseGammaRand::getFirstMarginalDistribution(UnivariateProbabilityDistribution &distribution) const
{
    double alpha = Y.getShape(), beta = Y.getRate();
    StudentTRand X(2 * Y.getShape(), mu, std::sqrt(alpha * lambda / beta));
    distribution = X;
}

void NormalInverseGammaRand::getSecondMarginalDistribution(UnivariateProbabilityDistribution &distribution) const
{
    distribution = Y;
}

