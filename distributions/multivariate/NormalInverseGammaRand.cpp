#include "NormalInverseGammaRand.h"
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

double NormalInverseGammaRand::f(double2d point) const
{
    double x = point.x, sigmaSq = point.y;
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

double NormalInverseGammaRand::F(double2d point) const
{
    double x = point.x, sigmaSq = point.y;
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

double2d NormalInverseGammaRand::variate() const
{
    double2d var;
    var.y = Y.variate();
    double coef = std::sqrt(var.y) / lambda;
    var.x = mu + coef * NormalRand::standardVariate();
    return var;
}

double2d NormalInverseGammaRand::Mean() const
{
    double2d mean;
    mean.x = mu;
    mean.y = Y.Mean();
    return mean;
}

bool NormalInverseGammaRand::Covariance(Matrix & matrix) const
{
    int n = matrix.height();
    int m = matrix.width();
    if (n != m || n != 2)
        return false;
    double alpha = Y.getShape();
    double beta = Y.getRate();
    if (alpha <= 0.5)
        matrix(1, 1) = NAN;
    else if (alpha <= 1)
        matrix(1, 1) = INFINITY;
    else
        matrix(1, 1) = alpha * alpha * lambda / (beta * (alpha - 1));
    matrix(1, 2) = matrix(2, 1) = 0;
    matrix(2, 2) = Y.Variance();
    return true;
}

