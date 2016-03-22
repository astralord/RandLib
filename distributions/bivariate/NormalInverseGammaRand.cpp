#include "NormalInverseGammaRand.h"

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
    Y.setParameters(shape, 1.0 / rate);

    pdfCoef = std::sqrt(lambda / M_2_PI);
    double alpha = Y.getShape();
    double beta = Y.getRate();
    pdfCoef *= std::pow(beta, alpha) / std::tgamma(alpha);
}

double NormalInverseGammaRand::f(double2d grid) const
{
    double x = grid.x, sigmaSq = grid.y;
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

double NormalInverseGammaRand::F(double2d grid) const
{
    double sigmaSq = grid.y;
    if (sigmaSq <= 0)
        return 0.0;
    //TODO
    return NAN;
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

