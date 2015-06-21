#include "GammaRand.h"

GammaRand::GammaRand(double shape, double scale)
{
    setParameters(shape, scale);
}

void GammaRand::setParameters(double shape, double scale)
{
    k = qMax(shape, MIN_POSITIVE);
    theta = qMax(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;
    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * qPow(thetaInv, k);
}

void GammaRand::setShape(double shape)
{
    k = qMax(shape, MIN_POSITIVE);
    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * qPow(thetaInv, k);
}

void GammaRand::setScale(double scale)
{
    theta = qMax(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;
    pdfCoef = cdfCoef * qPow(thetaInv, k);
}

double GammaRand::pdf(double x)
{
    if (x < 0)
        return 0;
    double y = qPow(x, k - 1);
    y *= qExp(-x * thetaInv);
    return pdfCoef * y;
}

double GammaRand::cdf(double x)
{
    return cdfCoef * x;
}

double GammaRand::value()
{
    return 0;
}
