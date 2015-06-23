#include "GammaRand.h"

GammaRand::GammaRand(double shape, double scale)
{
    setParameters(shape, scale);
}

void GammaRand::setParameters(double shape, double scale)
{
    k = std::max(shape, MIN_POSITIVE);
    theta = std::max(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;
    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
}

void GammaRand::setShape(double shape)
{
    k = std::max(shape, MIN_POSITIVE);
    cdfCoef = 1.0 / std::tgamma(k);
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
}

void GammaRand::setScale(double scale)
{
    theta = std::max(scale, MIN_POSITIVE);
    thetaInv = 1.0 / theta;
    pdfCoef = cdfCoef * std::pow(thetaInv, k);
}

double GammaRand::f(double x) const
{
    if (x < 0)
        return 0;
    double y = std::pow(x, k - 1);
    y *= std::exp(-x * thetaInv);
    return pdfCoef * y;
}

double GammaRand::F(double x) const
{
    return cdfCoef * RandMath::lowerIncGamma(k, x * thetaInv);
}

double GammaRand::value()
{
    return 0;
}
