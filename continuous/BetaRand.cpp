#include "BetaRand.h"

BetaRand::BetaRand(double shape1, double shape2) :
    X(shape1, 1), Y(shape2, 1)
{
    setParameters(shape1, shape2);
}

void BetaRand::setParameters(double shape1, double shape2)
{
    alpha = std::max(shape1, MIN_POSITIVE);
    X.setShape(alpha);
    gammaA = std::tgamma(alpha);

    beta = std::max(shape2, MIN_POSITIVE);
    Y.setShape(beta);
    gammaB = std::tgamma(beta);

    pdfCoef = std::tgamma(alpha + beta) / (gammaA * gammaB);
}

void BetaRand::setAlpha(double shape1)
{
    alpha = std::max(shape1, MIN_POSITIVE);
    X.setShape(alpha);
    gammaA = std::tgamma(alpha);
    pdfCoef = std::tgamma(alpha + beta) / (gammaA * gammaB);
}

void BetaRand::setBeta(double shape2)
{
    beta = std::max(shape2, MIN_POSITIVE);
    Y.setShape(beta);
    gammaB = std::tgamma(beta);
    pdfCoef = std::tgamma(alpha + beta) / (gammaA * gammaB);
}

double BetaRand::f(double x) const
{
    if (x < 0 || x > 1)
        return 0;
    double rv = std::pow(x, alpha - 1);
    rv *= std::pow(1 - x, beta - 1);
    return pdfCoef * rv;
}

double BetaRand::F(double x) const
{
    if (x <= 0)
        return 0;
    if (x >= 1)
        return 1;
    return x;
}

double BetaRand::value()
{
    double x = X.value();
    return x / (x + Y.value());
}
