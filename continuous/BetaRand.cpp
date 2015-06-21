#include "BetaRand.h"

BetaRand::BetaRand(double shape1, double shape2) :
    X(shape1, 1), Y(shape2, 1)
{
    setParameters(shape1, shape2);
}

void BetaRand::setParameters(double shape1, double shape2)
{
    alpha = qMax(shape1, MIN_POSITIVE);
    X.setShape(alpha);
    gammaA = std::tgamma(alpha);

    beta = qMax(shape2, MIN_POSITIVE);
    Y.setShape(beta);
    gammaB = std::tgamma(beta);

    pdfCoef = std::tgamma(alpha + beta) / (gammaA + gammaB);
}

void BetaRand::setAlpha(double shape1)
{
    alpha = qMax(shape1, MIN_POSITIVE);
    X.setShape(alpha);
    gammaA = std::tgamma(alpha);
    pdfCoef = std::tgamma(alpha + beta) / (gammaA + gammaB);
}

void BetaRand::setBeta(double shape2)
{
    beta = qMax(shape2, MIN_POSITIVE);
    Y.setShape(beta);
    gammaB = std::tgamma(beta);
    pdfCoef = std::tgamma(alpha + beta) / (gammaA + gammaB);
}

double BetaRand::pdf(double x)
{
    double rv = qPow(x, alpha - 1);
    rv *= qPow(1 - x, beta - 1);
    return pdfCoef * rv;
}

double BetaRand::cdf(double x)
{
    return x;
}

double BetaRand::value()
{
    double x = X.value();
    return x / (x + Y.value());
}
