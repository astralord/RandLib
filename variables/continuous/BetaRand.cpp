#include "BetaRand.h"

BetaRand::BetaRand(double shape1, double shape2)
{
    setParameters(shape1, shape2);
}

void BetaRand::setParameters(double shape1, double shape2)
{
    alpha = std::max(shape1, MIN_POSITIVE);
    X.setParameters(alpha, 1);
    gammaA = std::tgamma(alpha);

    beta = std::max(shape2, MIN_POSITIVE);
    Y.setParameters(beta, 1);
    gammaB = std::tgamma(beta);

    pdfCoef = std::tgamma(alpha + beta) / (gammaA * gammaB);
}

void BetaRand::setAlpha(double shape1)
{
    alpha = std::max(shape1, MIN_POSITIVE);
    X.setParameters(alpha, 1);
    gammaA = std::tgamma(alpha);
    pdfCoef = std::tgamma(alpha + beta) / (gammaA * gammaB);
}

void BetaRand::setBeta(double shape2)
{
    beta = std::max(shape2, MIN_POSITIVE);
    Y.setParameters(beta, 1);
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

double BetaRand::variate()
{
    if (alpha == beta)
    {
        int iter = 0;
        do {
            double u1 = UniformRand::standardVariate();
            double u2 = UniformRand::standardVariate();
            if (u2 <= std::pow(4 * u1 * (1 - u1), alpha - 1))
                return u1;
        } while (++iter <= 1e9); /// one billion should be enough
        return 0; /// fail
    }
    double x = X.variate();
    return x / (x + Y.variate());
}
