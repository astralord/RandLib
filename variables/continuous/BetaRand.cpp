#include "BetaRand.h"

BetaRand::BetaRand(double shape1, double shape2)
{
    setParameters(shape1, shape2);
}

std::string BetaRand::name()
{
    return "Beta(" + toStringWithPrecision(getAlpha()) + ", " + toStringWithPrecision(getBeta()) + ")";
}

void BetaRand::setParameters(double shape1, double shape2)
{
    double alpha = std::max(shape1, MIN_POSITIVE);
    X.setParameters(alpha, 1);

    double beta = std::max(shape2, MIN_POSITIVE);
    Y.setParameters(beta, 1);

    pdfCoef = std::tgamma(alpha + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
}

void BetaRand::setAlpha(double shape1)
{
    double alpha = std::max(shape1, MIN_POSITIVE);
    X.setParameters(alpha, 1);
    pdfCoef = std::tgamma(alpha + Y.getShape()) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
}

void BetaRand::setBeta(double shape2)
{
    double beta = std::max(shape2, MIN_POSITIVE);
    Y.setParameters(beta, 1);
    pdfCoef = std::tgamma(X.getShape() + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
}

double BetaRand::f(double x) const
{
    if (x < 0 || x > 1)
        return 0;
    double rv = std::pow(x, X.getShape() - 1);
    rv *= std::pow(1 - x, Y.getShape() - 1);
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

double BetaRand::variate() const
{
    if (X.getShape() == Y.getShape())
        return variateForEqualParameters();
    return variateForDifferentParameters();
}

void BetaRand::sample(QVector<double> &outputData)
{
    if (X.getShape() == Y.getShape()) {
        for (double &var : outputData)
            var = variateForEqualParameters();
    }
    else {
        for (double &var : outputData)
            var = variateForDifferentParameters();
    }
}

double BetaRand::variateForEqualParameters() const
{
    double alpha = X.getShape();
    int iter = 0;
    do {
        double u1 = UniformRand::standardVariate();
        double u2 = UniformRand::standardVariate();
        if (u2 <= std::pow(4 * u1 * (1 - u1), alpha - 1))
            return u1;
    } while (++iter <= 1e9); /// one billion should be enough
    return 0; /// fail
}

double BetaRand::variateForDifferentParameters() const
{
    double x = X.variate();
    return x / (x + Y.variate());
}
