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
    setVariateConstants();
}

void BetaRand::setAlpha(double shape1)
{
    double alpha = std::max(shape1, MIN_POSITIVE);
    X.setParameters(alpha, 1);
    pdfCoef = std::tgamma(alpha + Y.getShape()) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    setVariateConstants();
}

void BetaRand::setBeta(double shape2)
{
    double beta = std::max(shape2, MIN_POSITIVE);
    Y.setParameters(beta, 1);
    pdfCoef = std::tgamma(X.getShape() + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    setVariateConstants();
}

double BetaRand::f(double x) const
{
    if (x < 0 || x > 1)
        return 0;
    double alpha = X.getShape(), beta = Y.getShape();
    if (std::fabs(alpha - beta) < MIN_POSITIVE)
        return pdfCoef * std::pow(x - x * x, alpha - 1);
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

double BetaRand::variate() const
{
    double alpha = X.getShape(), beta = Y.getShape();
    if (std::fabs(alpha - beta) > MIN_POSITIVE || alpha < 1)
        return variateForDifferentParameters();
    if (alpha == 1)
        return UniformRand::standardVariate();
    if (alpha <= edgeForGenerators)
        return variateForSmallEqualParameters();
    return variateForLargeEqualParameters();
}

void BetaRand::sample(QVector<double> &outputData)
{
    double alpha = X.getShape(), beta = Y.getShape();
    if (std::fabs(alpha - beta) > MIN_POSITIVE || alpha < 1) {
        for (double &var : outputData)
            var = variateForDifferentParameters();
    }
    else if (alpha == 1) {
        for (double &var : outputData)
            var = UniformRand::standardVariate();
    }
    else if (alpha <= edgeForGenerators) {
        for (double &var : outputData)
            var = variateForSmallEqualParameters();
    }
    else {
        for (double &var : outputData)
            var = variateForLargeEqualParameters();
    }
}

double BetaRand::variateForSmallEqualParameters() const
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

double BetaRand::variateForLargeEqualParameters() const
{
    int iter = 0;
    do {
        double n = N.variate();
        if (n < 0 || n > 1)
            continue;
        double u = UniformRand::standardVariate();
        if (u < N.f(n) / (variateCoef * f(n)))
            return n;
    } while (++iter <= 1e9); /// one billion should be enough
    return 0; /// fail
}

double BetaRand::variateForDifferentParameters() const
{
    double x = X.variate();
    return x / (x + Y.variate());
}

void BetaRand::setVariateConstants()
{
    double alpha = X.getShape();
    /// We need to storage variate coefficient only if alpha = beta and large enough
    if (alpha > edgeForGenerators && std::fabs(alpha - Y.getShape()) < MIN_POSITIVE)
    {
        double t = 1.0 / (alpha + alpha + 1);
        variateCoef = std::sqrt(0.5 * M_PI * M_E * M_E * M_E * t);
        variateCoef *= std::pow(0.25 * (1 - 3 * t), alpha - 1);
        variateCoef *= pdfCoef; /// /= Beta(alpha, alpha)

        N.setMean(0.5);
        N.setVar(0.25 * t);
    }
}
