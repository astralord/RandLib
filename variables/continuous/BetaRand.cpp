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

    if (alpha + beta > 30)
    {
        /// we use log(Gamma(x)) in order to avoid too big numbers
        double logGammaX = std::log(X.getInverseGammaFunction());
        double logGammaY = std::log(Y.getInverseGammaFunction());
        pdfCoef = std::lgamma(alpha + beta) + logGammaX + logGammaY;
        pdfCoef = std::exp(pdfCoef);
    }
    else {
        pdfCoef = std::tgamma(alpha + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    }
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
    return RandMath::integral([this] (double t)
    {
        return f(t);
    },
    0, x);
}

double BetaRand::variate() const
{
    double alpha = X.getShape(), beta = Y.getShape();
    if (!RandMath::areEqual(alpha, beta) || alpha < 1)
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
    if (!RandMath::areEqual(alpha, beta) || alpha < 1) {
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
        variateCoef = M_E * std::sqrt(0.5 * M_PI * M_E * t);
        variateCoef *= std::pow(0.25 - 0.75 * t, alpha - 1);
        variateCoef *= pdfCoef; /// /= Beta(alpha, alpha)

        N.setMean(0.5);
        N.setVar(0.25 * t);
    }
}

double BetaRand::E() const
{
    return X.getShape() / (X.getShape() + Y.getShape());
}

double BetaRand::Var() const
{
    double alpha = X.getShape();
    double beta = Y.getShape();
    double denominator = alpha + beta;
    denominator *= denominator * (denominator + 1);
    return alpha * beta / denominator;
}

double BetaRand::quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    double root = 0;
    RandMath::findRoot([this, p] (double x)
    {
        return F(x) - p;
    },
    0, 1, root);
    return root;
}

double BetaRand::Mode() const
{
    double alpha = X.getShape();
    double beta = Y.getShape();
    return (alpha - 1) / (alpha + beta - 2);
}

double BetaRand::Skewness() const
{
    double alpha = X.getShape();
    double beta = Y.getShape();
    double skewness = (alpha + beta + 1) / (alpha * beta);
    skewness = std::sqrt(skewness);
    skewness *= (alpha - beta);
    skewness /= (alpha + beta + 2);
    return skewness + skewness;
}

double BetaRand::ExcessKurtosis() const
{
    double alpha = X.getShape();
    double beta = Y.getShape();
    double sum = alpha + beta;
    double kurtosis = alpha - beta;
    kurtosis *= kurtosis;
    kurtosis *= (sum + 1);
    kurtosis /= (alpha * beta * (sum + 2));
    --kurtosis;
    kurtosis /= (sum + 3);
    return 6 * kurtosis;
}
