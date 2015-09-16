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
    alpha = shape1;
    if (alpha <= 0)
        alpha = MIN_POSITIVE;
    X.setParameters(alpha, 1);

    beta = shape2;
    if (beta <= 0)
        beta = MIN_POSITIVE;
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
    alpha = shape1;
    if (alpha <= 0)
        alpha = MIN_POSITIVE;
    X.setParameters(alpha, 1);
    pdfCoef = std::tgamma(alpha + Y.getShape()) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    setVariateConstants();
}

void BetaRand::setBeta(double shape2)
{
    beta = shape2;
    if (beta <= 0)
        beta = MIN_POSITIVE;
    Y.setParameters(beta, 1);
    pdfCoef = std::tgamma(X.getShape() + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    setVariateConstants();
}

double BetaRand::f(double x) const
{
    if (x < 0 || x > 1)
        return 0;
    if (RandMath::areEqual(alpha, beta))
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
        double y = BetaRand::f(t);
        if (std::isnan(y) || std::isinf(y)) /// kind of hack
            return 0.0;
        return y;
    },
    0, x);
}

double BetaRand::variateArcsine() const
{
    double x = std::sin(UniformRand::variate(-M_PI, M_PI));
    return x * x;
}

double BetaRand::variateForSmallEqualParameters() const
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
    /// We need to storage variate coefficient only if alpha = beta and large enough
    if (alpha > edgeForGenerators && RandMath::areEqual(alpha, beta))
    {
        double t = 1.0 / (alpha + alpha + 1);
        variateCoef = M_E * std::sqrt(0.5 * M_PI * M_E * t);
        variateCoef *= std::pow(0.25 - 0.75 * t, alpha - 1);
        variateCoef *= pdfCoef; /// /= Beta(alpha, alpha)

        N.setMean(0.5);
        N.setVar(0.25 * t);
    }
}

double BetaRand::variate() const
{
    if (RandMath::areEqual(alpha, beta) && alpha == 0.5)
        return variateArcsine();
    if (!RandMath::areEqual(alpha, beta) || alpha < 1)
        return variateForDifferentParameters();
    if (alpha == 1)
        return UniformRand::standardVariate();
    if (alpha <= edgeForGenerators)
        return variateForSmallEqualParameters();
    return variateForLargeEqualParameters();
}

void BetaRand::sample(QVector<double> &outputData) const
{
    if (RandMath::areEqual(alpha, beta) && alpha == 0.5) {
        for (double &var : outputData)
            var = variateArcsine();
    }
    else if (!RandMath::areEqual(alpha, beta) || alpha < 1) {
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

double BetaRand::Mean() const
{
    return alpha / (alpha + beta);
}

double BetaRand::Variance() const
{
    double denominator = alpha + beta;
    denominator *= denominator * (denominator + 1);
    return alpha * beta / denominator;
}

double BetaRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    double root = p;
    if (RandMath::findRoot([this, p] (double x)
    {
        return BetaRand::F(x) - p;
    },
    0, 1, root))
        return root;
    return NAN;
}

double BetaRand::Median() const
{
    if (RandMath::areEqual(alpha, beta))
        return 0.5;
    return Quantile(0.5);
}

double BetaRand::Mode() const
{
    if (alpha > 1)
    {
        if (beta > 1)
            return (alpha - 1) / (alpha + beta - 2);
        return 1.0;
    }
    if (beta > 1)
        return 0.0;
    return (signed)RandGenerator::variate() < 0 ? 0 : 1;
}

double BetaRand::Skewness() const
{
    double skewness = (alpha + beta + 1) / (alpha * beta);
    skewness = std::sqrt(skewness);
    skewness *= (alpha - beta);
    skewness /= (alpha + beta + 2);
    return skewness + skewness;
}

double BetaRand::ExcessKurtosis() const
{
    double sum = alpha + beta;
    double kurtosis = alpha - beta;
    kurtosis *= kurtosis;
    kurtosis *= (sum + 1);
    kurtosis /= (alpha * beta * (sum + 2));
    --kurtosis;
    kurtosis /= (sum + 3);
    return 6 * kurtosis;
}


BaldingNicholsRand::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    setParameters(fixatingIndex, frequency);
}

std::string BaldingNicholsRand::name()
{
    return "Balding-Nichols(" + toStringWithPrecision(getFixatingIndex()) + ", " + toStringWithPrecision(getFrequency()) + ")";
}

void BaldingNicholsRand::setParameters(double fixatingIndex, double frequency)
{
    F = fixatingIndex;
    if (F <= 0)
        F = MIN_POSITIVE;
    if (F >= 1)
        F = 1.0 - MIN_POSITIVE;

    p = frequency;
    if (p <= 0)
        p = MIN_POSITIVE;
    if (p >= 1)
        p = 1.0 - MIN_POSITIVE;

    double frac = (1.0 - F) / F;
    BetaRand::setParameters(frac * p, frac * (1 - p));
}
