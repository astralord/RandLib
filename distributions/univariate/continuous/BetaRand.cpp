#include "BetaRand.h"
#include "../discrete/BernoulliRand.h"
#include "UniformRand.h"

BetaRand::BetaRand(double shape1, double shape2, double minValue, double maxValue)
{
    setShapes(shape1, shape2);
    setSupport(minValue, maxValue);
}

std::string BetaRand::name()
{
    return "Beta(" + toStringWithPrecision(getAlpha()) + ", "
                   + toStringWithPrecision(getBeta()) + ", "
                   + toStringWithPrecision(getMin()) + ", "
                   + toStringWithPrecision(getMax()) + ")";
}

void BetaRand::setShapes(double shape1, double shape2)
{
    alpha = shape1;
    if (alpha <= 0)
        alpha = 1.0;
    X.setParameters(alpha, 1);

    beta = shape2;
    if (beta <= 0)
        beta = 1.0;
    Y.setParameters(beta, 1);

    if (alpha + beta > 30)
    {
        /// we use log(Gamma(x)) in order to avoid too big numbers
        double logGammaX = std::lgamma(alpha);
        double logGammaY = std::lgamma(beta);
        cdfCoef = std::lgamma(alpha + beta) - logGammaX - logGammaY;
        cdfCoef = std::exp(cdfCoef);
    }
    else {
        cdfCoef = std::tgamma(alpha + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    }
    pdfCoef = cdfCoef / bma;

    setVariateConstants();
}

void BetaRand::setSupport(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    if (a >= b)
        b = a + 1.0;

    bma = b - a;
    pdfCoef = cdfCoef / bma;
}

void BetaRand::setAlpha(double shape1)
{
    alpha = shape1;
    if (alpha <= 0)
        alpha = 1.0;
    X.setParameters(alpha, 1);
    cdfCoef = std::tgamma(alpha + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    pdfCoef = cdfCoef / bma;
    setVariateConstants();
}

void BetaRand::setBeta(double shape2)
{
    beta = shape2;
    if (beta <= 0)
        beta = 1.0;
    Y.setParameters(beta, 1);
    cdfCoef = std::tgamma(alpha + beta) * X.getInverseGammaFunction() * Y.getInverseGammaFunction();
    pdfCoef = cdfCoef / bma;
    setVariateConstants();
}

double BetaRand::f(double x) const
{
    if (x <= a || x >= b)
        return 0;

    /// Standardize
    x -= a;
    x /= bma;

    if (alpha == beta)
        return pdfCoef * std::pow(x - x * x, alpha - 1);
    double rv = std::pow(x, alpha - 1);
    rv *= std::pow(1 - x, beta - 1);
    return pdfCoef * rv;
}

double BetaRand::F(double x) const
{
    if (x <= a)
        return 0;
    if (x >= b)
        return 1;

    /// Standardize
    x -= a;
    x /= bma;

    if (alpha == beta && beta == 0.5)
        return M_2_PI * std::asin(std::sqrt(x));
    return cdfCoef * RandMath::incompleteBetaFun((x - a) / bma, alpha, beta);
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
    return NAN; /// fail
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
    return NAN; /// fail
}

double BetaRand::variateForDifferentParameters() const
{
    double x = X.variate();
    return x / (x + Y.variate());
}

void BetaRand::setVariateConstants()
{
    /// We need to storage variate coefficient only if alpha = beta and large enough
    if (alpha > edgeForGenerators && RandMath::areClose(alpha, beta))
    {
        double t = 1.0 / (alpha + alpha + 1);
        variateCoef = M_E * std::sqrt(0.5 * M_PI * M_E * t);
        variateCoef *= std::pow(0.25 - 0.75 * t, alpha - 1);
        variateCoef *= cdfCoef; /// /= Beta(alpha, alpha)

        N.setLocation(0.5);
        N.setVariance(0.25 * t);
    }
}

double BetaRand::variate() const
{
    double var;
    if (alpha == beta && alpha == 0.5)
        var = variateArcsine();
    else if (!(alpha == beta) || alpha < 1)
        var =  variateForDifferentParameters();
    else if (alpha == 1)
        var = UniformRand::standardVariate();
    else if (alpha <= edgeForGenerators)
        var = variateForSmallEqualParameters();
    else
        var = variateForLargeEqualParameters();
    return a + bma * var;
}

void BetaRand::sample(std::vector<double> &outputData) const
{
    if (alpha == beta && alpha == 0.5) {
        for (double &var : outputData)
            var = variateArcsine();
    }
    else if (!(alpha == beta) || alpha < 1) {
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

    /// Shift and scale
    for (double &var : outputData)
        var = a + bma * var;
}

double BetaRand::Mean() const
{
    double mean = alpha / (alpha + beta);
    return a + bma * mean;
}

double BetaRand::Variance() const
{
    double var = alpha + beta;
    var *= var * (var + 1);
    var = alpha * beta / var;
    return bma * bma * var;
}

double BetaRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    double root = p;
    if (alpha == beta && beta == 0.5)
    {
        double x = std::sin(0.5 * M_PI * p);
        return a + bma * x * x;
    }
    if (RandMath::findRoot([this, p] (double x)
    {
        return BetaRand::F(x) - p;
    },
    a, b, root))
        return root;
    return NAN;
}

double BetaRand::Median() const
{
    return (alpha == beta) ? 0.5 : Quantile(0.5);
}

double BetaRand::Mode() const
{
    double mode;
    if (alpha > 1)
        mode = (beta > 1) ? (alpha - 1) / (alpha + beta - 2) : 1.0;
    else
        mode = (beta > 1) ? 0.0 : BernoulliRand::standardVariate(); // WRONG!
    return a + bma * mode;
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


ArcsineRand::ArcsineRand(double shape, double minValue, double maxValue)
{
    setShape(shape);
    setSupport(minValue, maxValue);
}

std::string ArcsineRand::name()
{
    return "Arcsine(" + toStringWithPrecision(getMin()) + ", "
                      + toStringWithPrecision(getMax()) + ", "
                      + toStringWithPrecision(getShape()) + ")";
}

void ArcsineRand::setShape(double shape)
{
    BetaRand::setShapes(1.0 - shape, shape);
}


BaldingNicholsRand::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    setShapes(fixatingIndex, frequency);
}

std::string BaldingNicholsRand::name()
{
    return "Balding-Nichols(" + toStringWithPrecision(getFixatingIndex()) + ", " + toStringWithPrecision(getFrequency()) + ")";
}

void BaldingNicholsRand::setParameters(double fixatingIndex, double frequency)
{
    F = fixatingIndex;
    if (F <= 0 || F >= 1)
        F = 0.5;

    p = frequency;
    if (p <= 0 || p >= 1)
        p = 0.5;

    double frac = (1.0 - F) / F, fracP = frac * p;
    BetaRand::setShapes(fracP, frac - fracP);
}
