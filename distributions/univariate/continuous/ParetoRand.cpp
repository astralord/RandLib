#include "ParetoRand.h"
#include "UniformRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    SetShape(shape);
    SetScale(scale);
}

String ParetoRand::Name() const
{
    return "Pareto(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void ParetoRand::SetShape(double shape)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Pareto distribution: shape should be positive");
    alpha = shape;
    logAlpha = std::log(alpha);
}

void ParetoRand::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Pareto distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

double ParetoRand::f(const double & x) const
{
    return (x < sigma) ? 0.0 : std::exp(logf(x));
}

double ParetoRand::logf(const double & x) const
{
    if (x < sigma)
        return -INFINITY;
    double logX = std::log(x);
    double y = logSigma - logX;
    y *= alpha;
    y -= logX;
    y += logAlpha;
    return y;
}

double ParetoRand::F(const double & x) const
{
    return (x > sigma) ? -std::expm1(alpha * std::log(sigma / x)) : 0.0;
}

double ParetoRand::S(const double & x) const
{
    return (x > sigma) ? std::pow(sigma / x, alpha) : 1.0;
}

double ParetoRand::variateForAlphaOne()
{
    return 1.0 / UniformRand::StandardVariate();
}

double ParetoRand::variateForAlphaTwo()
{
    return 1.0 / std::sqrt(UniformRand::StandardVariate());
}

double ParetoRand::variateForGeneralAlpha(double shape)
{
    return std::exp(ExponentialRand::StandardVariate() / shape);
}

double ParetoRand::StandardVariate(double shape)
{
    if (RandMath::areClose(shape, 1.0))
        return variateForAlphaOne();
    if (RandMath::areClose(shape, 2.0))
        return variateForAlphaTwo();
    return variateForGeneralAlpha(shape);
}

double ParetoRand::Variate(double shape, double scale)
{
    return (shape <= 0.0 || scale <= 0.0) ? NAN : scale * StandardVariate(shape);
}

double ParetoRand::Variate() const
{
    return sigma * StandardVariate(alpha);
}

void ParetoRand::Sample(std::vector<double> &outputData) const
{
    if (RandMath::areClose(alpha, 1.0)) {
        for (double &var : outputData)
            var = sigma * variateForAlphaOne();
    }
    else if (RandMath::areClose(alpha, 2.0)) {
        for (double &var : outputData)
            var = sigma * variateForAlphaTwo();
    }
    else {
        for (double &var : outputData)
            var = sigma * variateForGeneralAlpha(alpha);
    }
}

double ParetoRand::Mean() const
{
    return (alpha > 1) ? alpha * sigma / (alpha - 1) : INFINITY;
}

double ParetoRand::Variance() const
{
    if (alpha > 2)
    {
        double var = sigma / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

double ParetoRand::Median() const
{
    double y = alpha * logSigma + M_LN2;
    return std::exp(y / alpha);
}

double ParetoRand::Mode() const
{
    return sigma;
}

double ParetoRand::Skewness() const
{
    if (alpha <= 3)
        return INFINITY;
    double skewness = (alpha - 2.0) / alpha;
    skewness = std::sqrt(skewness);
    skewness *= (1 + alpha) / (alpha - 3);
    return skewness + skewness;
}

double ParetoRand::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    double numerator = alpha + 1;
    numerator *= alpha;
    numerator -= 6;
    numerator *= alpha;
    numerator -= 2;
    double denominator = alpha * (alpha - 3) * (alpha - 4);
    return 6.0 * numerator / denominator;
}

double ParetoRand::quantileImpl(double p) const
{
    double y = logSigma - std::log1p(-p) / alpha;
    return std::exp(y);
}

double ParetoRand::quantileImpl1m(double p) const
{
    double y = logSigma - std::log(p) / alpha;
    return std::exp(y);
}

double ParetoRand::Entropy() const
{
    return logSigma - logAlpha + 1.0 / alpha + 1;
}

void ParetoRand::Fit(const std::vector<double> &sample)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, "All elements in the sample should be positive"));

    /// Calculate alpha
    double logAverage = 0.0;
    for (double var : sample)
        logAverage += std::log(var / minVar);
    logAverage = sample.size() / logAverage;

    SetShape(logAverage);
    SetScale(minVar);
}
