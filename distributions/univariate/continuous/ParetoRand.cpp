#include "ParetoRand.h"
#include "UniformRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    SetShape(shape);
    SetScale(scale);
}

std::string ParetoRand::Name() const
{
    return "Pareto(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void ParetoRand::SetShape(double shape)
{
    alpha = (shape > 0.0) ? shape : 1.0;
    logAlpha = std::log(alpha);
}

void ParetoRand::SetScale(double scale)
{
    xm = (scale > 0.0) ? scale : 1.0;
    logXm = std::log(xm);
}

double ParetoRand::f(double x) const
{
    return (x < xm) ? 0.0 : std::exp(logf(x));
}

double ParetoRand::logf(double x) const
{
    if (x < xm)
        return -INFINITY;
    double logX = std::log(x);
    double y = logXm - logX;
    y *= alpha;
    y -= logX;
    y += logAlpha;
    return y;
}

double ParetoRand::F(double x) const
{
    return (x > xm) ? -std::expm1(alpha * std::log(xm / x)) : 0.0;
}

double ParetoRand::S(double x) const
{
    return (x > xm) ? std::pow(xm / x, alpha) : 1.0;
}

double ParetoRand::variateForAlphaOne()
{
    return 1.0 / UniformRand::StandardVariate();
}

double ParetoRand::variateForAlphaTwo()
{
    return 1.0 / std::sqrt(UniformRand::StandardVariate());
}

double ParetoRand::variateForCommonAlpha(double shape)
{
    return std::exp(ExponentialRand::StandardVariate() / shape);
}

double ParetoRand::StandardVariate(double shape)
{
    if (RandMath::areClose(shape, 1.0))
        return variateForAlphaOne();
    if (RandMath::areClose(shape, 2.0))
        return variateForAlphaTwo();
    return variateForCommonAlpha(shape);
}

double ParetoRand::Variate(double shape, double scale)
{
    return scale * StandardVariate(shape);
}

double ParetoRand::Variate() const
{
    return xm * StandardVariate(alpha);
}

void ParetoRand::Sample(std::vector<double> &outputData) const
{
    if (RandMath::areClose(alpha, 1.0)) {
        for (double &var : outputData)
            var = xm * variateForAlphaOne();
    }
    else if (RandMath::areClose(alpha, 2.0)) {
        for (double &var : outputData)
            var = xm * variateForAlphaTwo();
    }
    else {
        for (double &var : outputData)
            var = xm * variateForCommonAlpha(alpha);
    }
}

double ParetoRand::Mean() const
{
    return (alpha > 1) ? alpha * xm / (alpha - 1) : INFINITY;
}

double ParetoRand::Variance() const
{
    if (alpha > 2)
    {
        double var = xm / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

double ParetoRand::Median() const
{
    double y = alpha * logXm + M_LN2;
    return std::exp(y / alpha);
}

double ParetoRand::Mode() const
{
    return xm;
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
    double y = alpha * logXm - std::log1p(-p);
    return std::exp(y / alpha);
}

double ParetoRand::quantileImpl1m(double p) const
{
    double y = alpha * logXm - std::log(p);
    return std::exp(y / alpha);
}

double ParetoRand::Entropy() const
{
    return logXm - logAlpha + 1.0 / alpha + 1;
}

bool ParetoRand::FitMLE(const std::vector<double> &sample)
{
    double n = sample.size();
    if (n <= 0)
        return false;

    /// Calculate xm
    double minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        return false;

    /// Calculate alpha
    double logAverage = 0.0;
    for (double var : sample)
        logAverage += std::log(var / minVar);
    logAverage = n / logAverage;

    SetShape(logAverage);
    SetScale(minVar);
    return true;
}
