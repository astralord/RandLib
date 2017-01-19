#include "ParetoRand.h"
#include "UniformRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    SetParameters(shape, scale);
}

std::string ParetoRand::Name() const
{
    return "Pareto(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void ParetoRand::SetParameters(double shape, double scale)
{
    alpha = (shape > 0.0) ? shape : 1.0;
    alphaInv = 1.0 / alpha;
    xm = (scale > 0.0) ? scale : 1.0;
    alphaLogXm = alpha * std::log(xm);
}

void ParetoRand::SetShape(double shape)
{
    alpha = (shape > 0.0) ? shape : 1.0;
    alphaInv = 1.0 / alpha;
    alphaLogXm = alpha * std::log(xm);
}

void ParetoRand::SetScale(double scale)
{
    xm = (scale > 0.0) ? scale : 1.0;
    alphaLogXm = alpha * std::log(xm);
}

double ParetoRand::f(double x) const
{
    if (x < xm)
        return 0.0;
    double y = alphaLogXm - (alpha + 1) * std::log(x);
    return alpha * std::exp(y);
}

double ParetoRand::F(double x) const
{
    return (x > xm) ? 1 - std::pow(xm / x, alpha) : 0;
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
    return std::exp(ExponentialRand::Variate(shape));
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

double ParetoRand::quantileImpl(double p) const
{
    double y = alphaLogXm - std::log1p(-p);
    return std::exp(alphaInv * y);
}

double ParetoRand::quantileImpl1m(double p) const
{
    double y = alphaLogXm - std::log(p);
    return std::exp(alphaInv * y);
}

double ParetoRand::Median() const
{
    double y = alphaLogXm + M_LN2;
    return std::exp(alphaInv * y);
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

double ParetoRand::Entropy() const
{
    return std::log(xm * alphaInv) + alphaInv + 1;
}

bool ParetoRand::FitMLE(const std::vector<double> &sample)
{
    double n = sample.size();
    if (n <= 0)
        return false;

    /// Calculate xm
    double minVar = sample.at(0);
    if (minVar <= 0)
        return false;
    for (double var : sample)
    {
        if (minVar < var)
        {
            if (minVar <= 0)
                return false;
            minVar = var;
        }
    }

    /// Calculate alpha
    double logAverage = 0.0;
    for (double var : sample)
        logAverage += std::log(var / minVar);
    logAverage = n / logAverage;

    SetParameters(logAverage, minVar);
    return true;
}
