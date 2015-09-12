#include "ParetoRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    setParameters(shape, scale);
}

std::string ParetoRand::name()
{
    return "Pareto(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void ParetoRand::setParameters(double shape, double scale)
{
    alpha = std::max(shape, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    xm = std::max(scale, MIN_POSITIVE);
    pdfCoef = alpha * std::pow(xm, alpha);
}

void ParetoRand::setShape(double shape)
{
    alpha = std::max(shape, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    pdfCoef = alpha * std::pow(xm, alpha);
}

void ParetoRand::setScale(double scale)
{
    xm = std::max(scale, MIN_POSITIVE);
    pdfCoef = alpha * std::pow(xm, alpha);
}

double ParetoRand::f(double x) const
{
    return (x >= xm) ? pdfCoef / std::pow(x, alpha + 1) : 0;
}

double ParetoRand::F(double x) const
{
    return (x > xm) ? 1 - std::pow(xm / x, alpha) : 0;
}

double ParetoRand::variate() const
{
    return ParetoRand::variate(xm, alpha);
}

double ParetoRand::variate(double shape, double scale)
{
    if (shape == 1)
        return scale * variateForAlphaOne();
    if (shape == 2)
        return scale * variateForAlphaTwo();
    return scale * variateForCommonAlpha(shape);
}

double ParetoRand::variateForAlphaOne()
{
    return 1.0 / UniformRand::standardVariate();
}

double ParetoRand::variateForAlphaTwo()
{
    return 1.0 / std::sqrt(UniformRand::standardVariate());
}

double ParetoRand::variateForCommonAlpha(double shape)
{
    return std::exp(ExponentialRand::variate(shape));
}

void ParetoRand::sample(QVector<double> &outputData)
{
    if (alpha == 1) {
        for (double &var : outputData)
            var = xm * variateForAlphaOne();
    }
    else if (alpha == 2) {
        for (double &var : outputData)
            var = xm * variateForAlphaTwo();
    }
    else {
        for (double &var : outputData)
            var = xm * variateForCommonAlpha(alpha);
    }
}

double ParetoRand::E() const
{
    return (alpha > 1) ? alpha * xm / (alpha - 1) : INFINITY;
}

double ParetoRand::Var() const
{
    if (alpha > 2)
    {
        double var = xm / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

double ParetoRand::quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 1)
        return INFINITY;
    if (p == 0)
        return -INFINITY;
    return xm / std::pow(1 - p, alphaInv);
}

double ParetoRand::Median() const
{
    return xm * std::pow(2.0, alphaInv);
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

bool ParetoRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate xm
    double minVar = sample.at(0);
    for (double var : sample)
    {
        minVar = std::min(minVar, var);
        if (minVar <= 0)
            return false;
    }

    /// Calculate alpha
    double alpha = 0.0;
    for (double var : sample)
        alpha += std::log(var / minVar);
    alpha = sample.size() / alpha;

    setParameters(minVar, alpha);
    return true;
}
