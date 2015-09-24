#include "ParetoRand.h"
#include "UniformRand.h"

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
    alpha = shape;
    if (alpha <= 0)
        alpha = 1.0;
    alphaInv = 1.0 / alpha;
    
    xm = scale;
    if (xm <= 0)
        xm = 1.0;
    pdfCoef = alpha * std::pow(xm, alpha);
}

void ParetoRand::setShape(double shape)
{
    alpha = shape;
    if (alpha <= 0)
        alpha = 1.0;
    alphaInv = 1.0 / alpha;
    pdfCoef = alpha * std::pow(xm, alpha);
}

void ParetoRand::setScale(double scale)
{
    xm = scale;
    if (xm <= 0)
        xm = 1.0;
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
    // TODO: decide which is faster
    //return std::pow(UniformRand::standardVariate(), -1.0 / shape);
    return std::exp(ExponentialRand::variate(shape));
}

double ParetoRand::standardVariate(double shape)
{
    if (shape == 1)
        return variateForAlphaOne();
    if (shape == 2)
        return variateForAlphaTwo();
    return variateForCommonAlpha(shape);
}

double ParetoRand::variate(double shape, double scale)
{
    return scale * standardVariate(shape);
}

double ParetoRand::variate() const
{
    return xm * standardVariate(alpha);
}

void ParetoRand::sample(QVector<double> &outputData) const
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

double ParetoRand::Quantile(double p) const
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
