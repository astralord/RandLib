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
    xm = std::max(scale, MIN_POSITIVE);
    alpha = std::max(shape, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
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

double ParetoRand::Median() const
{
    return xm * std::pow(2.0, alphaInv);
}

double ParetoRand::Mode() const
{
    return xm;
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
