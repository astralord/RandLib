#include "ParetoRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    setParameters(shape, scale);
}

void ParetoRand::setName()
{
    nameStr = "Pareto(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void ParetoRand::setParameters(double shape, double scale)
{
    xm = std::max(shape, MIN_POSITIVE);
    alpha = std::max(scale, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    pdfCoef = alpha * std::pow(xm, alpha);
    setName();
}

void ParetoRand::setShape(double shape)
{
    xm = std::max(shape, MIN_POSITIVE);
    pdfCoef = alpha * std::pow(xm, alpha);
    setName();
}

void ParetoRand::setScale(double scale)
{
    alpha = std::max(scale, MIN_POSITIVE);
    alphaInv = 1.0 / alpha;
    pdfCoef = alpha * std::pow(xm, alpha);
    setName();
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
    return xm * std::exp(ExponentialRand::variate(alpha));
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
