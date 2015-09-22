#include "ArcsineRand.h"
#include "../BasicRandGenerator.h"

ArcsineRand::ArcsineRand(double minValue, double maxValue, double shape)
{
    setSupport(minValue, maxValue);
    setShape(shape);
}

std::string ArcsineRand::name()
{
    return "Arcsine(" + toStringWithPrecision(getMin()) + ", "
                      + toStringWithPrecision(getMax()) + ", "
                      + toStringWithPrecision(getShape()) + ")";
}

void ArcsineRand::setSupport(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    if (a > b)
        SWAP(a, b);

    if (RandMath::areEqual(b, a))
        b = a + 1.0;

    bma = b - a;
}

void ArcsineRand::setShape(double shape)
{
    BetaRand::setParameters(1.0 - shape, shape);
    pdfCoef = std::sin(M_PI * beta) * M_1_PI;
}

double ArcsineRand::f(double x) const
{
    return BetaRand::f((x - a) / bma) / bma;
}

double ArcsineRand::F(double x) const
{
    if (x <= a)
        return 0;
    if (x >= b)
        return 1;
    if (RandMath::areEqual(beta, 0.5))
    {
        double y = (x - a) / (b - a);
        y = std::sqrt(y);
        return M_2_PI * std::asin(y);
    }
    return BetaRand::F((x - a) / bma);
}

double ArcsineRand::variate() const
{
    return a + bma * BetaRand::variate();
}

void ArcsineRand::sample(QVector<double> &outputData) const
{
    BetaRand::sample(outputData);
    for (double & var : outputData)
        var = a + bma * var;
}

double ArcsineRand::Mean() const
{
    return a + bma * BetaRand::Mean();
}

double ArcsineRand::Variance() const
{
    return bma * bma* BetaRand::Variance();
}

double ArcsineRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (RandMath::areEqual(beta, 0.5))
    {
        double x = std::sin(0.5 * M_PI * p);
        return a + bma * x * x;
    }
    return a + bma * BetaRand::Quantile(p);
}

double ArcsineRand::Median() const
{
    return a + bma * BetaRand::Median();
}

double ArcsineRand::Mode() const
{
    /// x \in {a, b}
    return (signed)RandGenerator::variate() < 0 ? a : b;
}

