#include "ArcsineRand.h"
#include "../BasicRandGenerator.h"

ArcsineRand::ArcsineRand(double minValue = 0, double maxValue = 1, double shape = 0.5)
{
    setSupport(minValue, maxValue);
    setShape(shape);
}

std::string ArcsineRand::name()
{
    return "Arcsine(" + toStringWithPrecision(getMin()) + ", "
                      + toStringWithPrecision(getMax()) + ", ";
                      + toStringWithPrecision(getShape()) + ")";
}

void ArcsineRand::setSupport(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    if (a > b)
        SWAP(a, b);

    if (RandMath::areEqual(b, a))
        b = a + MIN_POSITIVE;
}

void ArcsineRand::setShape(double shape)
{
    BetaRand::setParameters(1.0 - shape, shape);
    pdfCoef = std::sin(M_PI * beta) * M_1_PI;
}

double ArcsineRand::f(double x) const
{
    if (x < a || x > b)
        return 0;
    if (x == a || x == b)
        return INFINITY;
    if (RandMath::areEqual(beta, 0.5))
    {
        double y = (x - a) * (x - b);
        y = std::sqrt(y);
        return M_1_PI / y;
    }
    double y = std::pow(x, beta);
    y *= std::pow(1.0 - x, 1.0 - beta);
    return pdfCoef / y;
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
    return BetaRand::F((x - a) / (b - a));
}

double ArcsineRand::variate() const
{
    return a + (b - a) * BetaRand::variate();
}

double ArcsineRand::Mean() const
{
    return a + (b - a) * BetaRand::Mean();
}

double ArcsineRand::Variance() const
{
    return (b - a) * (b - a) * BetaRand::Variance();
}

double ArcsineRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (RandMath::areEqual(beta, 0.5))
    {
        double x = std::sin(0.5 * M_PI * p);
        return a + (b - a) * x * x;
    }
    return a + (b - a) * BetaRand::Quantile(p);
}

double ArcsineRand::Mode() const
{
    /// x \in {a, b}
    return RandGenerator::variate() < 0 ? a : b;
}

double ArcsineRand::Skewness() const
{
    return BetaRand::Skewness();
}

double ArcsineRand::ExcessKurtosis() const
{
    return BetaRand::ExcessKurtosis();
}

