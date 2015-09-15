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
    + toStringWithPrecision(getAlpha()) + ")";
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
    alpha = BetaRand::getBeta();
    pdfCoef = std::sin(M_PI * alpha) * M_1_PI;
}

double ArcsineRand::f(double x) const
{
    if (x < a || x > b)
        return 0;
    if (x == a || x == b)
        return INFINITY;
    if (RandMath::areEqual(shape, 0.5))
    {
        double y = (x - a) * (x - b);
        y = std::sqrt(y);
        return M_1_PI / y;
    }
    double y = std::pow(x, alpha);
    y *= std::pow(1 - x, 1.0 - alpha);
    return pdfCoef / y;
}

double ArcsineRand::F(double x) const
{
    if (x <= 0)
        return 0;
    if (x >= 1)
        return 1;
    return M_2_PI * std::asin(std::sqrt(x));
}

double ArcsineRand::Mean() const
{
    return 0.5;
}

double ArcsineRand::Variance() const
{

}

double ArcsineRand::Median() const
{
    return 0.5;
}

double ArcsineRand::Mode() const
{
    /// x \in {0, 1}
    return RandGenerator::variate() < 0 ? 0 : 1;
}

double ArcsineRand::Skewness() const
{
    return 0.0;
}

double ArcsineRand::ExcessKurtosis() const
{
    return -1.5;
}

