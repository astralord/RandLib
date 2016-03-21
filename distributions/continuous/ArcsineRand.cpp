#include "ArcsineRand.h"
#include "../discrete/BernoulliRand.h"

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

void ArcsineRand::setShape(double shape)
{
    BetaRand::setShapes(1.0 - shape, shape);
    cdfCoef = std::sin(M_PI * beta) * M_1_PI;
}

double ArcsineRand::F(double x) const
{
    if (beta == 0.5)
    {
        double y = (x - a) / bma;
        y = std::sqrt(y);
        return M_2_PI * std::asin(y);
    }
    return BetaRand::F(x);
}

double ArcsineRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (beta == 0.5)
    {
        double x = std::sin(0.5 * M_PI * p);
        return a + bma * x * x;
    }
    return BetaRand::Quantile(p);
}

double ArcsineRand::Mode() const
{
    /// x \in {a, b}
    return BernoulliRand::standardVariate() ? a : b;
}

