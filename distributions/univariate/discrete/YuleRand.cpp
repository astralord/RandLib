#include "YuleRand.h"

YuleRand::YuleRand(double shape) :
X(shape, 1.0)
{
    SetShape(shape);
}

std::string YuleRand::Name() const
{
    return "Yule(" + toStringWithPrecision(GetShape()) + ")";
}

void YuleRand::SetShape(double shape)
{
    ro = (shape <= 0.0) ? 1.0 : shape;
    lgamma1pRo = std::lgamma(ro + 1);
    X.SetShape(ro);
}

double YuleRand::P(int k) const
{
    if (k < 1)
        return 0;
    double y = lgamma1pRo;
    y += std::lgamma(k);
    y -= std::lgamma(k + ro + 1);
    y = std::exp(y);
    return ro * y;
}

double YuleRand::F(int k) const
{
    if (k < 1)
        return 0;
    double y = lgamma1pRo;
    y += std::lgamma(k);
    y -= std::lgamma(k + ro + 1);
    y = std::exp(y);
    return 1 - k * y;
}

int YuleRand::Variate() const
{
    double prob = 1.0 / X.Variate();
    return GeometricRand::Variate(prob) + 1;
}

int YuleRand::Variate(double shape)
{
    double prob = 1.0 / ParetoRand::StandardVariate(shape);
    return GeometricRand::Variate(prob) + 1;
}

double YuleRand::Mean() const
{
    return (ro <= 1) ? INFINITY : ro / (ro - 1);
}

double YuleRand::Variance() const
{
    if (ro <= 2)
        return INFINITY;
    double aux = ro / (ro - 1);
    return aux * aux / (ro - 2);
}

int YuleRand::Mode() const
{
    return 1;
}

double YuleRand::Skewness() const
{
    if (ro <= 3)
        return INFINITY;
    double skewness = ro + 1;
    skewness *= skewness;
    skewness *= std::sqrt(ro - 2);
    return skewness / (ro * (ro - 3));
}

double YuleRand::ExcessKurtosis() const
{
    if (ro <= 4)
        return INFINITY;
    double numerator = 11 * ro * ro - 49;
    numerator *= ro;
    numerator -= 22;
    double denominator = ro * (ro - 4) * (ro - 3);
    return ro + 3 + numerator / denominator;
}
