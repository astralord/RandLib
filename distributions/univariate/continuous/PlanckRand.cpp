#include "PlanckRand.h"

PlanckRand::PlanckRand(double shape, double scale)
{
    setParameters(shape, scale);
}

std::string PlanckRand::name()
{
    return "Planck(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void PlanckRand::setParameters(double shape, double scale)
{
    a = shape;
    if (a <= 0)
        a = 1.0;
    b = scale;
    if (b <= 0)
        b = 1.0;

    Z.setExponent(a + 1);
    G.setParameters(a + 1, 1.0);

    pdfCoef = std::log(Z.getInverseZetaFunction());
    pdfCoef += (a + 1) * std::log(b);
    pdfCoef -= std::lgamma(a + 1);
    pdfCoef = std::exp(pdfCoef);
}

double PlanckRand::f(double x) const
{
    if (x <= 0)
        return 0;
    return pdfCoef * std::pow(x, a) / std::expm1(b * x);
}

double PlanckRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    return RandMath::integral([this] (double t)
    {
        return f(t);
    },
    0, x);
}

double PlanckRand::variate() const
{
    double g = G.variate();
    double z = Z.variate();
    return g / (b * z);
}

double PlanckRand::Mean() const
{
    double mode = a / b, x = mode;
    do {
        x += mode;
    } while (f(x) > 1e-10);
    return RandMath::integral([this] (double t)
    {
        return t * f(t);
    },
    0, x);
}

double PlanckRand::Variance() const
{
    double mode = a / b, x = mode;
    do {
        x += mode;
    } while (f(x) > 1e-10);
    double mean = RandMath::integral([this] (double t)
    {
        return t * f(t);
    },
    0, x);
    double secondMoment = RandMath::integral([this] (double t)
    {
        return t * t * f(t);
    },
    0, x);
    return secondMoment - mean * mean;
}

double PlanckRand::Mode() const
{
    return a / b;
}
