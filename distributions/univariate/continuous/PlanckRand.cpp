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

    pdfCoef = Z.getInverseZetaFunction();
    pdfCoef *= std::pow(b, a + 1);
    pdfCoef /= std::tgamma(a + 1);
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
    //TODO:
    return NAN;
}

double PlanckRand::Variance() const
{
    //TODO:
    return NAN;
}

double PlanckRand::Mode() const
{
    return a / b;
}
