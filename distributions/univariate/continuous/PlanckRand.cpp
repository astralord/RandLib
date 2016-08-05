#include "PlanckRand.h"

PlanckRand::PlanckRand(double shape, double scale)
{
    setParameters(shape, scale);
}

std::string PlanckRand::name() const
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
    pdfCoef -= G.getLogGammaFunction();
}

double PlanckRand::f(double x) const
{
    if (x < 0)
        return 0;
    if (x == 0)
    {
        if (a > 1)
            return 0.0;
        if (a == 1)
            return std::exp(pdfCoef) / b;
        return INFINITY;
    }
    double y = pdfCoef + a * std::log(x);
    return std::exp(y) / std::expm1(b * x);
}

double PlanckRand::F(double x) const
{
    if (x <= 0)
        return 0.0;

    if (a >= 1) {
        return RandMath::integral([this] (double t)
        {
            return f(t);
        },
        0, x);
    }

    /// split F(x) on two integrals
    double y = pdfCoef + a * std::log(x);
    double integral1 = std::exp(y) / (b * a);
    double integral2 = RandMath::integral([this] (double t)
    {
        if (t <= 0)
            return 0.0;
        double logT = std::log(t);
        double y = pdfCoef + a * logT;
        double z = y - logT;
        y = std::exp(y);
        y /= std::expm1(b * t);
        z = std::exp(z);
        z /= b;
        return y - z;
    },
    0, x);
    return integral1 + integral2;
}

double PlanckRand::variate() const
{
    double g = G.variate();
    double z = Z.variate();
    return g / (b * z);
}

double PlanckRand::Mean() const
{
    return RandMath::integral([this] (double t)
    {
        if (t <= 0 || t >= 1)
            return 0.0;
        double denom = 1.0 / (1.0 - t);
        double p = t * denom;
        double y = p * f(p);
        denom *= denom;
        return y * denom;
    },
    0, 1);
}

double PlanckRand::Variance() const
{
    double mean = Mean();
    double secondMoment = RandMath::integral([this] (double t)
    {
        if (t <= 0 || t >= 1)
            return 0.0;
        double denom = 1.0 / (1.0 - t);
        double p = t * denom;
        double y = p * p * f(p);
        denom *= denom;
        return y * denom;
    },
    0, 1);
    return secondMoment - mean * mean;
}

double PlanckRand::Mode() const
{
    return a / b;
}
