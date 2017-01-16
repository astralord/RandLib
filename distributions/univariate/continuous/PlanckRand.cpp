#include "PlanckRand.h"

PlanckRand::PlanckRand(double shape, double scale)
{
    SetParameters(shape, scale);
}

std::string PlanckRand::Name() const
{
    return "Planck(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void PlanckRand::SetParameters(double shape, double scale)
{
    a = (shape > 0) ? shape : 1.0;
    b = (scale > 0) ? scale : 1.0;

    double ap1 = a + 1;
    Z.SetExponent(ap1);
    G.SetParameters(ap1, b);

    pdfCoef = std::log(Z.GetInverseZetaFunction());
    pdfCoef += ap1 * std::log(b);
    pdfCoef -= G.GetLogGammaFunction();
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
    double aux = pdfCoef + a * std::log(x);
    double integral1 = std::exp(aux) / (b * a);
    double integral2 = RandMath::integral([this] (double t)
    {
        if (t <= 0)
            return 0.0;
        double logT = std::log(t);
        double y = pdfCoef + a * logT;
        double expY = std::exp(y);
        double bt = b * t;
        double z = 1.0 / std::expm1(bt) - 1.0 / bt;
        return expY * z;
    },
    0, x);
    return integral1 + integral2;
}

double PlanckRand::Variate() const
{
    return G.Variate() / Z.Variate();
}

void PlanckRand::Sample(std::vector<double> &outputData) const
{
    G.Sample(outputData);
    for (double & var : outputData)
        var /= Z.Variate();
}

double PlanckRand::Mean() const
{
    double y = (a + 1) / b;
    y *= RandMath::zetaRiemann(a + 2);
    return Z.GetInverseZetaFunction() * y;
}

double PlanckRand::Variance() const
{
    double mean = Mean();
    double y = (a + 1) * (a + 2);
    y /= (b * b);
    y *= RandMath::zetaRiemann(a + 3);
    y *= Z.GetInverseZetaFunction();
    return y - mean * mean;
}

double PlanckRand::Mode() const
{
    if (a < 1)
        return 0.0;
    double y = -a * std::exp(-a);
    y = RandMath::W0Lambert(y);
    return (y + a) / b;
}
