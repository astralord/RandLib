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
    pdfCoef -= G.GetLogGammaShape();
}

double PlanckRand::leveledPdf(double t) const
{
    if (t <= 0)
        return 0.0;
    double y = pdfCoef + a * std::log(t);
    double expY = std::exp(y);
    double bt = b * t;
    double z = 1.0 / std::expm1(bt) - 1.0 / bt;
    return expY * z;
}

double PlanckRand::f(const double & x) const
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

double PlanckRand::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0)
    {
        if (a > 1)
            return -INFINITY;
        if (a == 1)
            return pdfCoef - G.GetLogRate();
        return INFINITY;
    }
    double y = pdfCoef + a * std::log(x);
    return y - RandMath::logexpm1(b * x);
}

double PlanckRand::F(const double & x) const
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

    /// split F(x) by two integrals
    double aux = pdfCoef + a * std::log(x);
    double integral1 = std::exp(aux) / (b * a);
    double integral2 = RandMath::integral([this] (double t)
    {
        return leveledPdf(t);
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

// TODO: implement also skewness and kurtosis, as analytical solution is known

double PlanckRand::Mode() const
{
    if (a <= 1)
        return 0.0;
    double y = -a * std::exp(-a);
    y = RandMath::W0Lambert(y);
    return (y + a) / b;
}

std::complex<double> PlanckRand::CFImpl(double t) const
{
    if (a >= 1)
        return ContinuousDistribution::CFImpl(t);

    /// We have singularity point at 0 for real part,
    /// so we split the integral in two intervals:
    /// First one from 0 to 1, for which we integrate
    /// numerically leveled pdf and add known solution for level.
    /// Second one from 1 to infinity, for which we use
    /// simple expected value for the rest of the function
    double re1 = RandMath::integral([this, t] (double x)
    {
        return std::cos(t * x) * leveledPdf(x);
    },
    0.0, 1.0);

    double re2 = ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    },
    1.0, INFINITY);

    double re3 = t * RandMath::integral([this, t] (double x)
    {
        if (x <= 0.0)
            return 0.0;
        return std::sin(t * x) * std::pow(x, a);
    },
    0.0, 1.0);

    re3 += std::cos(t);
    re3 *= std::exp(pdfCoef) / (b * a);

    double im = ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    },
    0.0, INFINITY);

    double re = re1 + re2 + re3;
    return std::complex<double>(re, im);
}
