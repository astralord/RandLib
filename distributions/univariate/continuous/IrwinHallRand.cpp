#include "IrwinHallRand.h"

IrwinHallRand::IrwinHallRand(int number)
{
    SetNumber(number);
}

std::string IrwinHallRand::Name() const
{
    return "Irwin-Hall(" + toStringWithPrecision(GetNumber()) + ")";
}

void IrwinHallRand::SetNumber(int number)
{
    n = std::max(1, number);
    cdfCoef = 0.5 / RandMath::factorial(n);
    pdfCoef = n * cdfCoef;
}

double IrwinHallRand::f(double x) const
{
    if (x < 0 || x > n)
        return 0.0;

    double sum = 0.0;
    for (int i = 0; i <= n; ++i)
    {
        double y = x - i;
        double add = RandMath::binomialCoef(n, i);
        add *= std::pow(y, n - 1);
        add *= RandMath::sign(y);
        sum += (i % 2) ? -add : add;
    }
    return pdfCoef * sum;
}

double IrwinHallRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    if (x >= n)
        return 1.0;

    double sum = 0.0;
    for (int i = 0; i <= n; ++i)
    {
        double y = x - i;
        double add = RandMath::binomialCoef(n, i);
        add *= std::pow(y, n);
        add *= RandMath::sign(y);
        sum += (i % 2) ? -add : add;
    }

    return 0.5 + cdfCoef * sum;
}

double IrwinHallRand::Variate() const
{
    double sum = 0.0;
    for (int i = 0; i != n; ++i)
        sum += U.Variate();
    return sum;
}

double IrwinHallRand::Mean() const
{
    return 0.5 * n;
}

double IrwinHallRand::Variance() const
{
    return n / 12;
}

std::complex<double> IrwinHallRand::CF(double t) const
{
    return std::pow(U.CF(t), n);
}

double IrwinHallRand::quantileImpl(double p) const
{
    double root = p * n;
    if (RandMath::findRoot([this, p] (double x)
    {
        return F(x) - p;
    },
    0, n, root))
        return root;
    return NAN;
}

double IrwinHallRand::quantileImpl1m(double p) const
{
    double root = n - p * n;
    if (RandMath::findRoot([this, p] (double x)
    {
        double y = F(x) - 1;
        return y + p;
    },
    0, n, root))
        return root;
    return NAN;
}

double IrwinHallRand::Median() const
{
    return 0.5 * n;
}

double IrwinHallRand::Mode() const
{
    return 0.5 * n;
}

double IrwinHallRand::Skewness() const
{
    return 0.0;
}

double IrwinHallRand::ExcessKurtosis() const
{
    return -1.2 / n;
}

