#include "IrwinHallRand.h"

IrwinHallRand::IrwinHallRand(int number)
{
    setNumber(number);
}

std::string IrwinHallRand::name()
{
    return "Irwin-Hall(" + toStringWithPrecision(getNumber()) + ")";
}

void IrwinHallRand::setNumber(int number)
{
    n = std::max(1, number);
    cdfCoef = 0.5 / RandMath::factorial(n);
    pdfCoef = n * cdfCoef;
}

int IrwinHallRand::getNumber()
{
    return n;
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
        add *= SIGN(y);
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
        add *= SIGN(y);
        sum += (i % 2) ? -add : add;
    }

    return 0.5 + cdfCoef * sum;
}

double IrwinHallRand::variate() const
{
    double sum = 0.0;
    for (int i = 0; i != n; ++i)
        sum += U.variate();
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

