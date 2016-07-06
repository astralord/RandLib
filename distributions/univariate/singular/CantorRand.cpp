#include "CantorRand.h"
#include "../discrete/BernoulliRand.h"
#include "../continuous/UniformRand.h"

double CantorRand::table[CantorRand::n] = {0};
const bool CantorRand::dummy = CantorRand::setupTable();

bool CantorRand::setupTable()
{
    table[0] = 0.33333333333333333333;
    for (int i = 1; i != n; ++i)
        table[i] = table[i - 1] / 3.0;
    return true;
}

CantorRand::CantorRand()
{
}

std::string CantorRand::name() const
{
    return "Cantor";
}

double CantorRand::F(double x) const
{
    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    double a = 0.0, b = 1.0;
    double v = 0.5, d = 0.5, delta = 1.0;
    while (delta > MIN_POSITIVE)
    {
        double delta = (b - a) / 3.0;
        if (x < a + delta)
        {
            b = a + delta;
            d *= 0.5;
            v -= d;
        }
        else if (x > b - delta)
        {
            a = b - delta;
            d *= 0.5;
            v += d;
        }
        else
            return v;
    }
    return v;
}

double CantorRand::variate() const
{
    long double sum = 0.0;
    for (int i = 0; i != n; ++i) {
        sum += table[i] * BernoulliRand::standardVariate();
    }
    return sum + sum;
}

double CantorRand::Mean() const
{
    return 0.5;
}

double CantorRand::Variance() const
{
    return 0.125;
}

double CantorRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;

    double root = 0.5;
    if (RandMath::findRoot([this, p] (double x)
    {
        return F(x) - p;
    }, 0.0, 1.0, root))
        return root;
    return NAN;
}

std::complex<double> CantorRand::CF(double t) const
{
    if (t == 0)
        return 1;
    double prod = 1.0;
    for (int i = 0; i != n; ++i)
        prod *= std::cos(table[i]);
    std::complex<double> y(0.0, 0.5 * t);
    y = std::exp(y);
    return y * prod;
}

double CantorRand::Median() const
{
    double U = UniformRand::standardVariate();
    return (U + 1.0) / 3.0;
}

double CantorRand::Skewness() const
{
    return 0.0;
}

double CantorRand::ExcessKurtosis() const
{
    return -1.6;
}

