#include "CantorRand.h"
#include "../discrete/BernoulliRand.h"
#include "../continuous/UniformRand.h"

double CantorRand::table[CantorRand::n] = {0};
const bool CantorRand::dummy = CantorRand::SetupTable();

bool CantorRand::SetupTable()
{
    table[0] = 0.33333333333333333333;
    for (int i = 1; i != n; ++i)
        table[i] = table[i - 1] / 3.0;
    return true;
}

CantorRand::CantorRand()
{
}

String CantorRand::Name() const
{
    return "Cantor";
}

double CantorRand::F(const double & x) const
{
    if (x <= 0.0)
        return 0.0;
    if (x >= 1.0)
        return 1.0;
    double a = 0.0, b = 1.0;
    double v = 0.5, d = 0.5, delta = 1.0;
    while (delta > MIN_POSITIVE)
    {
        delta = (b - a) / 3.0;
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

double CantorRand::Variate() const
{
    long double sum = 0.0;
    for (int i = 0; i != n; ++i) {
        sum += table[i] * BernoulliRand::StandardVariate(this->localRandGenerator);
    }
    return sum + sum;
}

long double CantorRand::Mean() const
{
    return 0.5;
}

long double CantorRand::Variance() const
{
    return 0.125;
}

double CantorRand::quantileImpl(double p, double initValue) const
{
    if (RandMath::findRoot<double>([this, p] (double x)
    {
        return F(x) - p;
    }, 0.0, 1.0, initValue))
        return initValue;
    return NAN;
}

double CantorRand::quantileImpl(double p) const
{
    return quantileImpl(p, p);
}

double CantorRand::quantileImpl1m(double p, double initValue) const
{
    if (RandMath::findRoot<double>([this, p] (double x)
    {
        return S(x) - p;
    }, 0.0, 1.0, initValue))
        return initValue;
    return NAN;
}

double CantorRand::quantileImpl1m(double p) const
{
    return quantileImpl1m(p, 1.0 - p);
}

std::complex<double> CantorRand::CFImpl(double t) const
{
    double prod = 1.0;
    for (int i = 0; i != n; ++i)
        prod *= std::cos(table[i]);
    std::complex<double> y(0.0, 0.5 * t);
    y = std::exp(y);
    return y * prod;
}

double CantorRand::Median() const
{
    return 1.0 / 3;
}

long double CantorRand::Skewness() const
{
    return 0.0l;
}

long double CantorRand::ExcessKurtosis() const
{
    return -1.6l;
}
