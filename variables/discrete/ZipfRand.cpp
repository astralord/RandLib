#include "ZipfRand.h"

ZipfRand::ZipfRand(double exponent, size_t number)
{
    setParameters(exponent, number);
}

std::string ZipfRand::name()
{
    return "Zipf(" + toStringWithPrecision(getExponent()) + ", "
                   + toStringWithPrecision(getNumber()) + ")";
}

double ZipfRand::harmonicNumber(double exponent, size_t number)
{
    if (exponent < 1 || number < 1)
        return 0;
    double res = 1.0;
    for (size_t i = 1; i != number; ++i)
        res += std::pow(i + 1, -exponent);
    return res;
}

void ZipfRand::setParameters(double exponent, size_t number)
{
    s = std::max(exponent, 1.0);
    N = number < 1 ? 1 : number;

    denominator = harmonicNumber(s, N);
}

double ZipfRand::P(int k) const
{
    return std::pow(k, -s) * denominator;
}

double ZipfRand::F(double x) const
{
    return harmonicNumber(s, std::floor(x)) * denominator;
}

double ZipfRand::variate() const
{
    return 1.0;
}

