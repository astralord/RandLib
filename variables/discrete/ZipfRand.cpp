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

void ZipfRand::setParameters(double exponent, size_t number)
{
    s = std::max(exponent, 1.0);
    N = number < 1 ? 1 : number;

    invHarmonicNumber = RandMath::harmonicNumber(s, N);
}

double ZipfRand::P(int k) const
{
    return std::pow(k, -s) * invHarmonicNumber;
}

double ZipfRand::F(double x) const
{
    return RandMath::harmonicNumber(s, std::floor(x)) * invHarmonicNumber;
}

double ZipfRand::variate() const
{
    return 1.0;
}

double ZipfRand::E() const
{
    return RandMath::harmonicNumber(s - 1, N) * invHarmonicNumber;
}

double ZipfRand::Var() const
{
    double numerator = RandMath::harmonicNumber(s - 1, N);
    numerator *= numerator;
    numerator = RandMath::harmonicNumber(s - 2, N) * RandMath::harmonicNumber(s, N) - numerator;
    return numerator * invHarmonicNumber * invHarmonicNumber;
}

double ZipfRand::Mode() const
{
    return 1.0;
}
