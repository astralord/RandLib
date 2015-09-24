#include "ZipfRand.h"
#include "../continuous/ParetoRand.h"

ZipfRand::ZipfRand(double exponent, int number)
{
    setParameters(exponent, number);
}

std::string ZipfRand::name()
{
    return "Zipf(" + toStringWithPrecision(getExponent()) + ", "
                   + toStringWithPrecision(getNumber()) + ")";
}

void ZipfRand::setParameters(double exponent, int number)
{
    s = std::max(exponent, 1.0);
    N = number < 1 ? 1 : number;

    invHarmonicNumber = RandMath::harmonicNumber(s, N);
    b = std::pow(2.0, s - 1.0);
}

double ZipfRand::P(int k) const
{
    if (k < 1 || k > N)
        return 0.0;
    return std::pow(k, -s) * invHarmonicNumber;
}

double ZipfRand::F(double x) const
{
    if (x < 1.0)
        return 0.0;
    if (x >= N)
        return 1.0;
    int k = static_cast<int>(std::floor(x));
    return RandMath::harmonicNumber(s, k) * invHarmonicNumber;
}

double ZipfRand::variate() const
{
    //TODO - maybe rejection from zeta
    return -1.0;
}

double ZipfRand::Mean() const
{
    return RandMath::harmonicNumber(s - 1, N) * invHarmonicNumber;
}

double ZipfRand::Variance() const
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

std::complex<double> ZipfRand::CF(double t) const
{
    std::complex<double> sum(0.0, 0.0);
    for (int i = 1; i <= N; ++i)
    {
        std::complex<double> numerator(0, i * t);
        numerator = std::exp(numerator);
        double denominator = std::pow(i, s);
        sum += numerator / denominator;
    }
    return invHarmonicNumber * sum;
}

double ZipfRand::Skewness() const
{
    double harmonic0 = 1.0 / invHarmonicNumber;
    double harmonic1 = RandMath::harmonicNumber(s - 1, N);
    double harmonic2 = RandMath::harmonicNumber(s - 2, N);
    double harmonic3 = RandMath::harmonicNumber(s - 3, N);

    double first = harmonic3 * harmonic0 * harmonic0;
    double harmonic2prod0 = harmonic2 * harmonic0;
    double second = -3 * harmonic2 * harmonic2prod0;
    double harmonic1Sq = harmonic1 * harmonic1;
    double third = harmonic1 * harmonic1Sq;
    third += third;

    double numerator = first + second + third;
    double denominator = harmonic2prod0 - harmonic1Sq;
    denominator *= std::sqrt(denominator);

    return numerator / denominator;
}

double ZipfRand::ExcessKurtosis() const
{
    double harmonic0 = 1.0 / invHarmonicNumber;
    double harmonic1 = RandMath::harmonicNumber(s - 1, N);
    double harmonic2 = RandMath::harmonicNumber(s - 2, N);
    double harmonic3 = RandMath::harmonicNumber(s - 3, N);
    double harmonic4 = RandMath::harmonicNumber(s - 4, N);

    double harmonic2prod0 = harmonic2 * harmonic0;
    double harmonic0Sq = harmonic0 * harmonic0;
    double harmonic1Sq = harmonic1 * harmonic1;

    double denominator = harmonic2prod0 - harmonic1Sq;
    denominator *= denominator;

    double numerator = harmonic0Sq * harmonic0 * harmonic4;
    numerator -= 4 * harmonic0Sq * harmonic1 * harmonic3;
    numerator += 6 * harmonic2prod0 * harmonic1Sq;
    numerator -= 3 * harmonic1Sq * harmonic1Sq;

    return numerator / denominator - 3;
}
