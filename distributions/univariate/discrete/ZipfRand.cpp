#include "ZipfRand.h"
#include "../continuous/UniformRand.h"

ZipfRand::ZipfRand(double exponent, int number)
{
    SetParameters(exponent, number);
}

std::string ZipfRand::Name() const
{
    return "Zipf(" + toStringWithPrecision(GetExponent()) + ", "
                   + toStringWithPrecision(GetNumber()) + ")";
}

void ZipfRand::SetParameters(double exponent, int number)
{
    s = exponent;
    if (s <= 1.0)
        s = 2.0;
    n = number < 1 ? 1 : number;

    invHarmonicNumber = 1.0 / RandMath::harmonicNumber(s, n);

    // WARNING: we calculate pow here and in invHarmonic number
    hashedVarNum = n > tableSize ? tableSize : n;
    table[0] = 1.0;
    for (int i = 1; i < hashedVarNum - 1; ++i)
        table[i] = table[i - 1] + std::pow(i + 1, -s);
    if (hashedVarNum == n)
        table[hashedVarNum - 1] = 1.0 / invHarmonicNumber;
    else
        table[hashedVarNum - 1] = table[hashedVarNum - 2] + std::pow(hashedVarNum, -s);
    for (int i = 0; i < hashedVarNum; ++i)
        table[i] *= invHarmonicNumber;
}

double ZipfRand::P(int k) const
{
    if (k < 1 || k > n)
        return 0.0;
    return std::pow(k, -s) * invHarmonicNumber;
}

double ZipfRand::F(int k) const
{
    if (k < 1.0)
        return 0.0;
    if (k >= n)
        return 1.0;
    return RandMath::harmonicNumber(s, k) * invHarmonicNumber;
}

int ZipfRand::Variate() const
{
    double U = UniformRand::StandardVariate();
    int k = 1;

    /// if we didn't manage to hash values for such U
    if (U > table[hashedVarNum - 1]) {
        k = hashedVarNum;
        double sum = table[hashedVarNum - 1];
        do {
            ++k;
            sum += std::pow(k, -s) * invHarmonicNumber;
        } while (sum < U);
    }
    else {
        while (k <= hashedVarNum && table[k - 1] < U)
            ++k;
    }
    return k;
}

double ZipfRand::Mean() const
{
    return RandMath::harmonicNumber(s - 1, n) * invHarmonicNumber;
}

double ZipfRand::Variance() const
{
    double numerator = RandMath::harmonicNumber(s - 1, n);
    numerator *= numerator;
    numerator = RandMath::harmonicNumber(s - 2, n) * RandMath::harmonicNumber(s, n) - numerator;
    return numerator * invHarmonicNumber * invHarmonicNumber;
}

int ZipfRand::Mode() const
{
    return 1.0;
}

std::complex<double> ZipfRand::CF(double t) const
{
    if (t == 0)
        return 1;
    std::complex<double> sum(0.0, 0.0);
    for (int i = 1; i <= n; ++i)
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
    double harmonic1 = RandMath::harmonicNumber(s - 1, n);
    double harmonic2 = RandMath::harmonicNumber(s - 2, n);
    double harmonic3 = RandMath::harmonicNumber(s - 3, n);

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
    double harmonic1 = RandMath::harmonicNumber(s - 1, n);
    double harmonic2 = RandMath::harmonicNumber(s - 2, n);
    double harmonic3 = RandMath::harmonicNumber(s - 3, n);
    double harmonic4 = RandMath::harmonicNumber(s - 4, n);

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
