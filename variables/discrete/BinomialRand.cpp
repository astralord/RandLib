#include "BinomialRand.h"

BinomialRand::BinomialRand(int number, double probability)
{
    setNumber(number);
    setProbability(probability);
}

std::string BinomialRand::name()
{
    return "Binomial(" + toStringWithPrecision(getNumber()) + ", " + toStringWithPrecision(getProbability()) + ")";
}

void BinomialRand::setNumber(int number)
{
    n = std::max(number, 1);
}

void BinomialRand::setProbability(double probability)
{
    p = std::min(std::max(probability, MIN_POSITIVE), 1.0);
    B.setProbability(p);
}

double BinomialRand::P(int k) const
{
    if (k < 0 || k > n)
        return 0;
    return RandMath::binomialCoef(n, k) * std::pow(p, k) * std::pow(1 - p, n - k);
}

double BinomialRand::F(double x) const
{
    if (x < 0)
        return 0;
    return x;
}

double BinomialRand::variate() const
{
    unsigned sum = 0;
    for (int i = 0; i != n; ++i)
        sum += B.variate();
    return sum;
}

std::complex<double> BinomialRand::CF(double t) const
{
    return std::pow(B.CF(t), n);
}
