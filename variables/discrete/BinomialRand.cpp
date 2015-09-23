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
    p = std::min(probability, 1.0);
    p = std::max(p, 0.0);
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
    // TODO: find max for k such that F(x) == 1, or there can occure numerical error
    int k = std::floor(x);
    return RandMath::regularizedBetaFun(1 - p, n - k, 1 + k);
}

double BinomialRand::variate() const
{
    int sum = 0;
    for (int i = 0; i != n; ++i)
        sum += B.variate();
    return sum;
}

std::complex<double> BinomialRand::CF(double t) const
{
    return std::pow(B.CF(t), n);
}

double BinomialRand::Median() const
{
    return std::round(n * p);
}

double BinomialRand::Mode() const
{
    double mode = (n + 1) * p - 1;
    return std::round(mode);
}

double BinomialRand::Skewness() const
{
    return (1 - p - p) / std::sqrt(n * p * (1 - p));
}

double BinomialRand::ExcessKurtosis() const
{
    return (1 - 6 * p * (1 - p)) / (n * p * (1 - p));
}
