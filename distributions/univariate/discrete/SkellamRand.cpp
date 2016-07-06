#include "SkellamRand.h"

SkellamRand::SkellamRand(double mean1, double mean2)
{
    setMeans(mean1, mean2);
}

std::string SkellamRand::name() const
{
    return "Skellam(" + toStringWithPrecision(getFirstMean()) + ", " + toStringWithPrecision(getSecondMean()) + ")";
}

void SkellamRand::setMeans(double mean1, double mean2)
{
    mu1 = mean1;
    X.setRate(mu1);

    mu2 = mean2;
    Y.setRate(mu2);

    pmfCoef1 = std::exp(-mu1 - mu2);
    pmfCoef2 = std::sqrt(mu1 / mu2);
    pmfCoef3 = 2.0 * std::sqrt(mu1 * mu2);
}

double SkellamRand::P(int k) const
{
    double y = pmfCoef1;
    y *= std::pow(pmfCoef2, k);
    y *= RandMath::modifiedBesselFirstKind(pmfCoef3, k);
    return y;
}

double SkellamRand::F(int k) const
{
    int i = std::max(0, -k);
    double sum = 0, summand = 0.0, xCDF = 0;
    do {
        xCDF = X.F(k + i);
        summand = xCDF * Y.P(i);
        sum += summand;
        ++i;
    } while (summand > MIN_POSITIVE || xCDF < 0.9);
    return sum;
}

int SkellamRand::variate() const
{
    return X.variate() - Y.variate();
}

double SkellamRand::Mean() const
{
    return mu1 - mu2;
}

double SkellamRand::Variance() const
{
    return mu1 + mu2;
}

std::complex<double> SkellamRand::CF(double t) const
{
    if (t == 0)
        return 1;
    std::complex<double> x(0, t);
    x = std::exp(x);
    x *= mu1;

    std::complex<double> y(0, -t);
    y = std::exp(y);
    y *= mu2;

    y = std::exp(y + x);
    return pmfCoef1 * y;
}

int SkellamRand::Mode() const
{
    return Mean();
}

double SkellamRand::Skewness() const
{
    return (mu1 - mu2) / std::pow(mu1 + mu2, 1.5);
}

double SkellamRand::ExcessKurtosis() const
{
    return 1.0 / (mu1 + mu2);
}
