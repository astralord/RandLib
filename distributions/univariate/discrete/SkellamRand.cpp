#include "SkellamRand.h"
#include "../continuous/NoncentralChiSquared.h"

SkellamRand::SkellamRand(double mean1, double mean2)
{
    SetMeans(mean1, mean2);
}

std::string SkellamRand::Name() const
{
    return "Skellam(" + toStringWithPrecision(GetFirstMean()) + ", " + toStringWithPrecision(GetSecondMean()) + ")";
}

void SkellamRand::SetMeans(double mean1, double mean2)
{
    mu1 = mean1;
    X.SetRate(mu1);

    mu2 = mean2;
    Y.SetRate(mu2);

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
    // TODO: implement marcum q-function
    // and use it instead of creating Y
    NoncentralChiSquared Y;
    if (k <= 0) {
        Y.SetParameters(-2 * k, 2 * mu1);
        return Y.F(2 * mu2);
    }
    Y.SetParameters(2 * k + 2, 2 * mu2);
    return 1.0 - Y.F(2 * mu1);
}

int SkellamRand::Variate() const
{
    return X.Variate() - Y.Variate();
}

void SkellamRand::Sample(std::vector<int> &outputData) const
{
    X.Sample(outputData);
    for (int & var : outputData)
        var -= Y.Variate();
}

double SkellamRand::Mean() const
{
    return mu1 - mu2;
}

double SkellamRand::Variance() const
{
    return mu1 + mu2;
}

std::complex<double> SkellamRand::CFImpl(double t) const
{
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
