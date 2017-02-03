#include "SkellamRand.h"
#include "../continuous/NoncentralChiSquaredRand.h"

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
    X.SetRate(mean1);
    mu1 = X.GetRate();

    Y.SetRate(mean2);
    mu2 = Y.GetRate();

    pmfCoef2 = std::log(mu1 / mu2);
    pmfCoef1 = 2.0 * std::sqrt(mu1 * mu2);
}

double SkellamRand::P(int k) const
{
    double y = RandMath::logModifiedBesselFirstKind(pmfCoef1, k);
    y += 0.5 * k * pmfCoef2;
    y -= mu1 + mu2;
    return std::exp(y);
}

double SkellamRand::F(int k) const
{
    return (k < 0) ? RandMath::MarcumP(-k, mu1, mu2) : RandMath::MarcumQ(k + 1, mu2, mu1);
}

double SkellamRand::S(int k) const
{
    return (k < 0) ? RandMath::MarcumQ(-k, mu1, mu2) : RandMath::MarcumP(k + 1, mu2, mu1);
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
    double cosT = std::cos(t), sinT = std::sin(t);
    double x = (cosT - 1) * (mu1 + mu2);
    double y = sinT * (mu1 - mu2);
    std::complex<double> z(x, y);
    return std::exp(z);
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
