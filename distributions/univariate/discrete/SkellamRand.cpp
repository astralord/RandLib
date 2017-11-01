#include "SkellamRand.h"
#include "../continuous/NoncentralChiSquaredRand.h"

SkellamRand::SkellamRand(double rate1, double rate2)
{
    SetRates(rate1, rate2);
}

String SkellamRand::Name() const
{
    return "Skellam(" + toStringWithPrecision(GetFirstRate()) + ", " + toStringWithPrecision(GetSecondRate()) + ")";
}

void SkellamRand::SetRates(double rate1, double rate2)
{
    if (rate1 <= 0.0 || rate2 <= 0.0)
        throw std::invalid_argument("Skellam distribution: rates should be positive");

    X.SetRate(rate1);
    mu1 = X.GetRate();
    logMu1 = std::log(mu1);
    sqrtMu1 = std::sqrt(mu1);

    Y.SetRate(rate2);
    mu2 = Y.GetRate();
    logMu2 = std::log(mu2);
    sqrtMu2 = std::sqrt(mu2);
}

double SkellamRand::P(const int & k) const
{
    return std::exp(logP(k));
}

double SkellamRand::logP(const int & k) const
{
    double y = RandMath::logBesselI(k, 2 * sqrtMu1 * sqrtMu2);
    y += 0.5 * k * (logMu1 - logMu2);
    y -= mu1 + mu2;
    return y;
}

double SkellamRand::F(const int & k) const
{
    return (k < 0) ? RandMath::MarcumP(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2) : RandMath::MarcumQ(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
}

double SkellamRand::S(const int & k) const
{
    return (k < 0) ? RandMath::MarcumQ(-k, mu1, mu2, sqrtMu1, sqrtMu2, logMu1, logMu2) : RandMath::MarcumP(k + 1, mu2, mu1, sqrtMu2, sqrtMu1, logMu2, logMu1);
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
