#include "LevyRand.h"
#include "NormalRand.h"

LevyRand::LevyRand(double location, double scale)
    : StableRand(0.5, 1, scale, location)
{
}

std::string LevyRand::Name() const
{
    return "Levy(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

double LevyRand::f(double x) const
{
    return StableRand::pdfLevy(x);
}

double LevyRand::F(double x) const
{
    return StableRand::cdfLevy(x);
}

double LevyRand::S(double x) const
{
    return StableRand::cdfLevyCompl(x);
}

double LevyRand::Variate() const
{
    double rv = NormalRand::StandardVariate();
    rv *= rv;
    rv = sigma / rv;
    return mu + rv;
}

double LevyRand::Variate(double location, double scale)
{
    return location + scale * StandardVariate();
}

double LevyRand::StandardVariate()
{
    double rv = NormalRand::StandardVariate();
    return 1.0 / (rv * rv);
}

std::complex<double> LevyRand::CFImpl(double t) const
{
    std::complex<double> y(0.0, -2 * sigma * t);
    y = -std::sqrt(y);
    y += std::complex<double>(0.0, mu * t);
    return std::exp(y);
}

double LevyRand::quantileImpl(double p) const
{
    double y = pdfCoef / RandMath::erfcinv(p);
    double guess = mu + M_PI * y * y;
    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

double LevyRand::quantileImpl1m(double p) const
{
    double y = pdfCoef / RandMath::erfinv(p);
    double guess = mu + M_PI * y * y;
    if (RandMath::findRoot([this, p] (double x)
    {
        double first = p - S(x);
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

bool LevyRand::FitScaleMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    long double invSum = 0.0;
    for (double var : sample)
        invSum += 1.0 / (var - mu);
    invSum = 1.0 / invSum;
    SetScale(sample.size() * invSum);
    return true;
}
