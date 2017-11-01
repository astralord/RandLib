#include "LevyRand.h"
#include "NormalRand.h"

LevyRand::LevyRand(double location, double scale)
    : StableDistribution(0.5, 1, scale, location)
{
}

String LevyRand::Name() const
{
    return "Levy(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

double LevyRand::f(const double & x) const
{
    return pdfLevy(x);
}

double LevyRand::logf(const double & x) const
{
    return logpdfLevy(x);
}

double LevyRand::F(const double & x) const
{
    return cdfLevy(x);
}

double LevyRand::S(const double & x) const
{
    return cdfLevyCompl(x);
}

double LevyRand::Variate() const
{
    double rv = NormalRand::StandardVariate();
    rv *= rv;
    rv = gamma / rv;
    return mu + rv;
}

double LevyRand::StandardVariate()
{
    double rv = NormalRand::StandardVariate();
    return 1.0 / (rv * rv);
}

double LevyRand::quantileImpl(double p) const
{
    return quantileLevy(p);
}

double LevyRand::quantileImpl1m(double p) const
{
    return quantileLevy1m(p);
}

std::complex<double> LevyRand::CFImpl(double t) const
{
    return cfLevy(t);
}

void LevyRand::FitScale(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNotLessThan(mu, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(mu)));
    long double invSum = 0.0;
    for (double var : sample)
        invSum += 1.0 / (var - mu);
    invSum = 1.0 / invSum;
    SetScale(sample.size() * invSum);
}
