#include "DegenerateRand.h"

DegenerateRand::DegenerateRand(double value)
{
    SetValue(value);
}

String DegenerateRand::Name() const
{
    return "Degenerate";
}

void DegenerateRand::SetValue(double value)
{
    a = value;
}

double DegenerateRand::f(const double & x) const
{
    return (x == a) ? INFINITY : 0.0;
}

double DegenerateRand::logf(const double & x) const
{
    return (x == a) ? INFINITY : -INFINITY;
}

double DegenerateRand::F(const double & x) const
{
    return (x < a) ? 0.0 : 1.0;
}

double DegenerateRand::Variate() const
{
    return a;
}

double DegenerateRand::Mean() const
{
    return a;
}

double DegenerateRand::Variance() const
{
    return 0.0;
}

std::complex<double> DegenerateRand::CFImpl(double t) const
{
    double re = std::cos(a * t);
    double im = std::sin(a * t);
    return std::complex<double>(re, im);
}

double DegenerateRand::quantileImpl(double) const
{
    return a;
}

double DegenerateRand::quantileImpl1m(double) const
{
    return a;
}

double DegenerateRand::Median() const
{
    return a;
}

double DegenerateRand::Mode() const
{
    return a;
}

double DegenerateRand::Skewness() const
{
    return NAN;
}

double DegenerateRand::ExcessKurtosis() const
{
    return NAN;
}

double DegenerateRand::Entropy() const
{
    return 0.0;
}

void DegenerateRand::Fit(const std::vector<double> &sample)
{
    auto sampleBegin = sample.begin();
    if (!std::equal(sampleBegin, sample.end(), sampleBegin))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, "All elements should be equal to each other"));
    SetValue(*sampleBegin);
}
