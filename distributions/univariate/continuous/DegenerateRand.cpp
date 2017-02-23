#include "DegenerateRand.h"

DegenerateRand::DegenerateRand(double value)
{
    SetValue(value);
}

std::string DegenerateRand::Name() const
{
    return "Degenerate";
}

void DegenerateRand::SetValue(double value)
{
    a = value;
}

double DegenerateRand::f(double x) const
{
    return (x == a) ? INFINITY : 0.0;
}

double DegenerateRand::logf(double x) const
{
    return (x == a) ? INFINITY : -INFINITY;
}

double DegenerateRand::F(double x) const
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

bool DegenerateRand::FitMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    SetValue(sample.at(0));
    return true;
}

bool DegenerateRand::FitMM(const std::vector<double> &sample)
{
    return FitMLE(sample);
}
