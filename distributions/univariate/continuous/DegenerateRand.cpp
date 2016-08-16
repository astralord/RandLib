#include "DegenerateRand.h"

DegenerateRand::DegenerateRand(double value)
{
    setValue(value);
}

std::string DegenerateRand::name() const
{
    return "Degenerate";
}

void DegenerateRand::setValue(double value)
{
    a = value;
}

double DegenerateRand::f(double x) const
{
    return (x == a) ? INFINITY : 0.0;
}

double DegenerateRand::F(double x) const
{
    return (x < a) ? 0.0 : 1.0;
}

double DegenerateRand::variate() const
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

std::complex<double> DegenerateRand::CF(double t) const
{
    double re = std::cos(a * t);
    double im = std::sin(a * t);
    return std::complex<double>(re, im);
}

double DegenerateRand::QuantileImpl(double) const
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

bool DegenerateRand::fitMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setValue(sample.at(0));
    return true;
}

bool DegenerateRand::fitMM(const std::vector<double> &sample)
{
    return fitMLE(sample);
}
