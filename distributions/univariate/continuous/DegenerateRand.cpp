#include "DegenerateRand.h"
#include <algorithm>
#include <functional>

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
    if (!std::isfinite(value))
        throw std::invalid_argument("Degenerate distribution: value should be finite");
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

long double DegenerateRand::Mean() const
{
    return a;
}

long double DegenerateRand::Variance() const
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

long double DegenerateRand::Skewness() const
{
    return NAN;
}

long double DegenerateRand::ExcessKurtosis() const
{
    return NAN;
}

long double DegenerateRand::Entropy() const
{
    return 0.0;
}

void DegenerateRand::Fit(const std::vector<double> &sample)
{
    if (sample.empty())
        throw std::invalid_argument(this->fitErrorDescription(this->TOO_FEW_ELEMENTS, "The sample is empty"));
    if (std::adjacent_find(sample.begin(), sample.end(), std::not_equal_to<>()) != sample.end())
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, "All elements should be equal to each other"));
    SetValue(sample.front());
}
