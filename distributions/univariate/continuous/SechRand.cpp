#include "SechRand.h"
#include "CauchyRand.h"

SechRand::SechRand()
{
}

String SechRand::Name() const
{
    return "Hyperbolic secant";
}

double SechRand::f(const double & x) const
{
    return 0.5 / std::cosh(M_PI_2 * x);
}

double SechRand::logf(const double & x) const
{
    return M_PI_2 * x - RandMath::log1pexp(M_PI * x);
}

double SechRand::F(const double & x) const
{
    double y = std::exp(M_PI_2 * x);
    return M_2_PI * RandMath::atan(y);
}

double SechRand::Variate() const
{
    double y = std::fabs(CauchyRand::StandardVariate());
    return M_2_PI * std::log(y);
}

double SechRand::Mean() const
{
    return 0.0;
}

double SechRand::Variance() const
{
    return 1.0;
}

std::complex<double> SechRand::CFImpl(double t) const
{
    return 1.0 / std::cosh(t);
}

double SechRand::quantileImpl(double p) const
{
    double x = M_PI_2 * p;
    x = std::tan(x);
    x = std::log(x);
    return M_2_PI * x;
}

double SechRand::quantileImpl1m(double p) const
{
    double x = M_PI_2 * p;
    x = std::tan(x);
    x = -std::log(x);
    return M_2_PI * x;
}

double SechRand::Median() const
{
    return 0.0;
}

double SechRand::Mode() const
{
    return 0.0;
}

double SechRand::Skewness() const
{
    return 0.0;
}

double SechRand::ExcessKurtosis() const
{
    return 2.0;
}

double SechRand::Entropy() const
{
    return 2.0 * M_2_PI * M_CATALAN;
}

