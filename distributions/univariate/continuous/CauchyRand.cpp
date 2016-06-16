#include "CauchyRand.h"
#include "UniformRand.h"

CauchyRand::CauchyRand(double location, double scale)
    : StableRand(1, 0, scale, location)
{
}

std::string CauchyRand::name() const
{
    return "Cauchy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

double CauchyRand::f(double x) const
{
    return StableRand::pdfCauchy(x);
}

double CauchyRand::F(double x) const
{
    return StableRand::cdfCauchy(x);
}

double CauchyRand::variate() const
{
    return mu + sigma * standardVariate();
}

double CauchyRand::variate(double location, double scale)
{
    return location + scale * standardVariate();
}

double CauchyRand::standardVariate()
{
    double x, y;
    do {
        x = UniformRand::variate(-1, 1);
        y = UniformRand::variate(-1, 1);
    } while (x * x + y * y > 1.0 || y == 0.0);
    return x / y;
}
    
std::complex<double> CauchyRand::CF(double t) const
{
    std::complex<double> x(-sigma * std::fabs(t), mu * t);
    return std::exp(x);
}

double CauchyRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
         return NAN;
    return mu + sigma * std::tan(M_PI * (p - 0.5));
}

double CauchyRand::Entropy() const
{
    return std::log(4 * sigma * M_PI);
}
