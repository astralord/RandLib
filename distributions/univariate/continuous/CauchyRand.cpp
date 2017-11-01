#include "CauchyRand.h"
#include "UniformRand.h"

CauchyRand::CauchyRand(double location, double scale)
    : StableDistribution(1, 0, scale, location)
{
}

String CauchyRand::Name() const
{
    return "Cauchy(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

double CauchyRand::f(const double & x) const
{
    return pdfCauchy(x);
}

double CauchyRand::F(const double & x) const
{
    return cdfCauchy(x);
}

double CauchyRand::S(const double & x) const
{
    return cdfCauchyCompl(x);
}

double CauchyRand::Variate() const
{
    return mu + gamma * StandardVariate();
}

double CauchyRand::StandardVariate()
{
    double x, y;
    do {
        x = UniformRand::Variate(-1, 1);
        y = UniformRand::Variate(-1, 1);
    } while (y == 0.0 || x * x + y * y > 1.0);
    return x / y;
}
    
std::complex<double> CauchyRand::CFImpl(double t) const
{
    return cfCauchy(t);
}

double CauchyRand::quantileImpl(double p) const
{
    return quantileCauchy(p);
}

double CauchyRand::quantileImpl1m(double p) const
{
    return quantileCauchy1m(p);
}

double CauchyRand::Entropy() const
{
    return 2 * M_LN2 + logGamma + M_LNPI;
}
