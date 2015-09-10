#include "CauchyRand.h"

CauchyRand::CauchyRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string CauchyRand::name()
{
    return "Cauchy(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void CauchyRand::setLocation(double location)
{
    x0 = location;
}

void CauchyRand::setScale(double scale)
{
    gamma = std::max(scale, MIN_POSITIVE);
    gammaInv = 1.0 / gamma;
}

double CauchyRand::f(double x) const
{
    double y = x - x0;
    y *= y;
    y *= gammaInv;
    y += gamma;
    return M_1_PI / y;
}

double CauchyRand::F(double x) const
{
    double y = x - x0;
    y *= gammaInv;
    y = std::atan(y);
    y *= M_1_PI;
    return y + .5;
}

double CauchyRand::variate() const
{
    return x0 + gamma * standardVariate();
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

double E() const
{ 
    return NAN;
}

double Var() const
{
    return INFINITY; 
}
    
double quantile(double p) const
{
    if (p < 0 || p > 1)
         return NAN;
    if (p == 0)
        return -INFINITY;
    if (p == 1)
        return INFINITY;
    return x0 + gamma * std::tan(M_PI * (p - 0.5));
}

double Median() const
{
    return x0;
}

double Mode() const
{
    return x0;
}

double Skewness() const
{
    return NAN;
}

double ExcessKurtosis() const
{
    return NAN;
}
    
double CauchyRand::Entropy() const
{
     return std::log(4 * gamma * M_PI);
}
