#include "HyperGeometricRand.h"

HyperGeometricRand::HyperGeometricRand()
{

}

std::string HyperGeometricRand::name()
{
    return "Hyper-Geometric(" + toStringWithPrecision(1) + ")";
}

double HyperGeometricRand::P(int k) const
{
    return 1.0 * k;
}

double HyperGeometricRand::F(double x) const
{
    return x;
}

double HyperGeometricRand::variate() const
{
    return 1.0;
}

