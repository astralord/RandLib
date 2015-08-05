#include "GeometricRand.h"

GeometricRand::GeometricRand(double probability)
{
    setProbability(probability);
}

void GeometricRand::setProbability(double probability)
{
    p = probability;
    Exp.setRate(-std::log(1 - p));
}

double GeometricRand::P(int k) const
{
    return p * std::pow(1 - p, k);
}

double GeometricRand::F(double x) const
{
    int k = static_cast<int>(x);
    return 1 - std::pow(1 - p, k + 1);
}

int GeometricRand::variate()
{
    double x = Exp.variate();
    return static_cast<int>(std::floor(x));
}

