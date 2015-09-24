#include "SingularRand.h"

SingularRand::SingularRand()
{

}

double SingularRand::Hazard(double) const
{
    return NAN;
}

double SingularRand::ExpectedValue(const std::function<double (double)> &, double) const
{
    return NAN;
}

double SingularRand::Mode() const
{
    return NAN;
}

