#include "VonMisesRand.h"

VonMisesRand::VonMisesRand(double location, double concentration)
{
    setLocation(location);
    setConcentration(concentration);
}

std::string VonMisesRand::name()
{
    return "von Mises(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getConcentration()) + ")";
}

void VonMisesRand::setLocation(double location)
{
    mu = location;
}

void VonMisesRand::setConcentration(double concentration)
{
    k = std::max(concentration, 0.0);
    I0kInv = 1.0 / RandMath::modifiedBesselFirstKind(k, 0);
}

double VonMisesRand::f(double x) const
{
    if (x < mu - M_PI || x > mu + M_PI)
        return 0.0;
    return 0.5 * M_1_PI * I0kInv * std::exp(k * std::cos(x - mu));
}

double VonMisesRand::F(double x) const
{
    if (x <= mu - M_PI)
        return 0;
    if (x >= mu + M_PI)
        return 1;
    return RandMath::integral([this] (double t)
    {
        return VonMisesRand::f(t);
    },
    0, x);
}

double VonMisesRand::variate() const
{
    //TODO:
    return 0.0;
}

double VonMisesRand::E() const
{
    return mu;
}

double VonMisesRand::Var() const
{
    //CIRCULAR!
    return 1.0 - RandMath::modifiedBesselFirstKind(k, 1) * I0kInv;
}

double VonMisesRand::Median() const
{
    return mu;
}

double VonMisesRand::Mode() const
{
    return mu;
}

