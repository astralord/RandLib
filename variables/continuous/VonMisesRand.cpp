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
    k = concentration;
    if (k <= 0)
        k = MIN_POSITIVE;
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
    mu - M_PI, x);
}

double VonMisesRand::variate() const
{
    //TODO:
    return 0.0;
}

double VonMisesRand::Mean() const
{
    return mu;
}

double VonMisesRand::Variance() const
{
    double m2 = RandMath::integral([this] (double t)
    {
        return t * t * VonMisesRand::f(t);
    },
    mu - M_PI, mu + M_PI);
    return m2 - mu * mu;
}

double VonMisesRand::Median() const
{
    return mu;
}

double VonMisesRand::Mode() const
{
    return mu;
}

