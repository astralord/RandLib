#include "VonMisesRand.h"
#include "UniformRand.h"

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
        k = 1.0;
    I0kInv = 1.0 / RandMath::modifiedBesselFirstKind(k, 0);
    if (k > 1.3)
        s = 1.0 / std::sqrt(k);
    else
        s = M_PI * std::exp(-k);
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
    /// Generating von Mises variates by the ratio-of-uniforms method
    /// Lucio Barabesi. Dipartimento di Metodi Quantitativi, Universiteta di Siena
    int iter = 0;
    do {
        double U = UniformRand::standardVariate(), V = UniformRand::variate(-1, 1);
        double theta = s * V / U;
        if ((std::fabs(theta) <= M_PI) &&
           ((k * theta * theta < 4.0 * (1.0 - U)) || (k * std::cos(theta) >= 2 * std::log(U) + k)))
            return mu + theta;
    } while (++iter < 1e9);
    return NAN; /// fail
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

std::complex<double> VonMisesRand::CF(double t) const
{
    // TODO: verify this one
    int n = std::fabs(std::floor(t));
    std::complex<double> y(0.0, n * mu);
    y = std::exp(y);
    return I0kInv * y * RandMath::modifiedBesselFirstKind(k, n);
}

double VonMisesRand::Median() const
{
    return mu;
}

double VonMisesRand::Mode() const
{
    return mu;
}

