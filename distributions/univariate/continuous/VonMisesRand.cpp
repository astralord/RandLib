#include "VonMisesRand.h"
#include "UniformRand.h"

VonMisesRand::VonMisesRand(double location, double concentration)
{
    SetLocation(location);
    SetConcentration(concentration);
}

std::string VonMisesRand::Name() const
{
    return "von Mises(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetConcentration()) + ")";
}

void VonMisesRand::SetLocation(double location)
{
    mu = location;
}

void VonMisesRand::SetConcentration(double concentration)
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

    double lowerBound, upperBound, shift, direction;
    if (x <= mu - M_PI_2) {
        lowerBound = mu - M_PI;
        upperBound = x;
        shift = 0.0;
        direction = 1.0;
    }
    else if (x <= mu) {
        lowerBound = x;
        upperBound = mu;
        shift = 0.5;
        direction = -1.0;
    }
    else if (x <= mu + M_PI_2) {
        lowerBound = mu;
        upperBound = x;
        shift = 0.5;
        direction = 1.0;
    }
    else {
        lowerBound = x;
        upperBound = mu + M_PI;
        shift = 1.0;
        direction = -1.0;
    }
    return shift + direction * RandMath::integral([this] (double t)
    {
        return VonMisesRand::f(t);
    },
    lowerBound, upperBound);
}

double VonMisesRand::Variate() const
{
    /// Generating von Mises variates by the ratio-of-uniforms method
    /// Lucio Barabesi. Dipartimento di Metodi Quantitativi, Universiteta di Siena
    int iter = 0;
    do {
        double U = UniformRand::StandardVariate(), V = UniformRand::Variate(-1, 1);
        double theta = s * V / U;
        if ((std::fabs(theta) <= M_PI) &&
           ((k * theta * theta < 4.0 * (1.0 - U)) || (k * std::cos(theta) >= 2 * std::log(U) + k)))
            return mu + theta;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail
}

double VonMisesRand::Mean() const
{
    return mu;
}

double VonMisesRand::Variance() const
{
    return 2 * RandMath::integral([this] (double t)
    {
        return t * t * VonMisesRand::f(t + mu);
    },
    0, M_PI);
}

std::complex<double> VonMisesRand::CFImpl(double t) const
{
    std::complex<double> y(0.0, t * mu);
    y = std::exp(y);
    return I0kInv * y * RandMath::modifiedBesselFirstKind(k, std::fabs(t));
}

double VonMisesRand::Median() const
{
    return mu;
}

double VonMisesRand::Mode() const
{
    return mu;
}

