#include "PoissonRand.h"

PoissonRand::PoissonRand(double rate)
{
    setRate(rate);
}

void PoissonRand::setName()
{
    nameStr = "Poisson(" + toStringWithPrecision(getRate()) + ")";
}

void PoissonRand::setRate(double rate)
{
    l = std::max(rate, MIN_POSITIVE);
    exp_l = std::exp(-l);

    setName();
}

double PoissonRand::P(int k) const
{
    return (k >= 0) ? exp_l * std::pow(l, k) / RandMath::fastFactorial(k) : 0;
}

double PoissonRand::F(double x) const
{
    if (x <= 0)
        return 0;
    int k = std::floor(x);
    return RandMath::upperIncGamma(k + 1, l) / RandMath::fastFactorial(k);
}

double PoissonRand::variate() const
{
    int y = 0;
    double u = UniformRand::standardVariate();
    double p = exp_l, s = p;
    while (s < u) {
        ++y;
        p *= l / y;
        s += p;
    }
    return y;
}
