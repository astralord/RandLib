#include "PoissonRand.h"

PoissonRand::PoissonRand(double rate) :
    U(0, 1)
{
    setRate(rate);
}

void PoissonRand::setRate(double rate)
{
    l = std::max(rate, MIN_POSITIVE);
    exp_l = std::exp(-l);
}

double PoissonRand::P(int k)
{
    return (k > 0) ? exp_l * std::pow(l, k) / factorial(k) : 0;
}

double PoissonRand::cdf(double x)
{
    return x;
}

double PoissonRand::value()
{
    double y = 0;
    double u = U.value();
    double p = exp_l, s = p;
    while (s < u) {
        ++y;
        p *= l / y;
        s += p;
    }
    return y;
}
