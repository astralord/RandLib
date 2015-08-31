#include "PoissonRand.h"

PoissonRand::PoissonRand(double rate)
{
    setRate(rate);
}

std::string PoissonRand::name()
{
    return "Poisson(" + toStringWithPrecision(getRate()) + ")";
}

void PoissonRand::setRate(double rate)
{
    lambda = std::max(rate, MIN_POSITIVE);
    expLambda = std::exp(-lambda);
}

double PoissonRand::P(int k) const
{
    return (k >= 0) ? expLambda * std::pow(lambda, k) / RandMath::factorial(k) : 0;
}

double PoissonRand::F(double x) const
{
    if (x <= 0)
        return 0;
    int k = std::floor(x);
    return RandMath::upperIncGamma(k + 1, lambda) / RandMath::factorial(k);
}

std::complex<double> PoissonRand::CF(double t) const
{
    std::complex<double> y(0.0, t);
    y = std::exp(y) - 1.0;
    y *= lambda;
    return std::exp(y);
}

double PoissonRand::variate() const
{
    int y = 0;
    double u = UniformRand::standardVariate();
    double p = expLambda, s = p;
    while (s < u) {
        ++y;
        p *= lambda / y;
        s += p;
    }
    return y;
}

double PoissonRand::variate(double rate)
{
    int y = -1;
    double s = 0;
    do {
        s += ExponentialRand::standardVariate();
        ++y;
    } while (s < rate);
    return y;
}
