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
    lambda = rate;
    if (lambda <= 0)
        lambda = 1.0;

    expmLambda = std::exp(-lambda);
    logLambda = std::log(lambda);

    floorLambda = std::floor(lambda);
    FFloorLambda = F(floorLambda);
    PFloorLambda = P(floorLambda);
}

double PoissonRand::P(int k) const
{
    if (k < 0)
        return 0;
    return (k < 20) ? expmLambda * std::pow(lambda, k) / RandMath::factorial(k) : std::exp(k * logLambda - lambda - std::lgamma(k + 1));
}

double PoissonRand::F(double x) const
{
    if (x < 0)
        return 0;
    double k = std::floor(x);
    return RandMath::upperIncGamma(k + 1, lambda) / RandMath::factorial(k);
}

double PoissonRand::variate() const
{
    double U = UniformRand::standardVariate();
    int y = floorLambda;
    double s = FFloorLambda, p = PFloorLambda;
    if (s < U)
    {
        do {
            ++y;
            p *= lambda / y;
            s += p;
        } while (s < U);
    }
    else
    {
        s -= p;
        while (s > U) {
            p /= lambda / y;
            --y;
            s -= p;
        }
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

std::complex<double> PoissonRand::CF(double t) const
{
    std::complex<double> y(0.0, t);
    y = std::exp(y) - 1.0;
    y *= lambda;
    return std::exp(y);
}

double PoissonRand::Median() const
{
    /// this value is approximate
    return std::floor(lambda + 1.0 / 3 - 0.02 / lambda);
}

double PoissonRand::Mode() const
{
    return std::floor(lambda);
}

double PoissonRand::Skewness() const
{
    return 1.0 / std::sqrt(lambda);
}

double PoissonRand::ExcessKurtosis() const
{
    return 1.0 / lambda;
}

bool PoissonRand::fitToData(const QVector<int> &sample)
{
    int N = sample.size();
    if (N == 0)
        return false;

    /// Calculate mu
    long double sum = 0.0L;
    for (int var : sample) {
        if (var < 0)
            return false;
        sum += var;
    }

    setRate(sum / N);
    return true;
}
