#include "PoissonRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

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

double PoissonRand::F(int k) const
{
    return (k < 0) ? 0 : RandMath::upperIncGamma(k + 1, lambda) / RandMath::factorial(k);
}

int PoissonRand::variate() const
{
    double U = UniformRand::standardVariate();
    int k = floorLambda;
    double s = FFloorLambda, p = PFloorLambda;
    if (s < U)
    {
        do {
            ++k;
            p *= lambda / k;
            s += p;
        } while (s < U);
    }
    else
    {
        s -= p;
        while (s > U) {
            p /= lambda / k;
            --k;
            s -= p;
        }
    }
    return k;
}

int PoissonRand::variate(double rate)
{
    int k = -1;
    double s = 0;
    do {
        s += ExponentialRand::standardVariate();
        ++k;
    } while (s < rate);
    return k;
}

double PoissonRand::Mean() const
{
    return lambda;
}

double PoissonRand::Variance() const
{
    return lambda;
}

std::complex<double> PoissonRand::CF(double t) const
{
    if (t == 0)
        return std::complex<double>(1, 0);
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

int PoissonRand::Mode() const
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

bool PoissonRand::checkValidity(const std::vector<double> &sample)
{
    for (int var : sample) {
        if (var < 0)
            return false;
    }
    return true;
}

bool PoissonRand::fitRateMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setRate(RandMath::sampleMean(sample));
    return true;
}

bool PoissonRand::fitRateBayes(const std::vector<double> &sample, GammaRand &priorDistribution)
{
    int n = sample.size();
    if (n <= 0 || !checkValidity(sample))
        return false;
    double alpha = priorDistribution.getShape();
    double beta = priorDistribution.getRate();
    priorDistribution.setParameters(alpha + RandMath::sum(sample), beta + n);
    setRate(priorDistribution.Mean());
    return true;
}
