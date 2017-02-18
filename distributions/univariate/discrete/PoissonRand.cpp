#include "PoissonRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

PoissonRand::PoissonRand(double rate)
{
    SetRate(rate);
}

std::string PoissonRand::Name() const
{
    return "Poisson(" + toStringWithPrecision(GetRate()) + ")";
}

void PoissonRand::SetRate(double rate)
{
    lambda = (rate <= 0.0) ? 1.0 : rate;

    logLambda = std::log(lambda);

    floorLambda = std::floor(lambda);
    FFloorLambda = F(floorLambda);
    PFloorLambda = P(floorLambda);
}

double PoissonRand::P(int k) const
{
    if (k < 0)
        return 0.0;
    double y = k * logLambda - lambda;
    if (k <= 20)
        return std::exp(y) / RandMath::factorial(k);
    y -= std::lgamma(k + 1);
    return std::exp(y);
}

double PoissonRand::F(int k) const
{
    return (k >= 0.0) ? RandMath::qgamma(k + 1, lambda) : 0.0;
}

double PoissonRand::S(int k) const
{
    return (k >= 0.0) ? RandMath::pgamma(k + 1, lambda) : 1.0;
}

int PoissonRand::Variate() const
{
    double U = UniformRand::StandardVariate();
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

int PoissonRand::Variate(double rate)
{
    int k = -1;
    double s = 0;
    do {
        s += ExponentialRand::StandardVariate();
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

std::complex<double> PoissonRand::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t) - 1.0, std::sin(t));
    return std::exp(lambda * y);
}

double PoissonRand::Median() const
{
    /// this value is approximate
    return std::floor(lambda + 1.0 / 3 - 0.02 / lambda);
}

int PoissonRand::Mode() const
{
    return floorLambda;
}

double PoissonRand::Skewness() const
{
    return 1.0 / std::sqrt(lambda);
}

double PoissonRand::ExcessKurtosis() const
{
    return 1.0 / lambda;
}

bool PoissonRand::FitMLE(const std::vector<int> &sample)
{
    if (!checkValidity(sample))
        return false;
    SetRate(sampleMean(sample));
    return true;
}

bool PoissonRand::FitMM(const std::vector<int> &sample)
{
    return FitMLE(sample);
}

bool PoissonRand::FitBayes(const std::vector<int> &sample, GammaRand &priorDistribution)
{
    if (!checkValidity(sample))
        return false;
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    priorDistribution.SetParameters(alpha + sampleSum(sample), beta + sample.size());
    SetRate(priorDistribution.Mean());
    return true;
}
