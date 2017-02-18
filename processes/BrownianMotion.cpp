#include "BrownianMotion.h"

BrownianMotion::BrownianMotion(double deltaT) :
    StableProcess(2.0, 0.0, 0.0, M_SQRT1_2, deltaT)
{

}

BrownianMotion::BrownianMotion(double drift, double volatility, double deltaT) :
    StableProcess(2.0, 0.0, drift, volatility * M_SQRT1_2, deltaT)
{
}

void BrownianMotion::NextImpl()
{
    currentValue += mu * dt + sigma * dtCoef * NormalRand::StandardVariate();
}

void BrownianMotion::ProbabilityDensityFunction(double t, const std::vector<double> &x, std::vector<double> &y) const
{
    double tau = t - currentTime;
    if (tau <= 0)
        return;
    NormalRand X(currentValue + mu * tau, sigma * sigma * tau);
    return X.ProbabilityDensityFunction(x, y);
}

void BrownianMotion::CumulativeDistributionFunction(double t, const std::vector<double> &x, std::vector<double> &y) const
{
    double tau = t - currentTime;
    if (tau <= 0)
        return;
    NormalRand X(currentValue + mu * tau, sigma * sigma * tau);
    return X.CumulativeDistributionFunction(x, y);
}

double BrownianMotion::Quantile(double t, double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    double tau = t - currentTime;
    double mean = currentValue + mu * tau;
    double scale = sigma * std::sqrt(tau);
    return mean - scale * M_SQRT2 * RandMath::erfcinv(2 * p);
}
