#include "GeometricBrownianMotion.h"

GeometricBrownianMotion::GeometricBrownianMotion(double deltaT, double drift, double volatility, double initialValue) :
    StochasticProcess(deltaT),
    mu(drift),
    sigma((volatility <= 0) ? 1.0 : volatility),
    S0(initialValue),
    mumSigma2_2(mu - 0.5 * sigma * sigma),
    B(deltaT)
{
    currentValue = S0;
}

void GeometricBrownianMotion::nextImpl()
{
    currentValue = B.next();
    currentValue *= sigma;
    currentValue += mumSigma2_2 * currentTime;
    currentValue = S0 * std::exp(currentValue);
}

void GeometricBrownianMotion::nextImpl(double deltaT)
{
    currentValue = B.next(deltaT);
    currentValue *= sigma;
    currentValue += mumSigma2_2 * currentTime;
    currentValue = S0 * std::exp(currentValue);
}

double GeometricBrownianMotion::MeanImpl(double t) const
{
    return currentValue * std::exp(mu * (t - currentTime));
}

double GeometricBrownianMotion::VarianceImpl(double t) const
{
    double var = std::exp(sigma * sigma * (t - currentTime));
    --var;
    var *= std::exp(2 * mu * t);
    return var * currentValue * currentValue;
}

double GeometricBrownianMotion::QuantileImpl(double t, double p) const
{
    double tau = t - currentTime;
    return currentValue * LogNormalRand::quantile(p, mumSigma2_2 * tau, sigma * std::sqrt(tau));
}
