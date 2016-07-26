#include "OrnsteinUhlenbeckProcess.h"

OrnsteinUhlenbeckProcess::OrnsteinUhlenbeckProcess(double mean, double reversionSpeed, double volatility, double initialValue, double deltaT) :
    StochasticProcess(deltaT, initialValue),
    alpha(mean > 0 ? mean : 1.0),
    beta(reversionSpeed > 0 ? reversionSpeed : 1.0),
    sigma(volatility > 0 ? volatility : 1.0),
    expmBetaDt(std::exp(-beta * dt)),
    X(alpha * (1 - expmBetaDt), 0.5 * sigma * sigma / beta * (1 - expmBetaDt * expmBetaDt))
{

}

void OrnsteinUhlenbeckProcess::nextImpl()
{
    currentValue *= expmBetaDt;
    currentValue += X.variate();
}

double OrnsteinUhlenbeckProcess::MeanImpl(double t) const
{
    double expmBetaDeltaT = std::exp(-beta * (t - currentTime));
    return expmBetaDeltaT * currentValue + alpha * (1.0 - expmBetaDeltaT);
}

double OrnsteinUhlenbeckProcess::VarianceImpl(double t) const
{
    double expmBetaDeltaT = std::exp(-beta * (t - currentTime));
    return 0.5 * sigma * sigma / beta * (1 - expmBetaDeltaT * expmBetaDeltaT);
}

double OrnsteinUhlenbeckProcess::Quantile(double t, double p) const
{
    double expmBetaDeltaT = std::exp(-beta * (t - currentTime));
    return NormalRand::quantile(p, expmBetaDeltaT * currentValue + alpha * (1 - expmBetaDeltaT), sigma * std::sqrt(0.5 / beta * (1 - expmBetaDeltaT * expmBetaDeltaT)));
}

