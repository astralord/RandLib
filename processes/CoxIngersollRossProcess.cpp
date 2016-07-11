#include "CoxIngersollRossProcess.h"

CoxIngersollRossProcess::CoxIngersollRossProcess(double drift, double reversionSpeed, double volatility, double initialValue, double deltaT) :
    StochasticProcess(deltaT, initialValue),
    alpha(drift > 0 ? drift : 1.0),
    beta(reversionSpeed > 0 ? reversionSpeed : 1.0),
    sigma(volatility > 0 ? volatility : 1.0),
    expmBetaDt(std::exp(-beta * dt)),
    c(4 * beta / (sigma * sigma * (1.0 - expmBetaDt))),
    degree(4 * alpha / (sigma * sigma))
{

}

void CoxIngersollRossProcess::nextImpl()
{
    currentValue = NoncentralChiSquared::variate(degree, c * currentValue * expmBetaDt) / c;
}

double CoxIngersollRossProcess::MeanImpl(double t) const
{
    double expmBetaDeltaT = std::exp(t - currentTime);
    return currentValue * expmBetaDeltaT + alpha / beta * (1.0 - expmBetaDeltaT);
}

double CoxIngersollRossProcess::VarianceImpl(double t) const
{
    double expmBetaDeltaT = std::exp(t - currentTime);
    double aux = 1.0 - expmBetaDeltaT;
    double var = currentValue * expmBetaDeltaT;
    var += 0.5 * alpha / beta * aux;
    var *= aux;
    var *= sigma * sigma / beta;
    return var;
}

double CoxIngersollRossProcess::QuantileImpl(double t, double p) const
{
    // very rough approximation by normal distribution
    // can work only for big t and existing mean and variance
    NormalRand X(Mean(t), Variance(t));
    return X.Quantile(p);
}
