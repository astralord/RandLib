#include "CoxIngersollRossProcess.h"

CoxIngersollRossProcess::CoxIngersollRossProcess(double mean, double reversionSpeed, double volatility, double initialValue, double deltaT) :
    StochasticProcess(deltaT, initialValue),
    alpha(mean > 0 ? mean : 1.0),
    beta(reversionSpeed > 0 ? reversionSpeed : 1.0),
    sigma(volatility > 0 ? volatility : 1.0),
    expmBetaDt(std::exp(-beta * dt)),
    c(4 * beta / (sigma * sigma * (1.0 - expmBetaDt))),
    degree(4 * alpha * beta / (sigma * sigma))
{

}

void CoxIngersollRossProcess::NextImpl()
{
    currentValue = NoncentralChiSquaredRand::Variate(degree, c * currentValue * expmBetaDt) / c;
}

double CoxIngersollRossProcess::MeanImpl(double t) const
{
    double expmBetaDeltaT = std::exp(t - currentTime);
    return currentValue * expmBetaDeltaT + alpha * (1.0 - expmBetaDeltaT);
}

double CoxIngersollRossProcess::VarianceImpl(double t) const
{
    double expmBetaDeltaT = std::exp(t - currentTime);
    double aux = 1.0 - expmBetaDeltaT;
    double var = currentValue * expmBetaDeltaT;
    var += 0.5 * alpha * aux;
    var *= aux;
    var *= sigma * sigma / beta;
    return var;
}
