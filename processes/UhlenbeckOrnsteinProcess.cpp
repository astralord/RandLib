#include "UhlenbeckOrnsteinProcess.h"

UhlenbeckOrnsteinProcess::UhlenbeckOrnsteinProcess(double deltaT, double mean, double volatility, double rate, double initialValue) :
    StochasticProcess(deltaT)
{
    setParameters(mean, volatility, rate, initialValue);
    lastRateTime = 0.0;
    lastExpRateTime = 1.0;
}

void UhlenbeckOrnsteinProcess::setParameters(double mean, double volatility, double rate, double initialValue)
{
    mu = mean;
    sigma = volatility;
    if (sigma <= 0)
        sigma = 1.0;
    theta = rate;
    if (theta <= 0)
        theta = 1.0;
    x0 = initialValue;

    W.setSigma(sigma / std::sqrt(theta + theta));
}

double UhlenbeckOrnsteinProcess::next() const
{
    return next(dt);
}

double UhlenbeckOrnsteinProcess::next(double deltaT) const
{
    // TODO: redo with faster implementation - use directly NormalRand instead of WienerProcess
    double thetaT = lastRateTime + theta * deltaT;
    double expThetaT = std::exp(thetaT);

    double wienerDeltaT = expThetaT * expThetaT - lastExpRateTime * lastExpRateTime;
    double y = W.next(wienerDeltaT);
    y += x0 - mu;
    y /= expThetaT;
    y += mu;

    lastRateTime = thetaT;
    lastExpRateTime = expThetaT;

    return y;
}

void UhlenbeckOrnsteinProcess::Mean(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
    {
        output[i] = std::exp(-theta * time[i]);
        output[i] *= x0 - mu;
        output[i] += mu;
    }
}

void UhlenbeckOrnsteinProcess::Variance(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    double coef = 0.5 * sigma * sigma / theta;
    for (int i = 0; i < size; ++i)
        output[i] = coef * (1.0 - std::exp(-2.0 * theta * time[i]));
}

