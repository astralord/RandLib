#include "GeometricBrownianMotion.h"

GeometricBrownianMotion::GeometricBrownianMotion(double deltaT, double drift, double volatility, double initialValue) :
    StochasticProcess(deltaT)
{
    setParameters(drift, volatility, initialValue);
}

void GeometricBrownianMotion::setParameters(double drift, double volatility, double initialValue)
{
    mu = drift;
    sigma = volatility;
    if (sigma <= 0)
        sigma = 1.0;
    S0 = initialValue;
    W.setMean(mu - 0.5 * sigma * sigma);
    W.setSigma(sigma);
}

double GeometricBrownianMotion::next() const
{
    return S0 * std::exp(W.next());
}

double GeometricBrownianMotion::next(double deltaT) const
{
    return S0 * std::exp(W.next(deltaT));
}

void GeometricBrownianMotion::Mean(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    for (int i = 0; i < size; ++i)
        output[i] = S0 * std::exp(mu * time[i]);
}

void GeometricBrownianMotion::Variance(const QVector<double> &time, QVector<double> &output) const
{
    int size = std::min(time.size(), output.size());
    double var = sigma * sigma;
    for (int i = 0; i < size; ++i)
    {
        output[i] = S0 * std::exp(mu * time[i]);
        output[i] *= output[i];
        output[i] *= std::expm1(var * time[i]);
    }
}
