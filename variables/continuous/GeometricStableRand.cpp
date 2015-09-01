#include "GeometricStableRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location) :
    S(exponent, skewness)
{
    setParameters(exponent, skewness, scale, location);
}

std::string GeometricStableRand::name()
{
    return "Geometric Stable("
            + toStringWithPrecision(getAlpha()) + ", "
            + toStringWithPrecision(getBeta()) + ", "
            + toStringWithPrecision(getSigma()) + ", "
            + toStringWithPrecision(getMu()) + ")";
}

void GeometricStableRand::setParameters(double exponent, double skewness, double scale, double location)
{
    sigma = std::max(scale, MIN_POSITIVE);
    mu = location;

    if (std::fabs(exponent - 1) < MIN_POSITIVE)
        S.setParameters(1.0, skewness, sigma, mu);
    else
        S.setParameters(exponent, skewness, 1, 0);
}

double GeometricStableRand::f(double x) const
{
    return x;
}

double GeometricStableRand::F(double x) const
{
    return x;
}

double GeometricStableRand::variate() const
{
    double e = ExponentialRand::standardVariate();
    double x = S.variate();
    double alphaInv = 1.0 / S.getAlpha();
    if (alphaInv == 1)
    {
        double rv = x + M_2_PI * S.getBeta() * std::log(sigma * e);
        return e * (mu + sigma * rv);
    }
    return mu * e + std::pow(e, alphaInv) * sigma * x;
}

void GeometricStableRand::sample(QVector<double> &outputData)
{
    S.sample(outputData);
    double alphaInv = 1.0 / S.getAlpha();
    if (alphaInv == 1) {
        double beta = S.getBeta();
        for (double &var : outputData) {
            double e = ExponentialRand::standardVariate();
            var += M_2_PI * beta * sigma * std::log(sigma * e);
            var *= e;
        }
    }
    else {
        for (double &var : outputData) {
            double e = ExponentialRand::standardVariate();
            var = mu * e + std::pow(e, alphaInv) * sigma * var;
        }
    }
}

