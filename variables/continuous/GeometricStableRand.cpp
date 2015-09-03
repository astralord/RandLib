#include "GeometricStableRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location) :
    StableRand(exponent, skewness, scale, location)
{
}

std::string GeometricStableRand::name()
{
    return "Geometric Stable("
            + toStringWithPrecision(getAlpha()) + ", "
            + toStringWithPrecision(getBeta()) + ", "
            + toStringWithPrecision(getSigma()) + ", "
            + toStringWithPrecision(getMu()) + ")";
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
    if (alpha == 2 && mu == 0 && beta == 0)
        return LaplaceRand::variate(0, sigma);
    double W = ExponentialRand::standardVariate();
    double Y = StableRand::variate();
    if (alphaInv == 1)
        return W * (Y + M_2_PI * beta * sigma * std::log(sigma * W));
    double W_adj = std::pow(W, alphaInv);
    return W_adj * Y + mu * (W - W_adj);
}

void GeometricStableRand::sample(QVector<double> &outputData)
{
    if (alpha == 2 && mu == 0 && beta == 0) {
        for (double &var : outputData) {
            var = LaplaceRand::variate(0, sigma);
        }
    }
    else {
        StableRand::sample(outputData);
        if (alphaInv == 1) {
            for (double &var : outputData) {
                double W = ExponentialRand::standardVariate();
                var += M_2_PI * beta * sigma * std::log(sigma * W);
                var *= W;
            }
        }
        else {
            for (double &var : outputData) {
                double W = ExponentialRand::standardVariate();
                double W_adj = std::pow(W, alphaInv);
                var *= W_adj;
                var += mu * (W - W_adj);
            }
        }
    }
}

