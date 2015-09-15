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

double GeometricStableRand::variateForAlphaEqualOne() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W1 = ExponentialRand::standardVariate();
    double W2 = ExponentialRand::standardVariate();
    double pi_2BetaU = M_PI_2 + beta * U;
    double X = logSigma;
    X += std::log(sigma * W2 * pi_2BetaU / (M_PI_2 * W1 * std::cos(U)));
    X *= beta;
    X += pi_2BetaU * std::tan(U);
    X *= M_2_PI * sigma;
    X += mu;
    X *= W2;
    return X;
}

double GeometricStableRand::variateForCommonAlpha() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::standardVariate();
    double alphaUB = alpha * U + B;
    double X = S * std::sin(alphaUB);
    double W_adj = W / std::cos(U - alphaUB);
    X *= W_adj;
    W = ExponentialRand::standardVariate();
    W_adj *= std::cos(U);
    W_adj = W / W_adj;
    X *= std::pow(W_adj, alphaInv);
    return X * sigma + mu * W;
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
    else if (alphaInv == 1) {
        for (double &var : outputData) {
            var = variateForAlphaEqualOne();
    }
    else {
        for (double &var : outputData) {
            var = variateForCommonAlpha();
        }
    }
}

std::complex<double> GeometricStableRand::CF(double t) const
{
    return 1.0 / (1.0 + psi(t));
}

double GeometricStableRand::ExcessKurtosis() const
{
    return (alpha == 2) ? 3.0 : NAN;
}
