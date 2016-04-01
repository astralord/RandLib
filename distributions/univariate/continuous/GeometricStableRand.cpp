#include "GeometricStableRand.h"
#include "ExponentialRand.h"
#include "LaplaceRand.h"
#include "UniformRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location),
      Z(exponent, skewness)
{
}

std::string GeometricStableRand::name()
{
    return "Geometric Stable("
            + toStringWithPrecision(getExponent()) + ", "
            + toStringWithPrecision(getSkewness()) + ", "
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}

void GeometricStableRand::setParameters(double exponent, double skewness)
{
    LimitingDistribution::setParameters(exponent, skewness);
    Z.setParameters(alpha, beta);
}

double GeometricStableRand::f(double x) const
{
    // TODO: find mode and separate integrals according to this peak
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 0.0;
        double temp = sigma * std::pow(z, alphaInv);
        double par = (x - mu * z) / temp;
        return std::exp(-z) / temp * Z.f(par);
    },
    0, 10); // ideal maximum border is infinity, but fortunately integrand decreases very fast
}

double GeometricStableRand::F(double x) const
{
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 0.0;
        double par = x - mu * z;
        par /= std::pow(z, alphaInv) * sigma;
        return std::exp(-z) * Z.F(par);
    },
    0, 10); // ideal maximum border is infinity, but fortunately integrand decreases very fast
}

double GeometricStableRand::variateForAlphaEqualOne() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W1 = ExponentialRand::standardVariate();
    double W2 = ExponentialRand::standardVariate();
    double pi_2BetaU = M_PI_2 + beta * U;
    double X = 2 * logSigma;
    X += std::log(W2 * pi_2BetaU / (M_PI_2 * W1 * std::cos(U)));
    X *= beta;
    X += pi_2BetaU * std::tan(U);
    X *= M_2_PI * sigma;
    X += mu;
    X *= W2;
    return X;
}

double GeometricStableRand::variateForAlphaEqualTwo() const
{
    double W = ExponentialRand::standardVariate();
    double X = NormalRand::standardVariate();
    return mu * W + std::sqrt(W) * sigma * X;
}

double GeometricStableRand::variateForCommonAlpha() const
{
    double U = UniformRand::variate(-M_PI_2, M_PI_2);
    double W = ExponentialRand::standardVariate();
    double alphaUB = alpha * U + B;
    double X = S * std::sin(alphaUB);
    double W_adj = W / std::cos(U - alphaUB);
    X *= W_adj;
    double Y = ExponentialRand::standardVariate();
    W_adj *= std::cos(U);
    W_adj = Y / W_adj;
    X *= std::pow(W_adj, alphaInv);
    return X * sigma + mu * Y;
}

double GeometricStableRand::variate() const
{
    if (alpha == 2)
    {
        if (mu == 0)
            return LaplaceRand::variate(0, sigma);
        return variateForAlphaEqualTwo();
    }
    return (alpha == 1) ? variateForAlphaEqualOne() : variateForCommonAlpha();
}

void GeometricStableRand::sample(std::vector<double> &outputData) const
{
    if (alpha == 2) {
        if (mu == 0) {
            for (double &var : outputData)
                var = LaplaceRand::variate(0, sigma);
        }
        else {
            for (double &var : outputData)
                var = variateForAlphaEqualTwo();
        }
    }
    else if (alpha == 1) {
        for (double &var : outputData)
            var = variateForAlphaEqualOne();
    }
    else {
        for (double &var : outputData)
            var = variateForCommonAlpha();
    }
}

std::complex<double> GeometricStableRand::CF(double t) const
{
    return 1.0 / (1.0 + psi(t));
}

double GeometricStableRand::Skewness() const
{
    if (alpha == 2) {
        double kSq = 1.0; //TODO!
        double k4 = kSq * kSq;
        double k6 = k4 * kSq;
        double z = (k4 + 1);
        z *= std::sqrt(z);
        double y = 2 * (1 - k6);
        return y / z;
    }
    return NAN;
}

double GeometricStableRand::ExcessKurtosis() const
{
    return (alpha == 2) ? 3.0 : NAN;
}
