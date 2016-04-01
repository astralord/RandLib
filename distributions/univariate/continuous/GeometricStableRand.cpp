#include "GeometricStableRand.h"
#include "ExponentialRand.h"
#include "LaplaceRand.h"
#include "UniformRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location) :
    StableRand(exponent, skewness, scale, location)
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

double GeometricStableRand::parameterTransform(double x, double z)
{
    double par = x - mu * z;
    par /= std::pow(z, alphaInv);
    return par + mu;
}

double GeometricStableRand::f(double x) const
{
    // TODO: find mode and separate integrals according to this peak
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 0.0;
        double par = parameterTransform(x, z);
        double coef = log(z) * alphaInv + z;
        return std::exp(-coef) * StableRand::f(par);
    },
    0, 10); // ideal maximum border is infinity, but fortunately integrand decreases very fast
}

double GeometricStableRand::F(double x) const
{
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 0.0;
        double par = parameterTransform(x, z);
        return std::exp(-z) * StableRand::F(par);
    },
    0, 10); // ideal maximum border is infinity, but fortunately integrand decreases very fast
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

double GeometricStableRand::variateForAlphaEqualTwo() const
{
    double Z = ExponentialRand::standardVariate();
    double X = NormalRand::standardVariate();
    return mu * Z + std::sqrt(Z) * sigma * X;
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

double GeometricStableRand::ExcessKurtosis() const
{
    return (alpha == 2) ? 3.0 : NAN;
}
