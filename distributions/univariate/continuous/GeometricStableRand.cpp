#include "GeometricStableRand.h"
#include "ExponentialRand.h"
#include "LaplaceRand.h"
#include "UniformRand.h"

GeometricStableRand::GeometricStableRand(double exponent, double skewness, double scale, double location)
    : LimitingDistribution(exponent, skewness, scale, location),
      Z(exponent, skewness)
{
    setParameters(exponent, skewness);
    setScale(scale);
    setLocation(location);
}

std::string GeometricStableRand::name() const
{
    return "Geometric Stable("
            + toStringWithPrecision(getExponent()) + ", "
            + toStringWithPrecision(getSkewness()) + ", "
            + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getLocation()) + ")";
}

void GeometricStableRand::setAsymmetry()
{
    /// Calculate value of asymmetry
    if (mu == 0)
        k = 1;
    else {
        double sigma2 = sigma + sigma;
        k = mu * mu + sigma2 * sigma2;
        k = std::sqrt(k);
        k -= mu;
        k /= sigma2;
    }
    kInv = 1.0 / k;
    kSq = k * k;
    pdfCoef = 1.0 / (sigma * (k + kInv));
    cdfCoef = 1.0 / (1 + kSq);
}

void GeometricStableRand::setParameters(double exponent, double skewness)
{
    LimitingDistribution::setParameters(exponent, skewness);
    Z.setParameters(alpha, beta);

    if (alpha == 2)
        setAsymmetry();
}

void GeometricStableRand::setLocation(double location)
{
    LimitingDistribution::setLocation(location);
    if (alpha == 2)
        setAsymmetry();
}

void GeometricStableRand::setScale(double scale)
{
    LimitingDistribution::setScale(scale);
    pdfCoef = 1.0 / (sigma * (k + kInv));
    if (alpha == 2)
        setAsymmetry();
}

double GeometricStableRand::pdfLaplace(double x) const
{
    double y = x / sigma;
    y *= (x < 0) ? kInv : -k;
    y = std::exp(y);
    return pdfCoef * y;
}

double GeometricStableRand::cdfLaplace(double x) const
{
    double y = x / sigma;
    if (x < 0) {
        y *= kInv;
        y = std::exp(y);
        return kSq * cdfCoef * y;
    }
    y *= -k;
    y = std::exp(y);
    return 1.0 - cdfCoef * y;
}

double GeometricStableRand::f(double x) const
{
    if (alpha == 2)
        return pdfLaplace(x);
    // TODO: find mode and separate integrals according to this peak
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 0.0;
        double temp = sigma * std::pow(z, alphaInv);
        double par = (x - mu * z) / temp;
        double y = std::exp(-z) / temp * Z.f(par);
        return y;
    },
    0, 10); // ideal maximum border is infinity, but fortunately integrand decreases very fast
}

double GeometricStableRand::F(double x) const
{
    if (alpha == 2)
        return cdfLaplace(x);
    return RandMath::integral([this, x] (double z)
    {
        if (z == 0)
            return 1.0;
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
    double Y = ExponentialRand::standardVariate();
    W_adj *= std::cos(U);
    W_adj = Y / W_adj;
    X *= std::pow(W_adj, alphaInv);
    return X * sigma + mu * Y;
}

double GeometricStableRand::variate() const
{
    if (alpha == 2)
        return (mu == 0) ? LaplaceRand::variate(0, sigma) : LaplaceRand::variate(0, sigma, k);
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
                var = LaplaceRand::variate(0, sigma, k);
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

double GeometricStableRand::Variance() const
{
    if (alpha == 2) {
        double y = sigma * sigma / kSq;
        return (1.0 + kSq * kSq) * y;
    }
    return INFINITY;
}

std::complex<double> GeometricStableRand::CF(double t) const
{
    return (t == 0) ? 1.0 : 1.0 / (1.0 + psi(t));
}

double GeometricStableRand::Median() const
{
    if (alpha == 2) {
        double y = 0.5 * (1.0 / kSq + 1.0);
        y = std::log(y);
        return sigma * k * y;
    }
    return ContinuousDistribution::Median();
}

double GeometricStableRand::Mode() const
{
    return (alpha == 2) ? 0.0 : ContinuousDistribution::Mode();
}

double GeometricStableRand::Skewness() const
{
    if (alpha == 2) {
        double kSq = k * k;
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
    if (alpha == 2)
    {
        if (k == 1)
            return 3.0;
        return NAN; // TODO!
    }
    return NAN;
}
