#include "RayleighRand.h"

RayleighRand::RayleighRand(double scale)
{
    setScale(scale);
}

std::string RayleighRand::name()
{
    return "Rayleigh(" + toStringWithPrecision(getScale()) + ")";
}

void RayleighRand::setScale(double scale)
{
    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    sigmaSqInv = 1.0 / (sigma * sigma);
}

double RayleighRand::f(double x) const
{
    if (x <= 0)
        return 0.0;
    double y = x * sigmaSqInv;
    return y * std::exp(-0.5 * x * y);
}

double RayleighRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    return 1.0 - std::exp(-0.5 * x * x * sigmaSqInv);
}

double RayleighRand::variate() const
{
    double rv = ExponentialRand::standardVariate();
    return sigma * std::sqrt(rv + rv);
}

double RayleighRand::Mean() const
{
    return sigma * M_SQRTPI * M_SQRT1_2;
}

double RayleighRand::Variance() const
{
    return sigma * sigma * (2.0 - M_PI_2);
}

double RayleighRand::Quantile(double p) const
{
    if (p < 0.0 || p > 1.0)
        return NAN;
    if (p == 1.0)
        return INFINITY;
    return sigma * std::sqrt(-2 * std::log(1 - p));
}

double RayleighRand::Median() const
{
    return sigma * std::sqrt(M_LN2 + M_LN2);
}

double RayleighRand::Mode() const
{
    return sigma;
}

double RayleighRand::Skewness() const
{
    return 2 * M_SQRTPI * (M_PI - 3) / std::pow(4.0 - M_PI, 1.5);
}

double RayleighRand::ExcessKurtosis() const
{
    return (6 * M_PI - 16.0 / (M_PI - 4)) / (M_PI - 4);
}

bool RayleighRand::fitToData(const std::vector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate sigma^2
    double sigmaEst = 0.0;
    int N = sample.size();
    for (double var : sample)
    {
        if (var < 0)
            return false;
        sigmaEst += var * var;
    }
    sigmaEst *= .5 / N;

    /// Calculate unbiased sigma
    sigmaEst = std::sqrt(sigmaEst);

    if (N > 30)
        setScale((1 + 0.1252 / N) * sigmaEst); /// err < 1e-6
    else
    {
        double coef = RandMath::factorial(N - 1);
        coef *= N * coef;
        coef *= M_1_SQRTPI * std::sqrt(static_cast<double>(N));
        coef /= RandMath::factorial(N << 1);
        int pow2N = 1 << N; /// < 2^31
        coef *= pow2N;
        coef *= pow2N;

        setScale(coef * sigmaEst);
    }
    return true;
}
