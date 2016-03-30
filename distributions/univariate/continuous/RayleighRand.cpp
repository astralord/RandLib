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

bool RayleighRand::checkValidity(const std::vector<double> &sample)
{
    for (double var : sample) {
        if (var < 0)
            return false;
    }
    return true;
}

bool RayleighRand::fitScaleMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double sigmaSq = 0.5 * RandMath::rawMoment(sample, 2);
    setScale(std::sqrt(sigmaSq));
    return true;
}

bool RayleighRand::fitScaleUMVU(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    size_t n = sample.size();
    if (n == 0)
        return false;

    double sigmaSq = 0.5 * RandMath::rawMoment(sample, 2);

    /// Calculate unbiased sigma
    double sigmaEst = std::sqrt(sigmaSq);

    if (n > 30)
        setScale((1 + 0.1252 / n) * sigmaEst); /// err < 1e-6
    else
    {
        double coef = RandMath::factorial(n - 1);
        coef *= n * coef;
        coef *= M_1_SQRTPI * std::sqrt(static_cast<double>(n));
        coef /= RandMath::factorial(n << 1);
        int pow2n = 1 << n; /// < 2^31
        coef *= pow2n;
        coef *= pow2n;

        setScale(coef * sigmaEst);
    }
    return true;
}
