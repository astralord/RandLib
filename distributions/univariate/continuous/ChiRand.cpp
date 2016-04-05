#include "ChiRand.h"

ChiRand::ChiRand(int degree, double scale)
{
    setDegree(degree);
    setScale(scale);
}

std::string ChiRand::name()
{
    return "Chi(" + toStringWithPrecision(getDegree()) + ", " + toStringWithPrecision(getScale()) +  ")";
}

void ChiRand::setDegree(int degree)
{
    X.setDegree(degree);
}

void ChiRand::setScale(double scale)
{
    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    sigmaSqInv = 1.0 / (sigma * sigma);
}

double ChiRand::f(double x) const
{
    if (x <= 0)
        return 0.0;
    return 2 * x * X.f(x * x * sigmaSqInv) * sigmaSqInv;
}

double ChiRand::F(double x) const
{
    if (x <= 0)
        return 0.0;
    return X.F(x * x * sigmaSqInv);
}

double ChiRand::variate() const
{
    return sigma * std::sqrt(X.variate());
}

void ChiRand::sample(std::vector<double> &outputData) const
{
    X.sample(outputData);
    for (double & var : outputData)
        var = sigma * std::sqrt(var);
}

double ChiRand::varianceImpl(double mean) const
{
    return (v - mean * mean) / sigmaSqInv;
}

double ChiRand::skewnessImpl(double mean, double sigma) const
{
    double variance = sigma * sigma;
    double y = mean * (1 - 2 * variance);
    return y / (sigma * variance);
}

double ChiRand::Mean() const
{
    double halfV = 0.5 * X.getDegree();
    double x = std::lgamma(halfV + 0.5);
    double y = std::lgamma(halfV);
    return sigma * M_SQRT2 * std::exp(x - y);
}

double ChiRand::Variance() const
{
    return varianceImpl(Mean());
}

double ChiRand::Mode() const
{
    double v = X.getDegree();
    return (v > 1) ? sigma * std::sqrt(v - 1) : 0;
}

double ChiRand::Skewness() const
{
    double mean = Mean();
    return skewnessImpl(mean, std::sqrt(varianceImpl(mean)));
}

double ChiRand::ExcessKurtosis() const
{
    double mean = Mean();
    double variance = varianceImpl(mean);
    double sigma = std::sqrt(variance);
    double skewness = skewnessImpl(mean, sigma);
    double y = 1 - mean * sigma * skewness;
    y /= variance;
    --y;
    return y + y;
}

/// MAXWELL-BOLTZMANN
MaxwellBoltzmannRand::MaxwellBoltzmannRand(double scale) :
    ChiRand(3, scale)
{
}

std::string MaxwellBoltzmannRand::name()
{
    return "Maxwell-Boltzmann(" + toStringWithPrecision(getScale()) + ")";
}

double MaxwellBoltzmannRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double x2 = x * x;
    double y = std::exp(-.5 * x2 * sigmaSqInv);
    return M_SQRT2 * M_1_SQRTPI * x2 * y * sigmaSqInv / sigma;
}

double MaxwellBoltzmannRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double xAdj = M_SQRT1_2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= M_SQRT2 * M_1_SQRTPI * x / sigma;
    return std::erf(xAdj) - y;
}

double MaxwellBoltzmannRand::Mean() const
{
    return 2 * M_1_SQRTPI * M_SQRT2 * sigma;
}

double MaxwellBoltzmannRand::Variance() const
{
    return (3 - 8.0 * M_1_PI) / sigmaSqInv;
}

double MaxwellBoltzmannRand::Mode() const
{
    return M_SQRT2 * sigma;
}

double MaxwellBoltzmannRand::Skewness() const
{
    double skewness = 3 * M_PI - 8;
    skewness = 2.0 / skewness;
    skewness *= std::sqrt(skewness);
    return (16 - 5 * M_PI) * skewness;
}

double MaxwellBoltzmannRand::ExcessKurtosis() const
{
    double numerator = 40 - 3 * M_PI;
    numerator *= M_PI;
    numerator -= 96;
    double denominator = 3 * M_PI - 8;
    denominator *= denominator;
    return 4 * numerator / denominator;
}


RayleighRand::RayleighRand(double scale) :
    ChiRand(2, scale)
{
}

std::string RayleighRand::name()
{
    return "Rayleigh(" + toStringWithPrecision(getScale()) + ")";
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

double RayleighRand::Mean() const
{
    return sigma * M_SQRTPI * M_SQRT1_2;
}

double RayleighRand::Variance() const
{
    return (2.0 - M_PI_2) / sigmaSqInv;
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
