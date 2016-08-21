#include "NakagamiRand.h"

NakagamiRand::NakagamiRand(double shape, double spread)
{
    SetParameters(shape, spread);
}

std::string NakagamiRand::Name() const
{
    return "Nakagami(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetSpread()) + ")";
}

void NakagamiRand::SetParameters(double shape, double spread)
{
    m = std::max(shape, 0.5);
    w = spread;
    if (w <= 0)
        w = 1.0;

    sigma = m / w;
    Y.SetParameters(m, sigma);
}

double NakagamiRand::f(double x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0)
        return (m > 0.5) ? 0.0 : M_SQRT2 / M_SQRTPI / std::sqrt(w);
    return (x < 0) ? 0.0 : 2 * x * Y.f(x * x);
}

double NakagamiRand::F(double x) const
{
    return (x <= 0) ? 0.0 : Y.F(x * x);
}

double NakagamiRand::Variate() const
{
    return std::sqrt(Y.Variate());
}

void NakagamiRand::Sample(std::vector<double> &outputData) const
{
    Y.Sample(outputData);
    for (double & var : outputData)
        var = std::sqrt(var);
}

double NakagamiRand::Mean() const
{
    double y = std::lgamma(m + 0.5);
    y -= Y.GetLogGammaFunction();
    y -= 0.5 * std::log(sigma);
    return std::exp(y);
}

double NakagamiRand::Variance() const
{
    double y = std::lgamma(m + 0.5);
    y -= Y.GetLogGammaFunction();
    y = std::exp(y + y);
    return w * (1 - y / m);
}

double NakagamiRand::Mode() const
{
    double mode = 0.5 * w / m;
    return std::sqrt(w - mode);
}


/// CHI
ChiRand::ChiRand(int degree, double scale)
{
    SetParameters(degree, scale);
}

std::string ChiRand::Name() const
{
    return "Chi(" + toStringWithPrecision(GetDegree()) + ", " + toStringWithPrecision(GetScale()) +  ")";
}

void ChiRand::SetParameters(int degree, double scale)
{
    sigma = scale;
    if (sigma <= 0)
        sigma = 1.0;
    sigmaSqInv = 1.0 / (sigma * sigma);

    NakagamiRand::SetParameters(0.5 * degree, degree * sigma * sigma);
}

double ChiRand::skewnessImpl(double mean, double sigma) const
{
    double variance = sigma * sigma;
    double y = mean * (1 - 2 * variance);
    return y / (sigma * variance);
}

double ChiRand::Skewness() const
{
    double mean = Mean();
    return skewnessImpl(mean, std::sqrt(Variance()));
}

double ChiRand::ExcessKurtosis() const
{
    double mean = Mean();
    double variance = Variance();
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

std::string MaxwellBoltzmannRand::Name() const
{
    return "Maxwell-Boltzmann(" + toStringWithPrecision(GetScale()) + ")";
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


/// RAYLEIGH
RayleighRand::RayleighRand(double scale) :
    ChiRand(2, scale)
{
}

std::string RayleighRand::Name() const
{
    return "Rayleigh(" + toStringWithPrecision(GetScale()) + ")";
}

void RayleighRand::SetScale(double scale)
{
    ChiRand::SetParameters(2, scale);
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

double RayleighRand::quantileImpl(double p) const
{
    return sigma * std::sqrt(-2 * std::log1p(-p));
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

bool RayleighRand::FitScaleMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    double sigmaSq = 0.5 * rawMoment(sample, 2);
    SetScale(std::sqrt(sigmaSq));
    return true;
}

bool RayleighRand::FitScaleUMVU(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    size_t n = sample.size();
    if (n == 0)
        return false;

    double sigmaSq = 0.5 * rawMoment(sample, 2);

    /// Calculate unbiased sigma
    double sigmaEst = std::sqrt(sigmaSq);

    if (n > 30)
        SetScale((1 + 0.1252 / n) * sigmaEst); /// err < 1e-6
    else
    {
        double coef = RandMath::factorial(n - 1);
        coef *= n * coef;
        coef *= M_1_SQRTPI * std::sqrt(static_cast<double>(n));
        coef /= RandMath::factorial(n << 1);
        int pow2n = 1 << n; /// < 2^31
        coef *= pow2n;
        coef *= pow2n;

        SetScale(coef * sigmaEst);
    }
    return true;
}

