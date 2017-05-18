#include "NakagamiRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

NakagamiDistribution::NakagamiDistribution(double shape, double spread)
{
    SetParameters(shape, spread);
}

void NakagamiDistribution::SetParameters(double shape, double spread)
{
    m = (shape > 0.0) ? shape : 1.0;
    w = (spread > 0.0) ? spread : 1.0;
    Y.SetParameters(m, m / w);
    lgammaShapeRatio = std::lgamma(m + 0.5) - Y.GetLogGammaShape();
}

double NakagamiDistribution::f(const double & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0) {
        if (m > 0.5)
            return 0.0;
        return (m < 0.5) ? INFINITY : std::sqrt(M_2_PI / w);
    }
    return 2 * x * Y.f(x * x);
}

double NakagamiDistribution::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0) {
        if (m > 0.5)
            return 0-INFINITY;
        return (m < 0.5) ? INFINITY : 0.5 * std::log(M_2_PI / w);
    }
    return std::log(2 * x) + Y.logf(x * x);
}

double NakagamiDistribution::F(const double & x) const
{
    return (x > 0.0) ? Y.F(x * x) : 0.0;
}

double NakagamiDistribution::S(const double & x) const
{
    return (x > 0.0) ? Y.S(x * x) : 1.0;
}

double NakagamiDistribution::Variate() const
{
    return std::sqrt(Y.Variate());
}

void NakagamiDistribution::Sample(std::vector<double> &outputData) const
{
    Y.Sample(outputData);
    for (double & var : outputData)
        var = std::sqrt(var);
}

double NakagamiDistribution::Mean() const
{
    double y = lgammaShapeRatio;
    y += 0.5 * std::log(w / m);
    return std::exp(y);
}

double NakagamiDistribution::Variance() const
{
    double y = lgammaShapeRatio;
    y = std::exp(2 * y);
    return w * (1 - y / m);
}

double NakagamiDistribution::Mode() const
{
    double mode = 0.5 * w / m;
    return std::sqrt(std::max(w - mode, 0.0));
}

double NakagamiDistribution::quantileImpl(double p) const
{
    return std::sqrt(Y.Quantile(p));
}

double NakagamiDistribution::quantileImpl1m(double p) const
{
    return std::sqrt(Y.Quantile1m(p));
}

/// NAKAGAMI
std::string NakagamiRand::Name() const
{
    return "Nakagami(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetSpread()) + ")";
}


/// CHI
std::string ChiRand::Name() const
{
    return "Chi(" + toStringWithPrecision(GetDegree()) +  ")";
}

void ChiRand::SetDegree(int degree)
{
    NakagamiDistribution::SetParameters((degree < 1) ? 0.5 : 0.5 * degree, degree);
}

double ChiRand::Skewness() const
{
    double mu = Mean();
    double sigmaSq = Variance();
    double skew = mu * (1 - 2 * sigmaSq);
    skew /= std::pow(sigmaSq, 1.5);
    return skew;
}

double ChiRand::ExcessKurtosis() const
{
    double mu = Mean();
    double sigmaSq = Variance();
    double sigma = std::sqrt(sigmaSq);
    double skew = Skewness();
    double kurt = 1.0 - mu * sigma * skew;
    kurt /= sigmaSq;
    --kurt;
    return 2 * kurt;
}


/// MAXWELL-BOLTZMANN
MaxwellBoltzmannRand::MaxwellBoltzmannRand(double scale)
{
    SetScale(scale);
}

std::string MaxwellBoltzmannRand::Name() const
{
    return "Maxwell-Boltzmann(" + toStringWithPrecision(GetScale()) + ")";
}

void MaxwellBoltzmannRand::SetScale(double scale)
{
    sigma = (scale > 0.0) ? scale : 1.0;
    NakagamiDistribution::SetParameters(1.5, 3 * sigma * sigma);
}

double MaxwellBoltzmannRand::f(const double & x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x / sigma;
    double xAdjSq = xAdj * xAdj;
    double y = std::exp(-0.5 * xAdjSq);
    return M_SQRT2 * M_1_SQRTPI * xAdjSq * y / sigma;
}

double MaxwellBoltzmannRand::F(const double & x) const
{
    if (x <= 0.0)
        return 0.0;
    double xAdj = M_SQRT1_2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * M_1_SQRTPI;
    return std::erf(xAdj) - y;
}

double MaxwellBoltzmannRand::S(const double & x) const
{
    if (x <= 0.0)
        return 1.0;
    double xAdj = M_SQRT1_2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * M_1_SQRTPI;
    return std::erfc(xAdj) + y;
}

double MaxwellBoltzmannRand::Variate() const
{
    double W = ExponentialRand::StandardVariate();
    double N = NormalRand::StandardVariate();
    return sigma * std::sqrt(2 * W + N * N);
}

void MaxwellBoltzmannRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

double MaxwellBoltzmannRand::Mean() const
{
    return 2 * M_1_SQRTPI * M_SQRT2 * sigma;
}

double MaxwellBoltzmannRand::Variance() const
{
    return (3 - 8.0 * M_1_PI) * sigma * sigma;
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
RayleighRand::RayleighRand(double scale)
{
    SetScale(scale);
}

std::string RayleighRand::Name() const
{
    return "Rayleigh(" + toStringWithPrecision(GetScale()) + ")";
}

void RayleighRand::SetScale(double scale)
{
    sigma = (scale > 0.0) ? scale : 1.0;
    NakagamiDistribution::SetParameters(1, 2 * sigma * sigma);
}

double RayleighRand::f(const double & x) const
{
    if (x <= 0)
        return 0.0;
    double y = x / (sigma * sigma);
    return y * std::exp(-0.5 * x * y);
}

double RayleighRand::F(const double & x) const
{
    if (x <= 0)
        return 0.0;
    double xAdj = x / sigma;
    return -std::expm1(-0.5 * xAdj * xAdj);
}

double RayleighRand::S(const double & x) const
{
    if (x <= 0)
        return 1.0;
    double xAdj = x / sigma;
    return std::exp(-0.5 * xAdj * xAdj);
}

double RayleighRand::Variate() const
{
    double W = ExponentialRand::StandardVariate();
    return sigma * std::sqrt(2 * W);
}

void RayleighRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

double RayleighRand::Mean() const
{
    return sigma * M_SQRTPI * M_SQRT1_2;
}

double RayleighRand::Variance() const
{
    return (2.0 - M_PI_2) * sigma * sigma;
}

double RayleighRand::quantileImpl(double p) const
{
    return sigma * std::sqrt(-2 * std::log1p(-p));
}

double RayleighRand::quantileImpl1m(double p) const
{
    return sigma * std::sqrt(-2 * std::log(p));
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

void RayleighRand::FitScaleMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double sigmaSq = 0.5 * rawMoment(sample, 2);
    SetScale(std::sqrt(sigmaSq));
}

void RayleighRand::FitScaleUMVU(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    size_t n = sample.size();
    double sigmaBiasedSq = 0.5 * rawMoment(sample, 2);
    /// Calculate unbiased sigma
    if (n > 100) {
        double coef = 1.0 / (640 * std::pow(n, 5));
        coef -= 1.0 / (192 * std::pow(n, 3));
        coef += 0.125 / n;
        SetScale((1.0 + coef) * std::sqrt(sigmaBiasedSq)); /// err ~ o(n^{-6.5}) < 1e-13
    }
    else if (n > 15) {
        double coef = RandMath::lfact(n);
        coef += RandMath::lfact(n - 1);
        coef += 2 * n * M_LN2;
        coef += 0.5 * std::log(n);
        coef -= RandMath::lfact(2 * n);
        coef -= 0.5 * M_LNPI;
        coef += 0.5 * std::log(sigmaBiasedSq);
        SetScale(std::exp(coef));
    }
    else
    {
        double scale = RandMath::lfact(n - 1);
        scale += RandMath::lfact(n);
        scale += 0.5 * std::log(n / M_PI * sigmaBiasedSq);
        scale += 2 * n * M_LN2;
        scale -= RandMath::lfact(2 * n);
        SetScale(scale);
    }
}

