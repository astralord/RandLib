#include "NakagamiRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

NakagamiDistribution::NakagamiDistribution(double shape, double spread)
{
    SetParameters(shape, spread);
}

void NakagamiDistribution::SetParameters(double shape, double spread)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Nakagami distribution: shape should be positive");
    if (spread <= 0.0)
        throw std::invalid_argument("Nakagami distribution: spread should be positive");
    mu = shape;
    omega = spread;
    Y.SetParameters(mu, mu / omega);
    lgammaShapeRatio = std::lgammal(mu + 0.5) - Y.GetLogGammaShape();
}

double NakagamiDistribution::f(const double & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0) {
        if (mu > 0.5)
            return 0.0;
        return (mu < 0.5) ? INFINITY : std::sqrt(M_2_PI / omega);
    }
    return 2 * x * Y.f(x * x);
}

double NakagamiDistribution::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0) {
        if (mu > 0.5)
            return 0-INFINITY;
        return (mu < 0.5) ? INFINITY : 0.5 * (M_LN2 - M_LNPI - std::log(omega));
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

void NakagamiDistribution::Reseed(unsigned long seed) const
{
    localRandGenerator.Reseed(seed);
    Y.Reseed(seed + 1);
}

long double NakagamiDistribution::Mean() const
{
    double y = lgammaShapeRatio;
    y -= 0.5 * Y.GetLogRate();
    return std::exp(y);
}

long double NakagamiDistribution::Variance() const
{
    double y = lgammaShapeRatio;
    y = std::exp(2 * y);
    return omega * (1 - y / mu);
}

double NakagamiDistribution::Median() const
{
    return std::sqrt(Y.Quantile(0.5));
}

double NakagamiDistribution::Mode() const
{
    double mode = 0.5 * omega / mu;
    return std::sqrt(std::max(omega - mode, 0.0));
}

long double NakagamiDistribution::Skewness() const
{
    double thirdMoment = lgammaShapeRatio;
    thirdMoment -= 1.5 * Y.GetLogRate();
    thirdMoment = (mu + 0.5) * std::exp(thirdMoment);
    double mean = Mean();
    double variance = Variance();
    return (thirdMoment - mean * (3 * variance + mean * mean)) / std::pow(variance, 1.5);
}

long double NakagamiDistribution::FourthMoment() const
{
    double fourthMoment = omega / mu;
    fourthMoment *= fourthMoment;
    fourthMoment *= mu * (mu + 1);
    return fourthMoment;
}

long double NakagamiDistribution::ExcessKurtosis() const
{
    double mean = Mean();
    double secondMoment = SecondMoment();
    double thirdMoment = ThirdMoment();
    double fourthMoment = FourthMoment();
    double meanSq = mean * mean;
    double variance = secondMoment - meanSq;
    double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

double NakagamiDistribution::quantileImpl(double p) const
{
    return std::sqrt(Y.Quantile(p));
}

double NakagamiDistribution::quantileImpl1m(double p) const
{
    return std::sqrt(Y.Quantile1m(p));
}

std::complex<double> NakagamiDistribution::CFImpl(double t) const
{
    if (mu >= 0.5)
        return ContinuousDistribution::CFImpl(t);

    double re = this->ExpectedValue([this, t] (double x)
    {
        if (x == 0.0)
            return 0.0;
        return std::cos(t * x) - 1.0;
    }, 0, INFINITY) + 1.0;

    double im = this->ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, 0, INFINITY);

    return std::complex<double>(re, im);
}

/// NAKAGAMI
String NakagamiRand::Name() const
{
    return "Nakagami(" + this->toStringWithPrecision(GetShape()) + ", " + this->toStringWithPrecision(GetSpread()) + ")";
}


/// CHI
ChiRand::ChiRand(int degree)
{
    SetDegree(degree);
}

String ChiRand::Name() const
{
    return "Chi(" + this->toStringWithPrecision(GetDegree()) +  ")";
}

void ChiRand::SetDegree(int degree)
{
    if (degree < 1)
        throw std::invalid_argument("Chi distribution: degree parameter should be positive");
    NakagamiDistribution::SetParameters(0.5 * degree, degree);
}

long double ChiRand::Skewness() const
{
    double mean = Mean();
    double sigmaSq = Variance();
    double skew = mean * (1 - 2 * sigmaSq);
    skew /= std::pow(sigmaSq, 1.5);
    return skew;
}

long double ChiRand::ExcessKurtosis() const
{
    double mean = Mean();
    double sigmaSq = Variance();
    double sigma = std::sqrt(sigmaSq);
    double skew = Skewness();
    double kurt = 1.0 - mean * sigma * skew;
    kurt /= sigmaSq;
    --kurt;
    return 2 * kurt;
}


/// MAXWELL-BOLTZMANN
MaxwellBoltzmannRand::MaxwellBoltzmannRand(double scale)
{
    SetScale(scale);
}

String MaxwellBoltzmannRand::Name() const
{
    return "Maxwell-Boltzmann(" + this->toStringWithPrecision(GetScale()) + ")";
}

void MaxwellBoltzmannRand::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Maxwell-Boltzmann distribution: scale should be positive");
    sigma = scale;
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
    double W = ExponentialRand::StandardVariate(localRandGenerator);
    double N = NormalRand::StandardVariate(localRandGenerator);
    return sigma * std::sqrt(2 * W + N * N);
}

void MaxwellBoltzmannRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

long double MaxwellBoltzmannRand::Mean() const
{
    return 2 * M_1_SQRTPI * M_SQRT2 * sigma;
}

long double MaxwellBoltzmannRand::Variance() const
{
    return (3 - 8.0 * M_1_PI) * sigma * sigma;
}

double MaxwellBoltzmannRand::Mode() const
{
    return M_SQRT2 * sigma;
}

long double MaxwellBoltzmannRand::Skewness() const
{
    double skewness = 3 * M_PI - 8;
    skewness = 2.0 / skewness;
    skewness *= std::sqrt(skewness);
    return (16 - 5 * M_PI) * skewness;
}

long double MaxwellBoltzmannRand::ExcessKurtosis() const
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

String RayleighRand::Name() const
{
    return "Rayleigh(" + this->toStringWithPrecision(GetScale()) + ")";
}

void RayleighRand::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Rayleigh distribution: scale should be positive");
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
    return -std::expm1l(-0.5 * xAdj * xAdj);
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
    double W = ExponentialRand::StandardVariate(localRandGenerator);
    return sigma * std::sqrt(2 * W);
}

void RayleighRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

long double RayleighRand::Mean() const
{
    return sigma * M_SQRTPI * M_SQRT1_2;
}

long double RayleighRand::Variance() const
{
    return (2.0 - M_PI_2) * sigma * sigma;
}

double RayleighRand::quantileImpl(double p) const
{
    return sigma * std::sqrt(-2 * std::log1pl(-p));
}

double RayleighRand::quantileImpl1m(double p) const
{
    return sigma * std::sqrt(-2 * std::log(p));
}

double RayleighRand::Median() const
{
    static constexpr double medianCoef = std::sqrt(M_LN2 + M_LN2);
    return sigma * medianCoef;
}

double RayleighRand::Mode() const
{
    return sigma;
}

long double RayleighRand::Skewness() const
{
    static constexpr long double skewness = 2 * M_SQRTPI * (M_PI - 3) / std::pow(4.0 - M_PI, 1.5);
    return skewness;
}

long double RayleighRand::ExcessKurtosis() const
{
    static constexpr long double temp = 4 - M_PI;
    static constexpr long double kurtosis = -(6 * M_PI * M_PI - 24 * M_PI + 16) / (temp * temp);
    return kurtosis;
}

void RayleighRand::FitScale(const std::vector<double> &sample, bool unbiased)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    size_t n = sample.size();
    DoublePair stats = GetSampleMeanAndVariance(sample);
    double rawMoment = stats.second + stats.first * stats.first;
    double sigmaBiasedSq = 0.5 * rawMoment;
    if (unbiased == false) {
        SetScale(std::sqrt(sigmaBiasedSq));
    }
    /// Calculate unbiased sigma
    else if (n > 100) {
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
