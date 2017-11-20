#include "NakagamiRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

template < typename RealType >
NakagamiDistribution<RealType>::NakagamiDistribution(double shape, double spread)
{
    SetParameters(shape, spread);
}

template < typename RealType >
void NakagamiDistribution<RealType>::SetParameters(double shape, double spread)
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

template < typename RealType >
double NakagamiDistribution<RealType>::f(const RealType & x) const
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

template < typename RealType >
double NakagamiDistribution<RealType>::logf(const RealType & x) const
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

template < typename RealType >
double NakagamiDistribution<RealType>::F(const RealType & x) const
{
    return (x > 0.0) ? Y.F(x * x) : 0.0;
}

template < typename RealType >
double NakagamiDistribution<RealType>::S(const RealType & x) const
{
    return (x > 0.0) ? Y.S(x * x) : 1.0;
}

template < typename RealType >
RealType NakagamiDistribution<RealType>::Variate() const
{
    return std::sqrt(Y.Variate());
}

template < typename RealType >
void NakagamiDistribution<RealType>::Sample(std::vector<RealType> &outputData) const
{
    Y.Sample(outputData);
    for (RealType & var : outputData)
        var = std::sqrt(var);
}

template < typename RealType >
void NakagamiDistribution<RealType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    Y.Reseed(seed + 1);
}

template < typename RealType >
long double NakagamiDistribution<RealType>::Mean() const
{
    long double y = lgammaShapeRatio;
    y -= 0.5 * Y.GetLogRate();
    return std::exp(y);
}

template < typename RealType >
long double NakagamiDistribution<RealType>::Variance() const
{
    long double y = lgammaShapeRatio;
    y = std::exp(2 * y);
    return omega * (1 - y / mu);
}

template < typename RealType >
RealType NakagamiDistribution<RealType>::Median() const
{
    return std::sqrt(Y.Quantile(0.5));
}

template < typename RealType >
RealType NakagamiDistribution<RealType>::Mode() const
{
    long double mode = 1.0 - 0.5 / mu;
    if (mode <= 0.0)
        return 0.0;
    return std::sqrt(omega * mode);
}

template < typename RealType >
long double NakagamiDistribution<RealType>::Skewness() const
{
    long double thirdMoment = lgammaShapeRatio;
    thirdMoment -= 1.5 * Y.GetLogRate();
    thirdMoment = (mu + 0.5) * std::exp(thirdMoment);
    long double mean = Mean();
    long double variance = Variance();
    return (thirdMoment - mean * (3 * variance + mean * mean)) / std::pow(variance, 1.5);
}

template < typename RealType >
long double NakagamiDistribution<RealType>::FourthMoment() const
{
    long double fourthMoment = omega / mu;
    fourthMoment *= fourthMoment;
    fourthMoment *= mu * (mu + 1);
    return fourthMoment;
}

template < typename RealType >
long double NakagamiDistribution<RealType>::ExcessKurtosis() const
{
    long double mean = Mean();
    long double secondMoment = this->SecondMoment();
    long double thirdMoment = this->ThirdMoment();
    long double fourthMoment = FourthMoment();
    long double meanSq = mean * mean;
    long double variance = secondMoment - meanSq;
    long double numerator = fourthMoment - 4 * thirdMoment * mean + 6 * secondMoment * meanSq - 3 * meanSq * meanSq;
    long double denominator = variance * variance;
    return numerator / denominator - 3.0;
}

template < typename RealType >
RealType NakagamiDistribution<RealType>::quantileImpl(double p) const
{
    return std::sqrt(Y.Quantile(p));
}

template < typename RealType >
RealType NakagamiDistribution<RealType>::quantileImpl1m(double p) const
{
    return std::sqrt(Y.Quantile1m(p));
}

template < typename RealType >
std::complex<double> NakagamiDistribution<RealType>::CFImpl(double t) const
{
    if (mu >= 0.5)
        return ContinuousDistribution<RealType>::CFImpl(t);

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

template class NakagamiDistribution<float>;
template class NakagamiDistribution<double>;
template class NakagamiDistribution<long double>;

/// NAKAGAMI
template < typename RealType >
String NakagamiRand<RealType>::Name() const
{
    return "Nakagami(" + this->toStringWithPrecision(this->GetShape()) + ", " + this->toStringWithPrecision(this->GetSpread()) + ")";
}

template class NakagamiRand<float>;
template class NakagamiRand<double>;
template class NakagamiRand<long double>;

/// CHI
template < typename RealType >
ChiRand<RealType>::ChiRand(int degree)
{
    SetDegree(degree);
}

template < typename RealType >
String ChiRand<RealType>::Name() const
{
    return "Chi(" + this->toStringWithPrecision(GetDegree()) +  ")";
}

template < typename RealType >
void ChiRand<RealType>::SetDegree(int degree)
{
    if (degree < 1)
        throw std::invalid_argument("Chi distribution: degree parameter should be positive");
    NakagamiDistribution<RealType>::SetParameters(0.5 * degree, degree);
}

template < typename RealType >
long double ChiRand<RealType>::Skewness() const
{
    long double mean = this->Mean();
    long double sigmaSq = this->Variance();
    long double skew = mean * (1 - 2 * sigmaSq);
    skew /= std::pow(sigmaSq, 1.5);
    return skew;
}

template < typename RealType >
long double ChiRand<RealType>::ExcessKurtosis() const
{
    long double mean = this->Mean();
    long double sigmaSq = this->Variance();
    long double sigma = std::sqrt(sigmaSq);
    long double skew = Skewness();
    long double kurt = 1.0 - mean * sigma * skew;
    kurt /= sigmaSq;
    --kurt;
    return 2 * kurt;
}

template class ChiRand<float>;
template class ChiRand<double>;
template class ChiRand<long double>;

/// MAXWELL-BOLTZMANN
template < typename RealType >
MaxwellBoltzmannRand<RealType>::MaxwellBoltzmannRand(double scale)
{
    SetScale(scale);
}

template < typename RealType >
String MaxwellBoltzmannRand<RealType>::Name() const
{
    return "Maxwell-Boltzmann(" + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void MaxwellBoltzmannRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Maxwell-Boltzmann distribution: scale should be positive");
    sigma = scale;
    NakagamiDistribution<RealType>::SetParameters(1.5, 3 * sigma * sigma);
}

template < typename RealType >
double MaxwellBoltzmannRand<RealType>::f(const RealType &x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x / sigma;
    double xAdjSq = xAdj * xAdj;
    double y = std::exp(-0.5 * xAdjSq);
    return M_SQRT2 * M_1_SQRTPI * xAdjSq * y / sigma;
}

template < typename RealType >
double MaxwellBoltzmannRand<RealType>::F(const RealType & x) const
{
    if (x <= 0.0)
        return 0.0;
    double xAdj = M_SQRT1_2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * M_1_SQRTPI;
    return std::erf(xAdj) - y;
}

template < typename RealType >
double MaxwellBoltzmannRand<RealType>::S(const RealType &x) const
{
    if (x <= 0.0)
        return 1.0;
    double xAdj = M_SQRT1_2 * x / sigma;
    double y = std::exp(-xAdj * xAdj);
    y *= 2 * xAdj * M_1_SQRTPI;
    return std::erfc(xAdj) + y;
}

template < typename RealType >
RealType MaxwellBoltzmannRand<RealType>::Variate() const
{
    RealType W = ExponentialRand::StandardVariate(this->localRandGenerator);
    RealType N = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
    return sigma * std::sqrt(2 * W + N * N);
}

template < typename RealType >
void MaxwellBoltzmannRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    for (RealType & var : outputData)
        var = this->Variate();
}

template < typename RealType >
long double MaxwellBoltzmannRand<RealType>::Mean() const
{
    return 2 * M_1_SQRTPI * M_SQRT2 * sigma;
}

template < typename RealType >
long double MaxwellBoltzmannRand<RealType>::Variance() const
{
    return (3 - 8.0 * M_1_PI) * sigma * sigma;
}

template < typename RealType >
RealType MaxwellBoltzmannRand<RealType>::Mode() const
{
    return M_SQRT2 * sigma;
}

template < typename RealType >
long double MaxwellBoltzmannRand<RealType>::Skewness() const
{
    long double skewness = 3 * M_PI - 8;
    skewness = 2.0 / skewness;
    skewness *= std::sqrt(skewness);
    return (16 - 5 * M_PI) * skewness;
}

template < typename RealType >
long double MaxwellBoltzmannRand<RealType>::ExcessKurtosis() const
{
    long double numerator = 40 - 3 * M_PI;
    numerator *= M_PI;
    numerator -= 96;
    long double denominator = 3 * M_PI - 8;
    denominator *= denominator;
    return 4 * numerator / denominator;
}

template class MaxwellBoltzmannRand<float>;
template class MaxwellBoltzmannRand<double>;
template class MaxwellBoltzmannRand<long double>;

/// RAYLEIGH
template < typename RealType >
RayleighRand<RealType>::RayleighRand(double scale)
{
    SetScale(scale);
}

template < typename RealType >
String RayleighRand<RealType>::Name() const
{
    return "Rayleigh(" + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void RayleighRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Rayleigh distribution: scale should be positive");
    NakagamiDistribution<RealType>::SetParameters(1, 2 * sigma * sigma);
}

template < typename RealType >
double RayleighRand<RealType>::f(const RealType & x) const
{
    if (x <= 0)
        return 0.0;
    double y = x / (sigma * sigma);
    return y * std::exp(-0.5 * x * y);
}

template < typename RealType >
double RayleighRand<RealType>::F(const RealType & x) const
{
    if (x <= 0)
        return 0.0;
    double xAdj = x / sigma;
    return -std::expm1l(-0.5 * xAdj * xAdj);
}

template < typename RealType >
double RayleighRand<RealType>::S(const RealType & x) const
{
    if (x <= 0)
        return 1.0;
    double xAdj = x / sigma;
    return std::exp(-0.5 * xAdj * xAdj);
}

template < typename RealType >
RealType RayleighRand<RealType>::Variate() const
{
    double W = ExponentialRand::StandardVariate(this->localRandGenerator);
    return sigma * std::sqrt(2 * W);
}

template < typename RealType >
void RayleighRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    for (RealType & var : outputData)
        var = this->Variate();
}

template < typename RealType >
long double RayleighRand<RealType>::Mean() const
{
    return sigma * M_SQRTPI * M_SQRT1_2;
}

template < typename RealType >
long double RayleighRand<RealType>::Variance() const
{
    return (2.0 - M_PI_2) * sigma * sigma;
}

template < typename RealType >
RealType RayleighRand<RealType>::quantileImpl(double p) const
{
    return sigma * std::sqrt(-2 * std::log1pl(-p));
}

template < typename RealType >
RealType RayleighRand<RealType>::quantileImpl1m(double p) const
{
    return sigma * std::sqrt(-2 * std::log(p));
}

template < typename RealType >
RealType RayleighRand<RealType>::Median() const
{
    static constexpr double medianCoef = std::sqrt(M_LN2 + M_LN2);
    return sigma * medianCoef;
}

template < typename RealType >
RealType RayleighRand<RealType>::Mode() const
{
    return sigma;
}

template < typename RealType >
long double RayleighRand<RealType>::Skewness() const
{
    static constexpr long double skewness = 2 * M_SQRTPI * (M_PI - 3) / std::pow(4.0 - M_PI, 1.5);
    return skewness;
}

template < typename RealType >
long double RayleighRand<RealType>::ExcessKurtosis() const
{
    static constexpr long double temp = 4 - M_PI;
    static constexpr long double kurtosis = -(6 * M_PI * M_PI - 24 * M_PI + 16) / (temp * temp);
    return kurtosis;
}

template < typename RealType >
void RayleighRand<RealType>::FitScale(const std::vector<RealType> &sample, bool unbiased)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    size_t n = sample.size();
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
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

template class RayleighRand<float>;
template class RayleighRand<double>;
template class RayleighRand<long double>;
