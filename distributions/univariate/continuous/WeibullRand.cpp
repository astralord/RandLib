#include "WeibullRand.h"
#include "ExponentialRand.h"

WeibullRand::WeibullRand(double scale, double shape)
{
    SetParameters(scale, shape);
}

String WeibullRand::Name() const
{
    return "Weibull(" + this->toStringWithPrecision(GetScale()) + ", " + this->toStringWithPrecision(GetShape()) + ")";
}

void WeibullRand::SetParameters(double scale, double shape)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Weibull distribution: scale should be positive");
    if (shape <= 0.0)
        throw std::invalid_argument("Weibull distribution: shape should be positive");
    lambda = scale;
    k = shape;
    kInv = 1.0 / k;
    logk_lambda = std::log(k / lambda);
}

double WeibullRand::f(const double & x) const
{
    if (x < 0)
        return 0;
    if (x == 0) {
        if (k == 1)
            return 1.0 / lambda;
        return (k > 1) ? 0.0 : INFINITY;
    }
    double xAdj = x / lambda;
    double xAdjPow = std::pow(xAdj, k - 1);
    return k / lambda * xAdjPow * std::exp(-xAdj * xAdjPow);
}

double WeibullRand::logf(const double & x) const
{
    if (x < 0)
        return -INFINITY;
    if (x == 0) {
        if (k == 1)
            return logk_lambda;
        return (k > 1) ? -INFINITY : INFINITY;
    }
    double xAdj = x / lambda;
    double xAdjPow = std::pow(xAdj, k - 1);
    return logk_lambda + (k - 1) * std::log(xAdj) - xAdj * xAdjPow;
}

double WeibullRand::F(const double & x) const
{
    return (x > 0.0) ? -std::expm1l(-std::pow(x / lambda, k)) : 0.0;
}

double WeibullRand::S(const double & x) const
{
    return (x > 0.0) ? std::exp(-std::pow(x / lambda, k)) : 1.0;
}

double WeibullRand::Variate() const
{
    return lambda * std::pow(ExponentialRand::StandardVariate(localRandGenerator), kInv);
}

long double WeibullRand::Mean() const
{
    return lambda * std::tgammal(1 + kInv);
}

long double WeibullRand::Variance() const
{
    double res = std::tgammal(1 + kInv);
    res *= -res;
    res += std::tgammal(1 + kInv + kInv);
    return lambda * lambda * res;
}

double WeibullRand::Median() const
{
    return lambda * std::pow(M_LN2, kInv);
}

double WeibullRand::Mode() const
{
    if (k <= 1)
        return 0;
    double y = std::log1pl(-kInv);
    y = std::exp(kInv * y);
    return lambda * y;
}

long double WeibullRand::Skewness() const
{
    long double mu = Mean();
    long double var = Variance();
    long double sigma = std::sqrt(var);
    long double numerator = std::tgammal(1 + 3.0 * kInv);
    numerator *= lambda * lambda * lambda;
    numerator -= 3 * mu * var;
    numerator -= mu * mu * mu;
    long double denominator = var * sigma;
    return numerator / denominator;
}

long double WeibullRand::ExcessKurtosis() const
{
    long double mu = Mean();
    long double var = Variance();
    long double sigma = std::sqrt(var);
    long double skewness = Skewness();
    long double numerator = lambda * lambda;
    numerator *= numerator;
    numerator *= std::tgammal(1 + 4.0 * kInv);
    numerator -= 4 * skewness * var * sigma * mu;
    long double mu2 = mu * mu;
    numerator -= 6 * mu2 * var;
    numerator -= mu2 * mu2;
    long double kurtosis = numerator / (var * var);
    return kurtosis - 3;
}

double WeibullRand::quantileImpl(double p) const
{
    double x = -std::log1pl(-p);
    x = std::pow(x, kInv);
    return lambda * x;
}

double WeibullRand::quantileImpl1m(double p) const
{
    double x = -std::log(p);
    x = std::pow(x, kInv);
    return lambda * x;
}

std::complex<double> WeibullRand::CFImpl(double t) const
{
    double lambdaT = lambda * t;
    if (k >= 1) {
        if (lambdaT > 0.5)
            return ContinuousDistribution::CFImpl(t);
        /// for λt < 0.5, the worst case scenario for series expansion is n ~ 70
        double re = 0.0, im = 0.0;
        double addon = 0.0;
        double logLambdaT = std::log(lambdaT);
        /// Series representation for real part
        int n = 0;
        do {
            int n2 = n + n;
            addon = n2 * logLambdaT;
            addon += std::lgammal(1.0 + n2 / k);
            addon -= std::lgammal(1.0 + n2);
            addon = std::exp(addon);
            re += (n & 1) ? -addon : addon;
            ++n;
        } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(re));
        /// Series representation for imaginary part
        n = 0;
        do {
            int n2p1 = n + n + 1;
            addon = n2p1 * logLambdaT;
            addon += std::lgammal(1.0 + n2p1 / k);
            addon -= std::lgammal(1.0 + n2p1);
            addon = std::exp(addon);
            im += (n & 1) ? -addon : addon;
            ++n;
        } while (std::fabs(addon) > MIN_POSITIVE * std::fabs(im));
        return std::complex<double>(re, im);
    }

    /// For real part with k < 1 we split the integral on two intervals
    double re1 = RandMath::integral([this, t] (double x)
    {
        if (x <= 0.0 || x > 1.0)
            return 0.0;
        double xAdj = x / lambda;
        double xAdjPow = std::pow(xAdj, k - 1);
        double y = k / lambda * xAdjPow * std::expm1l(-xAdj * xAdjPow);
        return std::cos(t * x) * y;
    },
    0.0, 1.0);

    double re2 = this->ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    },
    1.0, INFINITY);

    double re3 = t * RandMath::integral([this, t] (double x)
    {
        if (x <= 0.0)
            return 0.0;
        return std::sin(t * x) * std::pow(x, k);
    },
    0.0, 1.0);
    re3 += std::cos(t);
    re3 /= std::pow(lambda, k);

    double re = re1 + re2 + re3;

    double im = this->ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    },
    0.0, INFINITY);

    return std::complex<double>(re, im);
}

double WeibullRand::Entropy() const
{
    return M_EULER * (1.0 - kInv) + std::log(lambda * kInv) + 1.0;
}

double WeibullRand::getPowSampleMean(const std::vector<double> &sample) const
{
    long double sum = 0;
    for (const double &var : sample) {
        sum += std::pow(var, k);
    }
    return sum / sample.size();
}

void WeibullRand::FitScale(const std::vector<double> &sample)
{
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    double powScale = getPowSampleMean(sample);
    SetParameters(std::pow(powScale, kInv), k);
}

InverseGammaRand WeibullRand::FitScaleBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution, bool MAP)
{
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newRate = priorDistribution.GetRate() + n * getPowSampleMean(sample);
    InverseGammaRand posteriorDistribution(newShape, newRate);
    double powScale =  MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean();
    SetParameters(std::pow(powScale, kInv), k);
    return posteriorDistribution;
}
