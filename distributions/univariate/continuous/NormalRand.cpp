#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "../BasicRandGenerator.h"
#include "GammaRand.h"
#include "StudentTRand.h"

long double NormalRand::stairWidth[257] = {0};
long double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::SetupTables();

NormalRand::NormalRand(double mean, double var)
    : StableDistribution(2.0, 0.0, 1.0, mean)
{
    SetVariance(var);
}

String NormalRand::Name() const
{
    return "Normal(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(Variance()) + ")";
}

void NormalRand::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Normal distribution: scale should be positive");
    sigma = scale;
    StableDistribution::SetScale(sigma * M_SQRT1_2);
}

bool NormalRand::SetupTables()
{
    static constexpr long double A = 4.92867323399e-3l; /// area under rectangle
    /// coordinates of the implicit rectangle in base layer
    stairHeight[0] = 0.001260285930498597l; /// exp(-0.5 * x1 * x1);
    stairWidth[0] = 3.9107579595370918075l; /// A / stairHeight[0];
    /// implicit value for the top layer
    stairWidth[256] = 0.0l;
    stairWidth[1] = x1;
    stairHeight[1] = 0.002609072746106362l;
    for (size_t i = 2; i <= 255; ++i) {
        /// such y_i that f(x_{i+1}) = y_i
        stairWidth[i] = std::sqrt(-2 * std::log(stairHeight[i - 1]));
        stairHeight[i] = stairHeight[i - 1] + A / stairWidth[i];
    }
    return true;
}

void NormalRand::SetVariance(double var)
{
    if (var <= 0.0)
        throw std::invalid_argument("Variance of Normal distribution should be positive");
    SetScale(var);
}

double NormalRand::f(const double & x) const
{
    return pdfNormal(x);
}

double NormalRand::logf(const double & x) const
{
    return logpdfNormal(x);
}

double NormalRand::F(const double & x) const
{
    return cdfNormal(x);
}

double NormalRand::S(const double & x) const
{
    return cdfNormalCompl(x);
}

double NormalRand::Variate() const
{
    return mu + sigma * StandardVariate();
}

double NormalRand::StandardVariate()
{
    /// Ziggurat algorithm by George Marsaglia using 256 strips
    int iter = 0;
    do {
        unsigned long long B = RandGenerator::Variate();
        int stairId = B & 255;
        double x = UniformRand::StandardVariate() * stairWidth[stairId]; /// Get horizontal coordinate
        if (x < stairWidth[stairId + 1])
            return ((signed)B > 0) ? x : -x;
        if (stairId == 0) /// handle the base layer
        {
            static double z = -1;
            if (z > 0) /// we don't have to generate another exponential variable as we already have one
            {
                x = ExponentialRand::StandardVariate() / x1;
                z -= 0.5 * x * x;
            }
            if (z <= 0) /// if previous generation wasn't successful
            {
                do {
                    x = ExponentialRand::StandardVariate() / x1;
                    z = ExponentialRand::StandardVariate() - 0.5 * x * x; /// we storage this value as after acceptance it becomes exponentially distributed
                } while (z <= 0);
            }
            x += x1;
            return ((signed)B > 0) ? x : -x;
        }
        /// handle the wedges of other stairs
        if (UniformRand::Variate(stairHeight[stairId - 1], stairHeight[stairId]) < std::exp(-.5 * x * x))
            return ((signed)B > 0) ? x : -x;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail due to some error
}

void NormalRand::Sample(std::vector<double> &outputData) const
{
    for (double & var : outputData)
        var = this->Variate();
}

std::complex<double> NormalRand::CFImpl(double t) const
{
    return cfNormal(t);
}

double NormalRand::quantileImpl(double p) const
{
    return quantileNormal(p);
}

double NormalRand::quantileImpl1m(double p) const
{
    return quantileNormal1m(p);
}

double NormalRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return (n & 1) ? std::exp(n * this->GetLogScale() + RandMath::ldfact(n - 1)) : 0.0;
}

void NormalRand::FitLocation(const std::vector<double> &sample)
{
    SetLocation(GetSampleMean(sample));
}

void NormalRand::FitLocation(const std::vector<double> &sample, DoublePair &confidenceInterval, double significanceLevel)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(fitErrorDescription(WRONG_LEVEL, "Input level is equal to " + toStringWithPrecision(significanceLevel)));

    FitLocation(sample);

    int n = sample.size();
    NormalRand NormalRV(0, 1);
    double halfAlpha = 0.5 * significanceLevel;
    double interval = NormalRV.Quantile1m(halfAlpha) * sigma / std::sqrt(n);
    confidenceInterval.first = mu - interval;
    confidenceInterval.second = mu + interval;
}

void NormalRand::FitVariance(const std::vector<double> &sample)
{
    SetVariance(GetSampleVariance(sample, mu));
}

void NormalRand::FitVariance(const std::vector<double> &sample, DoublePair &confidenceInterval, double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(fitErrorDescription(WRONG_LEVEL, "Input level is equal to " + toStringWithPrecision(significanceLevel)));

    FitVariance(sample);

    size_t n = sample.size();
    double halfAlpha = 0.5 * significanceLevel;
    ChiSquaredRand ChiSqRV(n);
    double numerator = sigma * sigma * (unbiased ? n : (n - 1));
    confidenceInterval.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceInterval.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

void NormalRand::FitScale(const std::vector<double> &sample, bool unbiased)
{
    if (unbiased == true) {
        size_t n = sample.size();
        double halfN = 0.5 * n;
        double s = GetSampleVariance(sample, mu);
        s *= halfN;
        s = std::log(s);
        s += std::log(halfN);
        s *= 0.5;
        s += std::lgamma(halfN + 0.5);
        s -= std::lgamma(halfN);
        SetScale(std::exp(s));
    }
    else {
        FitVariance(sample);
    }
}

void NormalRand::Fit(const std::vector<double> &sample, bool unbiased)
{
    double adjustment = 1.0;
    if (unbiased == true) {
        size_t n = sample.size();
        if (n <= 1)
            throw std::invalid_argument(fitErrorDescription(TOO_FEW_ELEMENTS, "There should be at least 2 elements"));
        adjustment = static_cast<double>(n) / (n - 1);
    }
    DoublePair stats = GetSampleMeanAndVariance(sample);
    SetLocation(stats.first);
    SetVariance(stats.second * adjustment);
}

void NormalRand::Fit(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double significanceLevel, bool unbiased)
{
    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(fitErrorDescription(WRONG_LEVEL, "Input level is equal to " + toStringWithPrecision(significanceLevel)));

    Fit(sample, unbiased);

    size_t n = sample.size();
    double sigmaAdj = unbiased ? sigma : (sigma * n) / (n - 1);
    /// calculate confidence interval for mean
    double halfAlpha = 0.5 * significanceLevel;
    StudentTRand tRV(n - 1);
    double interval = tRV.Quantile1m(halfAlpha) * sigmaAdj / std::sqrt(n);
    confidenceIntervalForMean.first = mu - interval;
    confidenceIntervalForMean.second = mu + interval;

    /// calculate confidence interval for variance
    ChiSquaredRand ChiSqRV(n - 1);
    double numerator = (n - 1) * sigmaAdj * sigmaAdj;
    confidenceIntervalForVariance.first = numerator / ChiSqRV.Quantile1m(halfAlpha);
    confidenceIntervalForVariance.second = numerator / ChiSqRV.Quantile(halfAlpha);
}

NormalRand NormalRand::FitLocationBayes(const std::vector<double> &sample, const NormalRand &priorDistribution)
{
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = GetPrecision();
    double numerator = GetSampleSum(sample) * tau + tau0 * mu0;
    double denominator = sample.size() * tau + tau0;
    NormalRand posteriorDistribution(numerator / denominator, 1.0 / denominator);
    SetLocation(posteriorDistribution.Mean());
    return posteriorDistribution;
}

InverseGammaRand NormalRand::FitVarianceBayes(const std::vector<double> &sample, const InverseGammaRand &priorDistribution)
{
    double halfN = 0.5 * sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double alphaPosterior = alphaPrior + halfN;
    double betaPosterior = betaPrior + halfN * GetSampleVariance(sample, mu);
    InverseGammaRand posteriorDistribution(alphaPosterior, betaPosterior);
    SetVariance(posteriorDistribution.Mean());
    return posteriorDistribution;
}

NormalInverseGammaRand NormalRand::FitBayes(const std::vector<double> &sample, const NormalInverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    double alphaPrior = priorDistribution.GetShape();
    double betaPrior = priorDistribution.GetRate();
    double muPrior = priorDistribution.GetLocation();
    double lambdaPrior = priorDistribution.GetPrecision();
    DoublePair stats = GetSampleMeanAndVariance(sample);
    double lambdaPosterior = lambdaPrior + n;
    double muPosterior = (lambdaPrior * muPrior + n * stats.first) / lambdaPosterior;
    double halfN = 0.5 * n;
    double alphaPosterior = alphaPrior + halfN;
    double aux = muPrior - stats.first;
    double betaPosterior = betaPrior + halfN * (stats.second + lambdaPrior / lambdaPosterior * aux * aux);
    NormalInverseGammaRand posteriorDistribution(muPosterior, lambdaPosterior, alphaPosterior, betaPosterior);
    DoublePair mean = posteriorDistribution.Mean();
    SetLocation(mean.first);
    SetVariance(mean.second);
    return posteriorDistribution;
}
