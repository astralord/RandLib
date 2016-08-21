#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "../BasicRandGenerator.h"
#include "GammaRand.h"
#include "StudentTRand.h"

double NormalRand::stairWidth[257] = {0};
double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::SetupTables();

NormalRand::NormalRand(double mean, double var)
    : StableRand(2.0, 0.0, 1.0, mean)
{
    SetVariance(var);
}

std::string NormalRand::Name() const
{
    return "Normal(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(Variance()) + ")";
}

void NormalRand::SetScale(double scale)
{
    sigma0 = scale;
    if (sigma0 <= 0.0)
        sigma0 = 1.0;
    StableRand::SetScale(scale * M_SQRT1_2);
}

bool NormalRand::SetupTables()
{
    static constexpr long double A = 4.92867323399e-3l; /// area under rectangle

    /// coordinates of the implicit rectangle in base layer
    stairHeight[0] = std::exp(-.5 * x1 * x1);
    stairWidth[0] = A / stairHeight[0];
    /// implicit value for the top layer
    stairWidth[256] = 0;

    for (size_t i = 1; i <= 255; ++i)
    {
        /// such y_i that f(x_{i+1}) = y_i
        stairWidth[i] = std::sqrt(-2 * std::log(stairHeight[i - 1]));
        stairHeight[i] = stairHeight[i - 1] + A / stairWidth[i];
    }
    return true;
}

void NormalRand::SetVariance(double var)
{
    if (var <= 0)
        var = 1.0;
    SetScale(std::sqrt(var));
}

double NormalRand::f(double x) const
{
    return StableRand::pdfNormal(x);
}

double NormalRand::F(double x) const
{
    return StableRand::cdfNormal(x);
}

double NormalRand::Variate() const
{
    return mu + sigma0 * StandardVariate();
}

double NormalRand::Variate(double mean, double rootVar)
{
    return mean + rootVar * StandardVariate();
}

double NormalRand::StandardVariate()
{
    /// Ziggurat algorithm by George Marsaglia using 256 strips
    int iter = 0;
    do {
        unsigned long long B = RandGenerator::variate();
        int stairId = B & 255;
        double x = UniformRand::StandardVariate() * stairWidth[stairId]; /// Get horizontal coordinate

        if (x < stairWidth[stairId + 1])
            return ((signed)B > 0) ? x : -x;

        if (stairId == 0) /// handle the base layer
        {
            static double z = -1;

            if (z > 0) /// we don't have to generate another exponential variable as we already have one
            {
                x = ExponentialRand::Variate(x1);
                z -= 0.5 * x * x;
            }

            if (z <= 0) /// if previous generation wasn't successful
            {
                do {
                    x = ExponentialRand::Variate(x1);
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

std::complex<double> NormalRand::CF(double t) const
{
    if (t == 0)
        return 1;
    double sigma0T = sigma0 * t;
    return std::exp(std::complex<double>(-0.5 * sigma0T * sigma0T, mu * t));
}

double NormalRand::standardQuantile(double p)
{
    // TODO: redo by erfinv
    /// Abramowitz and Stegun improved approximation
    long double t = (p < 0.5) ? -std::log(p) : -std::log1p(-p);
    t = std::sqrt(t + t);
    static constexpr long double c[] = {2.653962002601684482l, 1.561533700212080345l, 0.061146735765196993l};
    static constexpr long double d[] = {1.904875182836498708l, 0.454055536444233510l, 0.009547745327068945l};
    long double numerator = c[2] * t;
    numerator += c[1];
    numerator *= t;
    numerator += c[0];
    long double denominator = d[2] * t;
    denominator += d[1];
    denominator *= t;
    denominator += d[0];
    denominator *= t;
    denominator += 1.0;
    double y = numerator / denominator - t;
    return (p < 0.5) ? y : -y;
}

double NormalRand::quantile(double p, double mean, double scale)
{
    if (p < 0 || p > 1)
        return NAN;
    return mean + scale * standardQuantile(p);
}

double NormalRand::quantileImpl(double p) const
{
    return mu + sigma0 * standardQuantile(p);
}

double NormalRand::quantileImpl1m(double p) const
{
    // TODO: redo by erfcinv
    return mu + sigma0 * standardQuantile(1.0 - p);
}

double NormalRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return (n & 1) ? std::pow(sigma0, n) * RandMath::doubleFactorial(n - 1) : 0;
}

bool NormalRand::FitMeanMLE(const std::vector<double> &sample)
{
    SetLocation(sampleMean(sample));
    return true;
}

bool NormalRand::FitVarianceMLE(const std::vector<double> &sample)
{
    SetVariance(sampleVariance(sample));
    return true;
}

bool NormalRand::FitMLE(const std::vector<double> &sample)
{
    return FitMeanMLE(sample) ? FitVarianceMLE(sample) : false;
}

bool NormalRand::FitMeanMM(const std::vector<double> &sample)
{
    return FitMeanMLE(sample);
}

bool NormalRand::FitVarianceMM(const std::vector<double> &sample)
{
    return FitVarianceMLE(sample);
}

bool NormalRand::FitMM(const std::vector<double> &sample)
{
    return FitMLE(sample);
}

bool NormalRand::FitMeanUMVU(const std::vector<double> &sample)
{
    return FitMeanMLE(sample);
}

bool NormalRand::FitVarianceUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 1)
        return false;
    double s = sampleVariance(sample, mu);
    SetVariance(n * s / (n - 1));
    return true;
}

bool NormalRand::FitUMVU(const std::vector<double> &sample)
{
    return FitMeanMLE(sample) ? FitVarianceUMVU(sample) : false;
}

bool NormalRand::FitUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha)
{
    size_t n = sample.size();
    if (n < 2 || alpha <= 0 || alpha > 1)
        return false;
    if (!FitUMVU(sample))
        return false;
    double p = 1.0 - 0.5 * alpha;
    size_t nm1 = n - 1;

    /// calculate confidence interval for mean
    StudentTRand t(nm1);
    double interval = t.Quantile(p) * sigma0 / std::sqrt(n);
    confidenceIntervalForMean.first = mu - interval;
    confidenceIntervalForMean.second = mu + interval;

    /// calculate confidence interval for variance
    ChiSquaredRand chi(nm1);
    double numerator = nm1 * sigma0 * sigma0;
    confidenceIntervalForVariance.first = numerator / chi.Quantile(p);
    confidenceIntervalForVariance.second = numerator / chi.Quantile(1.0 - p); // TODO: redo with Quantile1m
    return true;
}

bool NormalRand::FitMeanBayes(const std::vector<double> &sample, NormalRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double mu0 = priorDistribution.GetLocation();
    double tau0 = priorDistribution.GetPrecision();
    double tau = GetPrecision();
    double numerator = sampleSum(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    priorDistribution.SetLocation(numerator / denominator);
    priorDistribution.SetVariance(1.0 / denominator);
    SetLocation(priorDistribution.Mean());
    return true;
}

bool NormalRand::FitVarianceBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double newAlpha = alpha + 0.5 * n;
    double newBeta = beta + 0.5 * n * sampleVariance(sample, mu);
    priorDistribution.SetParameters(newAlpha, newBeta);
    SetVariance(priorDistribution.Mean());
    return true;
}

bool NormalRand::FitBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    double mu0 = priorDistribution.GetLocation();
    double lambda = priorDistribution.GetPrecision();
    double sum = sampleSum(sample), average = sum / n;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = sampleVariance(sample, average);
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    priorDistribution.SetParameters(newMu0, newLambda, newAlpha, newBeta);
    DoublePair mean = priorDistribution.Mean();
    SetLocation(mean.first);
    SetVariance(mean.second);
    return true;
}
