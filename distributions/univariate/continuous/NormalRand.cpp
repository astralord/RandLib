#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "../BasicRandGenerator.h"
#include "GammaRand.h"
#include "StudentTRand.h"


double NormalRand::stairWidth[257] = {0};
double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var)
    : StableRand(2.0, 0.0, 1.0, mean)
{
    setVariance(var);
}

std::string NormalRand::name()
{
    return "Normal(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(Variance()) + ")";
}

void NormalRand::setScale(double scale)
{
    sigma0 = scale;
    if (sigma0 <= 0.0)
        sigma0 = 1.0;
    StableRand::setScale(scale * M_SQRT1_2);
}

bool NormalRand::setupTables()
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

void NormalRand::setVariance(double var)
{
    if (var <= 0)
        var = 1.0;
    setScale(std::sqrt(var));
}

double NormalRand::f(double x) const
{
    return StableRand::pdfNormal(x);
}

double NormalRand::F(double x) const
{
    return StableRand::cdfNormal(x);
}

double NormalRand::variate() const
{
    return mu + sigma0 * standardVariate();
}

double NormalRand::variate(double mean, double rootVar)
{
    return mean + rootVar * standardVariate();
}

double NormalRand::standardVariate()
{
    /// Ziggurat algorithm by George Marsaglia using 256 strips
    int iter = 0;
    do {
        unsigned long long B = RandGenerator::variate();
        int stairId = B & 255;
        double x = UniformRand::standardVariate() * stairWidth[stairId]; /// get horizontal coordinate

        if (x < stairWidth[stairId + 1])
            return ((signed)B > 0) ? x : -x;

        if (stairId == 0) /// handle the base layer
        {
            static double z = -1;

            if (z > 0) /// we don't have to generate another exponential variable as we already have one
            {
                x = ExponentialRand::variate(x1);
                z -= 0.5 * x * x;
            }

            if (z <= 0) /// if previous generation wasn't successful
            {
                do {
                    x = ExponentialRand::variate(x1);
                    z = ExponentialRand::standardVariate() - 0.5 * x * x; /// we storage this value as after acceptance it becomes exponentially distributed
                } while (z <= 0);
            }

            x += x1;
            return ((signed)B > 0) ? x : -x;
        }

        /// handle the wedges of other stairs
        if (UniformRand::variate(stairHeight[stairId - 1], stairHeight[stairId]) < std::exp(-.5 * x * x))
            return ((signed)B > 0) ? x : -x;

    } while (++iter <= 1e9); /// one billion should be enough
    return NAN; /// fail due to some error
}

std::complex<double> NormalRand::CF(double t) const
{
    if (t == 0)
        return std::complex<double>(1, 0);
    double sigma0T = sigma0 * t;
    return std::exp(std::complex<double>(-0.5 * sigma0T * sigma0T, mu * t));
}

double NormalRand::standardQuantile(double p)
{
    /// Abramowitz and Stegun improved approximation
    if (p < 0.5)
        return -standardQuantile(1.0 - p);
    else
        p = 1.0 - p;
    if (p < 0.0 || p > 1.0)
        return NAN;
    long double t = -std::log(p);
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
    return t - numerator / denominator;
}

double NormalRand::Quantile(double p) const
{
    return mu + sigma0 * standardQuantile(p);
}

double NormalRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return (n & 1) ? std::pow(sigma0, n) * RandMath::doubleFactorial(n - 1) : 0;
}

bool NormalRand::fitMeanMLE(const std::vector<double> &sample)
{
    setLocation(RandMath::sampleMean(sample));
    return true;
}

bool NormalRand::fitVarianceMLE(const std::vector<double> &sample)
{
    setVariance(RandMath::sampleVariance(sample));
    return true;
}

bool NormalRand::fitMeanAndVarianceMLE(const std::vector<double> &sample)
{
    return fitMeanMLE(sample) ? fitVarianceMLE(sample) : false;
}

bool NormalRand::fitMeanMM(const std::vector<double> &sample)
{
    return fitMeanMLE(sample);
}

bool NormalRand::fitVarianceMM(const std::vector<double> &sample)
{
    return fitVarianceMLE(sample);
}

bool NormalRand::fitMeanAndVarianceMM(const std::vector<double> &sample)
{
    return fitMeanAndVarianceMLE(sample);
}

bool NormalRand::fitMeanUMVU(const std::vector<double> &sample)
{
    return fitMeanMLE(sample);
}

bool NormalRand::fitVarianceUMVU(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 1)
        return false;
    double s = RandMath::sampleVariance(sample, mu);
    setVariance(n * s / (n - 1));
    return true;
}

bool NormalRand::fitMeanAndVarianceUMVU(const std::vector<double> &sample)
{
    return fitMeanMLE(sample) ? fitVarianceUMVU(sample) : false;
}

bool NormalRand::fitMeanAndVarianceUMVU(const std::vector<double> &sample, DoublePair &confidenceIntervalForMean, DoublePair &confidenceIntervalForVariance, double alpha)
{
    size_t n = sample.size();
    if (n < 2 || alpha <= 0 || alpha > 1)
        return false;
    if (!fitMeanAndVarianceUMVU(sample))
        return false;
    double p = 1.0 - 0.5 * alpha;

    /// calculate confidence interval for mean
    StudentTRand t(n - 1);
    double interval = t.Quantile(p) * sigma0 / std::sqrt(n);
    confidenceIntervalForMean.first = mu - interval;
    confidenceIntervalForMean.second = mu + interval;

    /// calculate confidence interval for variance
    ChiSquaredRand chi(n - 1);
    double numerator = (n - 1) * sigma0 * sigma0;
    confidenceIntervalForVariance.first = numerator / chi.Quantile(p);
    confidenceIntervalForVariance.second = numerator / chi.Quantile(1.0 - p);
    return true;
}

bool NormalRand::fitMeanBayes(const std::vector<double> &sample, NormalRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double mu0 = priorDistribution.getLocation();
    double tau0 = priorDistribution.getPrecision();
    double tau = getPrecision();
    double numerator = RandMath::sum(sample) * tau + tau0 * mu0;
    double denominator = n * tau + tau0;
    priorDistribution.setLocation(numerator / denominator);
    priorDistribution.setVariance(1.0 / denominator);
    setLocation(priorDistribution.Mean());
    return true;
}

bool NormalRand::fitVarianceBayes(const std::vector<double> &sample, InverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.getShape();
    double beta = priorDistribution.getRate();
    double newAlpha = alpha + 0.5 * n;
    double newBeta = beta + 0.5 * n * RandMath::sampleVariance(sample, mu);
    priorDistribution.setParameters(newAlpha, newBeta);
    setVariance(priorDistribution.Mean());
    return true;
}

bool NormalRand::fitMeanAndVarianceBayes(const std::vector<double> &sample, NormalInverseGammaRand &priorDistribution)
{
    size_t n = sample.size();
    if (n == 0)
        return false;
    double alpha = priorDistribution.getShape();
    double beta = priorDistribution.getRate();
    double mu0 = priorDistribution.getLocation();
    double lambda = priorDistribution.getPrecision();
    double sum = RandMath::sum(sample), average = sum / n;
    double newLambda = lambda + n;
    double newMu0 = (lambda * mu0 + sum) / newLambda;
    double newAlpha = alpha + 0.5 * n;
    double variance = RandMath::sampleVariance(sample, average);
    double aux = mu0 - average;
    double newBeta = beta + 0.5 * n * (variance + lambda / newLambda * aux * aux);
    priorDistribution.setParameters(newMu0, newLambda, newAlpha, newBeta);
    DoublePair mean = priorDistribution.Mean();
    setLocation(mean.first);
    setVariance(mean.second);
    return true;
}
