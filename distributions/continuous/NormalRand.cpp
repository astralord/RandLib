#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"

double NormalRand::stairWidth[257] = {0};
double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var) : StableRand(2.0, 0.0, 1.0, mean)
{
    setVariance(var);
}

std::string NormalRand::name()
{
    return "Normal(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getVar()) + ")";
}

bool NormalRand::setupTables()
{
    static constexpr long double A = 4.92867323399e-3l; /// area under rectangle

    /// coordinates of the implicit rectangle in base layer
    stairHeight[0] = std::exp(-.5 * x1 * x1);
    stairWidth[0] = A / stairHeight[0];
    /// implicit value for the top layer
    stairWidth[256] = 0;

    for (unsigned i = 1; i <= 255; ++i)
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

double NormalRand::getVar() const
{
    return sigma * sigma;
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
    return mu + sigma * standardVariate();
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
            double y;

            if (z >= 0) /// we don't have to generate another exponential variable as we already have one
            {
                y = ExponentialRand::standardVariate();
                z = y - 0.5 * z * z;
            }

            if (z < 0) /// if previous generation wasn't successful
            {
                do {
                    x = ExponentialRand::variate(x1);
                    y = ExponentialRand::standardVariate();
                    z = y - 0.5 * x * x; /// we storage this value as after acceptance it becomes exponentially distributed
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

double NormalRand::Mean() const
{
    return mu;
}

double NormalRand::Variance() const
{
    return sigma * sigma;
}

std::complex<double> NormalRand::CF(double t) const
{
    double sigmaT = sigma * t;
    return std::exp(std::complex<double>(0.5 * sigmaT * sigmaT, mu * t));
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
    if (p == 0.0)
        return -INFINITY;
    if (p == 1.0)
        return INFINITY;
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
    return mu + sigma * standardQuantile(p);
}

double NormalRand::Median() const
{
    return mu;
}

double NormalRand::Mode() const
{
    return mu;
}

double NormalRand::Skewness() const
{
    return 0.0;
}

double NormalRand::ExcessKurtosis() const
{
    return 0.0;
}

double NormalRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return (n & 1) ? std::pow(sigma, n) * RandMath::doubleFactorial(n - 1) : 0;
}

bool NormalRand::fitMean_MLE(const QVector<double> &sample)
{
    setLocation(RandMath::sampleMean(sample));
    return true;
}

bool NormalRand::fitVariance_MLE(const QVector<double> &sample)
{
    setVariance(RandMath::sampleVariance(sample));
    return true;
}

bool NormalRand::fit_MLE(const QVector<double> &sample)
{
    return fitMean_MLE(sample) ? fitVariance_MLE(sample) : false;
}

bool NormalRand::fitMean_MM(const QVector<double> &sample)
{
    return fitMean_MLE(sample);
}

bool NormalRand::fitVariance_MM(const QVector<double> &sample)
{
    return fitVariance_MLE(sample);
}

bool NormalRand::fit_MM(const QVector<double> &sample)
{
    return fit_MLE(sample);
}

bool NormalRand::fitMean_UMVU(const QVector<double> &sample)
{
    return fitMean_MLE(sample);
}

bool NormalRand::fitVariance_UMVU(const QVector<double> &sample)
{
    int n = sample.size();
    if (n <= 1)
        return false;
    double s = RandMath::sampleVariance(sample, mu);
    setVariance(n * s / (n - 1));
    return true;
}

bool NormalRand::fit_UMVU(const QVector<double> &sample)
{
    return fitMean_MLE(sample) ? fitVariance_UMVU(sample) : false;
}
