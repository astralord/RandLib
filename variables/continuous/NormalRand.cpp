#include "NormalRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"

double NormalRand::stairWidth[257] = {0};
double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var)
{
    setMean(mean);
    setVar(var);
}

std::string NormalRand::name()
{
    return "Normal(" + toStringWithPrecision(getMean()) + ", " + toStringWithPrecision(getVar()) + ")";
}

bool NormalRand::setupTables()
{
    constexpr double A = 4.92867323399e-3; /// area under rectangle

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

void NormalRand::setMean(double mean)
{
    mu = mean;
}

void NormalRand::setSigma(double rootVar)
{
    sigma = rootVar;
    if (sigma <= 0)
        sigma = 1.0;
    sigmaSqrt2Inv = M_SQRT1_2 / sigma;
}

void NormalRand::setVar(double var)
{
    setSigma(std::sqrt(std::max(var, 0.0)));
}

double NormalRand::f(double x) const
{
    double y = x - mu;
    y *= sigmaSqrt2Inv;
    y *= y;
    y = std::exp(-y);
    return M_1_SQRTPI * sigmaSqrt2Inv * y;
}

double NormalRand::F(double x) const
{
    double y = x - mu;
    y *= sigmaSqrt2Inv;
    y = std::erf(y);
    ++y;
    return .5 * y;
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
    int iter = 0;
    do {
        unsigned long B = RandGenerator::variate();
        unsigned long stairId = B & 255;
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
    return 0; /// fail due to some error
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

double NormalRand::Quantile(double p) const
{
    if (p < 0.5)
        return -Quantile(1 - p);
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    if (p == 1)
        return INFINITY;
    double t = -std::log(p);
    t = std::sqrt(t + t);
    static constexpr double c[] = {2.515517, 0.802853, 0.010328};
    static constexpr double d[] = {1.432788, 0.189269, 0.001308};
    double numerator = (c[2] * t + c[1]) * t + c[0];
    double denominator = ((d[2] * t + d[1]) * t + d[0]) * t + 1.0;
    return mu + sigma * (numerator / denominator - t);
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
    return (n & 1) ? std::pow(sigma, n) * RandMath::doubleFactorial(n - 1): 0;
}

bool NormalRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate mu
    double average = 0.0;
    for (double var : sample) {
        average += var;
    }
    average /= sample.size();

    /// Calculate sigma
    double deviation = 0.0;
    for (double var : sample) {
        double currDev = (var - average);
        deviation += currDev * currDev;
    }
    deviation /= std::max(sample.size() - 1, 1);

    setMean(average);
    setSigma(std::sqrt(deviation));
    return true;
}
