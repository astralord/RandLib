#include "ExponentialRand.h"
#include "UniformRand.h"
#include "../BasicRandGenerator.h"

double ExponentialRand::stairWidth[257] = {0};
double ExponentialRand::stairHeight[256] = {0};
bool ExponentialRand::dummy = ExponentialRand::setupTables();

ExponentialRand::ExponentialRand(double rate) : GammaRand()
{
    setRate(rate);
}

std::string ExponentialRand::name()
{
    return "Exponential(" + toStringWithPrecision(getRate()) + ")";
}

bool ExponentialRand::setupTables()
{
    /// Set up ziggurat tables
    static constexpr long double A = 3.9496598225815571993e-3l; /// area under rectangle

    /// coordinates of the implicit rectangle in base layer
    stairHeight[0] = std::exp(-x1);
    stairWidth[0] = A / stairHeight[0];
    /// implicit value for the top layer
    stairWidth[256] = 0;

    for (unsigned i = 1; i <= 255; ++i)
    {
        /// such y_i that f(x_{i+1}) = y_i
        stairWidth[i] = -std::log(stairHeight[i - 1]);
        stairHeight[i] = stairHeight[i - 1] + A / stairWidth[i];
    }
    return true;
}

void ExponentialRand::setRate(double rate)
{
    thetaInv = rate;
    if (thetaInv <= 0)
        thetaInv = 1.0;
    setParameters(1, 1.0 / thetaInv);
}

double ExponentialRand::f(double x) const
{
    return (x > 0) ? thetaInv * std::exp(-thetaInv * x) : 0;
}

double ExponentialRand::F(double x) const
{
    return (x > 0) ? 1 - std::exp(-thetaInv * x) : 0;
}

double ExponentialRand::variate() const
{
    return theta * standardVariate();
}

double ExponentialRand::variate(double rate)
{
    return standardVariate() / rate;
}

double ExponentialRand::standardVariate()
{
    /// Ziggurat algorithm
    int iter = 0;
    do {
        int stairId = RandGenerator::variate() & 255;
        double x = UniformRand::standardVariate() * stairWidth[stairId]; /// get horizontal coordinate

        if (x < stairWidth[stairId + 1]) /// if we are under the upper stair - accept
            return x;

        if (stairId == 0) /// if we catch the tail
            return x1 + standardVariate();

        if (UniformRand::variate(stairHeight[stairId - 1], stairHeight[stairId]) < std::exp(-x)) /// if we are under the curve - accept
            return x;

        /// rejection - go back

    } while (++iter <= 1e9); /// one billion should be enough to be sure there is a bug
    return NAN; /// fail due to some error
}

double ExponentialRand::Mean() const
{
    return beta;
}

double ExponentialRand::Variance() const
{
    return beta * beta;
}

std::complex<double> ExponentialRand::CF(double t) const
{
    return std::complex<double>(1.0, -theta * t);
}

double ExponentialRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 1)
        return INFINITY;
    return -theta * std::log(1.0 - p);
}

double ExponentialRand::Median() const
{
    return theta * M_LN2;
}

double ExponentialRand::Entropy() const
{
    return 1 - std::log(lambda);
}

double ExponentialRand::Moment(int n) const
{
    if (n < 0)
        return 0;
    if (n == 0)
        return 1;
    return RandMath::factorial(n) / std::pow(thetaInv, n);
}

bool ExponentialRand::fit_MLE(const QVector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setRate(1.0 / RandMath::sampleMean(sample));
    return true;
}

bool ExponentialRand::fit_MM(const QVector<double> &sample)
{
    return fit_MLE(sample);
}

bool ExponentialRand::fit_UMVU(const QVector<double> &sample)
{   
    int n = sample.size();
    if (n <= 1 || !checkValidity(sample))
        return false;
    double mean = RandMath::sampleMean(sample);
    setRate((n - 1) * mean / n);
    return true;
}
