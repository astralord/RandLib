#include "ExponentialRand.h"

double ExponentialRand::stairWidth[257] = {0};
double ExponentialRand::stairHeight[256] = {0};
bool ExponentialRand::dummy = ExponentialRand::setupTables();

ExponentialRand::ExponentialRand(double rate)
{
    setRate(rate);
}

void ExponentialRand::setName()
{
    nameStr = "Exponential(" + toStringWithPrecision(getRate()) + ")";
}

bool ExponentialRand::setupTables()
{
    /// Set up ziggurat tables
    double constexpr A = 3.9496598225815571993e-3; /// area under rectangle

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
    lambda = std::max(rate, MIN_POSITIVE);
    beta = 1.0 / lambda;
    setName();
}

double ExponentialRand::f(double x) const
{
    return (x > 0) ? lambda * std::exp(-lambda * x) : 0;
}

double ExponentialRand::F(double x) const
{
    return (x > 0) ? 1 - std::exp(-lambda * x) : 0;
}

std::complex<double> ExponentialRand::CF(double t) const
{
    double rate2 = lambda * lambda;
    double denominator = rate2 + t * t;
    return std::complex<double>(rate2 / denominator, lambda * t / denominator);
}

double ExponentialRand::variate() const
{
    return beta * standardVariate();
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
        unsigned long stairId = BasicRandGenerator::variate() & 255;
        double x = UniformRand::standardVariate() * stairWidth[stairId]; /// get horizontal coordinate

        if (x < stairWidth[stairId + 1]) /// if we are under the upper stair - accept
            return x;

        if (stairId == 0) /// if we catch the tail
            return x1 + standardVariate();

        if (UniformRand::variate(stairHeight[stairId - 1], stairHeight[stairId]) < std::exp(-x)) /// if we are under the curve - accept
            return x;

        /// rejection - go back

    } while (++iter <= 1e9); /// one billion should be enough to be sure there is a bug
    return 0; /// fail due to some error
}

bool ExponentialRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate rate
    double average = 0.0;
    for (double var : sample) {
        if (var < 0)
            return false;
        average += var;
    }
    average /= sample.size();
    if (average <= 0)
        return false;

    setRate(1.0 / average);
    return true;
}
