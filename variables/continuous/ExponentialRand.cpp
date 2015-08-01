#include "ExponentialRand.h"

unsigned long ExponentialRand::stairWidth[256] = {0};
double ExponentialRand::horizontalCoeffs[256] = {0};
double ExponentialRand::stairHeight[256] = {0};
bool ExponentialRand::dummy = ExponentialRand::setupTables();

ExponentialRand::ExponentialRand(double rate)
{
    setRate(rate);
}

bool ExponentialRand::setupTables()
{
    /// Set up ziggurat tables
    const double m2 = 4294967296.;
    double de = 7.697117470131487, te = de, ve = 3.949659822581572e-3;

    stairHeight[0] = 1.;
    stairHeight[255] = std::exp(-de);

    double q = ve / stairHeight[255];
    stairWidth[0] = (de / q) * m2;
    stairWidth[1] = 0;

    horizontalCoeffs[0] = q / m2;
    horizontalCoeffs[255] = de / m2;

    for (unsigned i = 254; i >= 1; --i)
    {
       de = -std::log(ve / de + stairHeight[i + 1]);
       stairWidth[i + 1] = (de / te) * m2;
       te = de;
       stairHeight[i] = std::exp(-de);
       horizontalCoeffs[i] = de / m2;
    }
    return true;
}

void ExponentialRand::setRate(double rate)
{
    l = std::max(rate, MIN_POSITIVE);
    beta = 1.0 / l;
}

double ExponentialRand::f(double x) const
{
    return (x > 0) ? l * std::exp(-l * x) : 0;
}

double ExponentialRand::F(double x) const
{
    return (x > 0) ? 1 - std::exp(-l * x) : 0;
}

double ExponentialRand::variate()
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
        unsigned long B = BasicRandGenerator::getRand();
        unsigned long stairId = B & 255; /// get Uniform[0, 255]
        double x = B * horizontalCoeffs[stairId]; /// get horizontal coordinate

        if (B < stairWidth[stairId]) /// if we are under the stair - accept
            return x;

        if (stairId == 0) /// if we catch the zero stair - launch fallback for tail
            return (x1 + standardVariate());

        if (UniformRand::variate(stairHeight[stairId], stairHeight[stairId - 1]) < std::exp(-x)) /// if we are under the curve - accept
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
