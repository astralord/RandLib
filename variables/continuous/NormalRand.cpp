#include "NormalRand.h"

unsigned long NormalRand::stairWidth[128] = {0};
double NormalRand::horizontalCoeffs[128] = {0};
double NormalRand::stairHeight[128] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var)
{
    setMean(mean);
    setVar(var);
}

bool NormalRand::setupTables()
{
    const double m1 = 2147483648.0;
    double dn = 3.442619855899, tn = dn, vn = 9.91256303526217e-3;

    stairHeight[0] = 1.;
    stairHeight[127] = std::exp(-.5 * dn * dn);

    double q = vn / stairHeight[127];

    stairWidth[0] = (dn / q) * m1;
    stairWidth[1] = 0;

    horizontalCoeffs[0] = q / m1;
    horizontalCoeffs[127] = dn / m1;

    for (size_t i = 126; i >= 1; --i)
    {
        dn = -std::log(vn / dn + stairHeight[i + 1]);
        dn = std::sqrt(dn + dn);
        stairWidth[i + 1] = (dn / tn) * m1;
        tn = dn;
        stairHeight[i] = std::exp(-.5 * dn * dn);
        horizontalCoeffs[i] = dn / m1;
    }
    return true;
}

void NormalRand::setMean(double mean)
{
    mu = mean;
}

void NormalRand::setSigma(double rootVar)
{
    sigma = std::max(rootVar, MIN_POSITIVE);
    sigmaSqrt2Inv = M_SQRT1_2 / sigma;
}

double NormalRand::f(double x) const
{
    double y = x - mu; /// x - mu
    y *= sigmaSqrt2Inv; /// (x - mu) / (sigma * sqrt(2))
    y *= y; /// (((x - mu) / sigma) ^ 2) / 2
    y = std::exp(-y); /// exp((((x - mu) / sigma) ^ 2) / 2)
    return M_1_SQRTPI * sigmaSqrt2Inv * y; /// exp((((x - mu) / sigma) ^ 2) / 2) / (sigma * sqrt(2pi))
}

double NormalRand::F(double x) const
{
    double y = x - mu; /// x - mu
    y *= sigmaSqrt2Inv; /// (x - mu) / (sigma * sqrt(2))
    y = std::erf(y); /// erf((x - mu) / (sigma * sqrt(2)))
    ++y; /// 1 + erf((x - mu) / (sigma * sqrt(2)))
    return .5 * y; /// (1 + erf((x - mu) / (sigma * sqrt(2)))) / 2
}

double NormalRand::variate()
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
        long hz = BasicRandGenerator::getRand();
        unsigned long stripId = hz & 127;
        double x = hz * horizontalCoeffs[stripId];
        if ((unsigned)std::abs(hz) < stairWidth[stripId])
            return x;

        if (stripId == 0) /// handle the base strip
        {
            static double z = -1;

            if (z >= 0) /// we don't have to generate another exponential variable as we already have one
            {
                x = ExponentialRand::standardVariate() * 0.2904764; /// 1.0 / 3.44262
                z -= 0.5 * x * x;
            }

            if (z < 0) /// if previous generation wasn't successful
            {
                double y;
                do {
                    x = ExponentialRand::standardVariate() * 0.2904764; /// 1.0 / 3.44262
                    y = ExponentialRand::standardVariate();
                    z = y - 0.5 * x * x; /// we storage this value as after acceptance it becomes exponentially distributed
                } while (z <= 0);
            }

            //TODO: maybe we need more accuracy?
            x += 3.44262; /// + start of the right tail
            return (hz > 0) ? x : -x;
        }

        /// handle the wedges of other strips
        if (UniformRand::variate(stairHeight[stripId], stairHeight[stripId - 1]) < std::exp(-.5 * x * x))
            return x;
    } while (++iter <= 1e9); /// one billion should be enough
    return 0; /// fail due to some error
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
