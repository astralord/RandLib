#include "NormalRand.h"

double NormalRand::stairWidth[257] = {0};
double NormalRand::stairHeight[256] = {0};
const bool NormalRand::dummy = NormalRand::setupTables();

NormalRand::NormalRand(double mean, double var)
{
    setMean(mean);
    setVar(var);
}

void NormalRand::setName()
{
    nameStr = "Normal(" + toStringWithPrecision(getMean()) + ", " + toStringWithPrecision(getVar()) + ")";
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
    setName();
}

void NormalRand::setSigma(double rootVar)
{
    sigma = std::max(rootVar, MIN_POSITIVE);
    sigmaSqrt2Inv = M_SQRT1_2 / sigma;
    setName();
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
        unsigned long B = BasicRandGenerator::variate();
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
