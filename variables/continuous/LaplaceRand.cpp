#include "LaplaceRand.h"

LaplaceRand::LaplaceRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

void LaplaceRand::setLocation(double location)
{
    mu = location;
}

void LaplaceRand::setScale(double scale)
{
    b = std::max(scale, MIN_POSITIVE);
    bInv = 1.0 / b;
    X.setRate(bInv);
}

double LaplaceRand::f(double x) const
{
    double y = -std::fabs(x - mu);
    y *= bInv;
    y = std::exp(y);
    y *= bInv;
    return .5 * y;
}

double LaplaceRand::F(double x) const
{
    double y = x - mu;
    y *= bInv;
    if (x < mu)
        return .5 * std::exp(y);
    y = -.5 * std::exp(-y);
    return y + 1;
}

double LaplaceRand::variate()
{
    double e = X.variate();
    return mu + (((signed)BasicRandGenerator::getRand() > 0) ? e : -e);
}

bool LaplaceRand::fitToData(const QVector<double> &sample)
{
    if (sample.size() == 0)
        return false;

    /// Calculate median
    /// we use root-finding algorithm for median search
    /// it is better to use median-for-median algorithm
    double median = 0.0;
    double min = sample[0], max = min;
    for (double var : sample) {
        min = std::min(var, min);
        max = std::max(var, max);
    }
    RandMath::findRoot([sample] (double med)
    {
        double x = 0;
        for (double var : sample) {
            if (var > med)
                ++x;
            else if (var < med)
                --x;
        }
        return x;
    },
    min, max, median
    );


    /// Calculate scale
    double deviation = 0.0;
    for (double var : sample) {
        deviation += std::fabs(var - median);
    }
    deviation /= sample.size();

    setLocation(median);
    setScale(deviation);
    return true;
}

