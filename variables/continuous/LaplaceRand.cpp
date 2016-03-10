#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

LaplaceRand::LaplaceRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LaplaceRand::name()
{
    return "Laplace(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LaplaceRand::setLocation(double location)
{
    mu = location;
}

void LaplaceRand::setScale(double scale)
{
    b = scale;
    if (b <= 0)
        b = 1.0;
    bInv = 1.0 / b;
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

double LaplaceRand::variate() const
{
    return LaplaceRand::variate(mu, b);
}

double LaplaceRand::variate(double location, double scale)
{
    double e = scale * ExponentialRand::standardVariate();
    return location + (BernoulliRand::standardVariate() ? -e : e);
}

double LaplaceRand::Mean() const
{
    return mu;
}

double LaplaceRand::Variance() const
{
    return 2 * b * b;
}

std::complex<double> LaplaceRand::CF(double t) const
{
    double bt = b * t;
    double denominator = 1 + bt * bt;
    return std::complex<double>(std::cos(t) / denominator, std::sin(t) / denominator);
}

double LaplaceRand::Median() const
{
    return mu;
}

double LaplaceRand::Mode() const
{
    return mu;
}

double LaplaceRand::Skewness() const
{
    return 0.0;
}

double LaplaceRand::ExcessKurtosis() const
{
    return 3.0;
}

bool LaplaceRand::fitLocation_MLE(const QVector<double> &sample)
{
    int n = sample.size();
    if (n == 0)
        return false;

    /// Calculate median
    /// we use root-finding algorithm for median search
    /// but note, that it can be better to use median-for-median algorithm
    double median = 0.0;
    double minVar = sample[0], maxVar = minVar;
    for (double var : sample) {
        minVar = std::min(var, minVar);
        maxVar = std::max(var, maxVar);
        median += var;
    }
    median /= n;

    if (!RandMath::findRoot([sample] (double med)
    {
        /// sum of sign(x) - derivative of sum of abs(x)
        double x = 0.0;
        for (double var : sample) {
            if (var > med)
                ++x;
            else if (var < med)
                --x;
        }
        return x;
    },
    minVar, maxVar, median
    ))
        return false;

    setLocation(median);
    return true;
}

bool LaplaceRand::fitScale_MLE(const QVector<double> &sample)
{
    int n = sample.size();
    if (n == 0)
        return false;

    double deviation = 0.0;
    for (double var : sample) {
        deviation += std::fabs(var - mu);
    }
    deviation /= n;

    setScale(deviation);
    return true;
}
    
bool LaplaceRand::fit_MLE(const QVector<double> &sample)
{
    return fitLocation_MLE(sample) ? fitScale_MLE(sample) : false;
}

