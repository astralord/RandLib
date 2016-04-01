#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

LaplaceRand::LaplaceRand(double location, double scale, double asymmetry)
{
    setLocation(location);
    setScale(scale);
    setAsymmetry(asymmetry);
}

std::string LaplaceRand::name()
{
    return "Laplace(" + toStringWithPrecision(getLocation()) + ", "
                      + toStringWithPrecision(getScale()) + ", "
                      + toStringWithPrecision(getAsymmetry()) + ")";
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
    pdfCoef = bInv / (k + kInv);
}

void LaplaceRand::setAsymmetry(double asymmetry)
{
    k = asymmetry;
    if (k <= 0)
        k = 1.0;
    kInv = 1.0 / k;
    kSq = k * k;
    pdfCoef = bInv / (k + kInv);
    cdfCoef = 1.0 / (1 + kSq);
}

double LaplaceRand::f(double x) const
{
    double y = bInv * (x - mu);
    y *= (x < mu) ? kInv : -k;
    y = std::exp(y);
    return pdfCoef * y;
}

double LaplaceRand::F(double x) const
{
    double y = bInv * (x - mu);
    if (x < mu) {
        y *= kInv;
        y = std::exp(y);
        return 1.0 - cdfCoef * y;
    }
    y *= -k;
    y = std::exp(y);
    return kSq * cdfCoef * y;
}

double LaplaceRand::variate() const
{
    return LaplaceRand::variate(mu, b, k);
}

double LaplaceRand::variate(double location, double scale, double asymmetry)
{
    double x = ExponentialRand::standardVariate() / asymmetry;
    double y = ExponentialRand::standardVariate() * asymmetry;;
    return location + scale * (x - y);
}

double LaplaceRand::Mean() const
{
    return mu + (1 - kSq) * b * kInv;
}

double LaplaceRand::Variance() const
{
    double y = b * b / kSq;
    return (1.0 + kSq * kSq) * y;
}

std::complex<double> LaplaceRand::CF(double t) const
{
    double bt = b * t;
    double btSq = bt * bt;
    double denominator = (1 + kSq * btSq) * (1 + btSq / kSq);
    std::complex<double> y(std::cos(mu * t), std::sin(mu * t));
    std::complex<double> x(1, -k * bt), z(1, bt * kInv);
    return x * y * z / denominator;
}

double LaplaceRand::Median() const
{
    double y = 0.5 * (1.0 / kSq + 1.0);
    y = std::log(y);
    return mu + b * k * y;
}

double LaplaceRand::Mode() const
{
    return mu;
}

double LaplaceRand::Skewness() const
{
    double k4 = kSq * kSq;
    double k6 = k4 * kSq;
    double z = (k4 + 1);
    z *= std::sqrt(z);
    double y = 2 * (1 - k6);
    return y / z;
}

double LaplaceRand::ExcessKurtosis() const
{
    if (k == 1)
        return 3.0;
    return NAN; // TODO!!
}

double LaplaceRand::Entropy() const
{
    double y = 1 + kSq;
    y *= kInv * b;
    return log1p(y);
}

bool LaplaceRand::fitLocationMLE(const std::vector<double> &sample)
{
    if (k != 1)
        return false;
    int n = sample.size();
    if (n <= 0)
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

bool LaplaceRand::fitScaleMLE(const std::vector<double> &sample)
{
    if (k != 1)
        return false;
    int n = sample.size();
    if (n <= 0)
        return false;

    double deviation = 0.0;
    for (double var : sample) {
        deviation += std::fabs(var - mu);
    }
    deviation /= n;

    setScale(deviation);
    return true;
}
    
bool LaplaceRand::fitLocationAndScaleMLE(const std::vector<double> &sample)
{
    return fitLocationMLE(sample) ? fitScaleMLE(sample) : false;
}

bool LaplaceRand::fitLocationMM(const std::vector<double> &sample)
{
    if (k != 1)
        return false;
    setLocation(RandMath::sampleMean(sample));
    return true;
}

bool LaplaceRand::fitScaleMM(const std::vector<double> &sample)
{
    if (k != 1)
        return false;
    double var = RandMath::sampleVariance(sample, mu);
    setScale(std::sqrt(0.5 * var));
    return true;
}

bool LaplaceRand::fitLocationAndScaleMM(const std::vector<double> &sample)
{
    return fitLocationMM(sample) ? fitScaleMM(sample) : false;
}
