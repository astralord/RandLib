#include "LaplaceRand.h"
#include "../discrete/BernoulliRand.h"

LaplaceRand::LaplaceRand(double shift, double scale, double asymmetry)
    : GeometricStableRand(2.0, 0.0)
{
    setShift(shift);
    setScale(scale);
    setAsymmetry(asymmetry);
}

std::string LaplaceRand::name() const
{
    return "Laplace(" + toStringWithPrecision(getShift()) + ", "
                      + toStringWithPrecision(getScale()) + ", "
            + toStringWithPrecision(getAsymmetry()) + ")";
}

void LaplaceRand::setShift(double shift)
{
    m = shift;
}

void LaplaceRand::setAsymmetry(double asymmetry)
{
    k = asymmetry;
    if (k <= 0)
        k = 1.0;
    setLocation((1.0 - k * k) * sigma / k);
}

double LaplaceRand::f(double x) const
{
    return pdfLaplace(x - m);
}

double LaplaceRand::F(double x) const
{
    return cdfLaplace(x - m);
}

double LaplaceRand::variate() const
{
    return (k == 1) ? LaplaceRand::variate(m, sigma) : LaplaceRand::variate(m, sigma, k);
}

double LaplaceRand::variate(double location, double scale)
{
    bool sign = BernoulliRand::standardVariate();
    double W = scale * ExponentialRand::standardVariate();
    return location + (sign ? W : -W);
}

double LaplaceRand::variate(double location, double scale, double asymmetry)
{
    double x = ExponentialRand::standardVariate() / asymmetry;
    double y = ExponentialRand::standardVariate() * asymmetry;;
    return location + scale * (x - y);
}

void LaplaceRand::sample(std::vector<double> &outputData) const
{
    if (k == 1) {
        for (double & var : outputData)
            var = LaplaceRand::variate(m, sigma);
    }
    else {
        for (double & var : outputData)
            var = LaplaceRand::variate(m, sigma, k);
    }
}

double LaplaceRand::Mean() const
{
    return m + GeometricStableRand::Mean();
}

std::complex<double> LaplaceRand::CF(double t) const
{
    if (t == 0)
        return 1;
    double bt = sigma * t;
    double btSq = bt * bt;
    double denominator = (1 + kSq * btSq) * (1 + btSq / kSq);
    std::complex<double> y(std::cos(m * t), std::sin(m * t));
    std::complex<double> x(1, -k * bt), z(1, bt * kInv);
    return x * y * z / denominator;
}

double LaplaceRand::Median() const
{
    return m + GeometricStableRand::Median();
}

double LaplaceRand::Mode() const
{
    return m;
}

double LaplaceRand::Entropy() const
{
    double y = 1 + kSq;
    y *= kInv * sigma;
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

    setShift(median);
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
    setShift(RandMath::sampleMean(sample));
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
