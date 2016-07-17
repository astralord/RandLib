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
    kInv = 1.0 / k;
    kSq = k * k;
    pdfCoef = 1.0 / (sigma * (k + kInv));
    cdfCoef = 1.0 / (1 + kSq);
    setLocation((1.0 - kSq) * sigma / k);
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
    int n = sample.size();
    if (n <= 0)
        return false;

    /// Calculate median (considering asymmetry)
    /// we use root-finding algorithm for median search
    double median = 0.0;
    double minVar = sample[0], maxVar = minVar;
    for (double var : sample) {
        minVar = std::min(var, minVar);
        maxVar = std::max(var, maxVar);
        median += var;
    }
    median /= n; /// sample mean

    if (!RandMath::findRoot([this, sample] (double med)
    {
        double y = 0.0;
        for (double x : sample) {
            if (x > med)
                y -= kSq;
            else if (x < med)
                ++y;
        }
        return y;
    },
    minVar, maxVar, median
    ))
        return false;

    setShift(median);
    return true;
}

bool LaplaceRand::fitScaleMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0)
        return false;

    double deviation = 0.0;
    for (double x : sample) {
        if (x > m)
            deviation += kSq * (x - m);
        else
            deviation -= (x - m);
    }
    deviation /= (k * n);

    setScale(deviation);
    return true;
}
    
bool LaplaceRand::fitLocationAndScaleMLE(const std::vector<double> &sample)
{
    return fitLocationMLE(sample) ? fitScaleMLE(sample) : false;
}

bool LaplaceRand::fitLocationMM(const std::vector<double> &sample)
{
    double y = RandMath::sampleMean(sample);
    setShift(y - sigma * (1.0 / k - k));
    return true;
}

bool LaplaceRand::fitScaleMM(const std::vector<double> &sample)
{
    if (k == 1) {
        /// can't derive scale from mean
        double var = RandMath::sampleVariance(sample, m);
        setScale(std::sqrt(0.5 * var));
    }
    else {
        double y = RandMath::sampleMean(sample);
        setScale((y - m) * k / (1.0 - kSq));
    }
    return true;
}

bool LaplaceRand::fitLocationAndScaleMM(const std::vector<double> &sample)
{
    if (!fitLocationMM(sample))
        return false;
    double var = RandMath::sampleVariance(sample, m);
    setScale(std::sqrt(0.5 * var));
    return true;
}
