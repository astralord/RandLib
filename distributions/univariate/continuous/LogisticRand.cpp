#include "LogisticRand.h"

LogisticRand::LogisticRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string LogisticRand::name() const
{
    return "Logistic(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void LogisticRand::setLocation(double location)
{
    mu = location;
}

void LogisticRand::setScale(double scale)
{
    s = scale;
    if (s <= 0)
        s = 1.0;
}

double LogisticRand::f(double x) const
{
    double numerator = std::exp((mu - x) / s);
    double denominator = (1 + numerator);
    denominator *= denominator;
    denominator *= s;
    return numerator / denominator;
}

double LogisticRand::F(double x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

double LogisticRand::variate() const
{
    /// there can be used rejection method from Laplace or Cauchy (Luc Devroye, p. 471)
    /// or ziggurat
    return mu + s * std::log(1.0 / UniformRand::standardVariate() - 1);
}

double LogisticRand::Variance() const
{
    double sPi = s * M_PI;
    return sPi * sPi / 3;
}

std::complex<double> LogisticRand::CF(double t) const
{
    if (t == 0)
        return 1;
    double pist = M_PI * s * t;
    std::complex<double> y(0.0, t * mu);
    y = std::exp(y);
    y *= pist;
    y /= std::sinh(pist);
    return y;
}

double LogisticRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    return mu - s * std::log(1.0 / p - 1);
}

double LogisticRand::Median() const
{
    return mu;
}

double LogisticRand::Mode() const
{
    return mu;
}

double LogisticRand::Skewness() const
{
    return 0;
}

double LogisticRand::ExcessKurtosis() const
{
    return 1.2;
}

/// Method of moments
bool LogisticRand::fitLocationMM(const std::vector<double> &sample)
{
    setLocation(RandMath::sampleMean(sample));
    return true;
}

bool LogisticRand::fitScaleMM(const std::vector<double> &sample)
{
    double var = RandMath::sampleVariance(sample, mu);
    setScale(std::sqrt(3 * var) / M_PI);
    return true;
}

bool LogisticRand::fitMM(const std::vector<double> &sample)
{
    return fitLocationMM(sample) ? fitScaleMM(sample) : false;
}

bool LogisticRand::fitLocationMLE(const std::vector<double> &sample)
{
    double nHalf = 0.5 * sample.size();
    if (nHalf <= 0)
        return false;
    double root = 0;
    if (!RandMath::findRoot([this, sample, nHalf](double m)
    {
        double f1 = 0, f2 = 0;
        for (double x : sample)
        {
            double aux = std::exp((m - x) / s);
            double denom = 1.0 + aux;
            f1 += 1.0 / denom;
            denom *= denom;
            f2 -= aux / denom;
        }
        f1 -= nHalf;
        return DoublePair(f1, f2);
    }, root))
        return false;

    setLocation(root);
    return true;
}
