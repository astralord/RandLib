#include "LogisticRand.h"

LogisticRand::LogisticRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

std::string LogisticRand::Name() const
{
    return "Logistic(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void LogisticRand::SetLocation(double location)
{
    mu = location;
}

void LogisticRand::SetScale(double scale)
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

double LogisticRand::Variate() const
{
    /// there can be used rejection method from Laplace or Cauchy (Luc Devroye, p. 471)
    /// or ziggurat
    return mu + s * std::log(1.0 / UniformRand::StandardVariate() - 1);
}

double LogisticRand::Mean() const
{
    return mu;
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

double LogisticRand::quantileImpl(double p) const
{
    return mu - s * std::log(1.0 / p - 1);
}

double LogisticRand::quantileImpl1m(double p) const
{
    return mu - s * std::log(p / (1.0 - p));
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
bool LogisticRand::FitLocationMM(const std::vector<double> &sample)
{
    SetLocation(sampleMean(sample));
    return true;
}

bool LogisticRand::FitScaleMM(const std::vector<double> &sample)
{
    double var = sampleVariance(sample, mu);
    SetScale(std::sqrt(3 * var) / M_PI);
    return true;
}

bool LogisticRand::FitMM(const std::vector<double> &sample)
{
    return FitLocationMM(sample) ? FitScaleMM(sample) : false;
}

bool LogisticRand::FitLocationMLE(const std::vector<double> &sample)
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

    SetLocation(root);
    return true;
}
