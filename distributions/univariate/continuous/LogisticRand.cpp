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
    s = (scale <= 0.0) ? 1.0 : scale;
    logS = std::log(s);
}

double LogisticRand::f(const double & x) const
{
    double numerator = std::exp((mu - x) / s);
    double denominator = (1 + numerator);
    denominator *= denominator;
    denominator *= s;
    return numerator / denominator;
}

double LogisticRand::logf(const double & x) const
{
    double x0 = (mu - x) / s;
    double y = std::exp(x0);
    y = 2 * std::log1p(y);
    y += logS;
    return x0 - y;
}

double LogisticRand::F(const double & x) const
{
    double expX = std::exp((mu - x) / s);
    return 1.0 / (1 + expX);
}

double LogisticRand::S(const double & x) const
{
    double expX = std::exp((mu - x) / s);
    return expX / (1 + expX);
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

std::complex<double> LogisticRand::CFImpl(double t) const
{
    double pist = M_PI * s * t;
    std::complex<double> y(0.0, t * mu);
    y = std::exp(y);
    y *= pist;
    y /= std::sinh(pist);
    return y;
}

double LogisticRand::Entropy() const
{
    return 2 + logS;
}

double LogisticRand::quantileImpl(double p) const
{
    return mu - s * (std::log1p(-p) - std::log(p));
}

double LogisticRand::quantileImpl1m(double p) const
{
    return mu - s * (std::log(p) - std::log1p(-p));
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
void LogisticRand::FitLocationMM(const std::vector<double> &sample)
{
    SetLocation(sampleMean(sample));
}

void LogisticRand::FitScaleMM(const std::vector<double> &sample)
{
    double var = sampleVariance(sample, mu);
    SetScale(std::sqrt(3 * var) / M_PI);
}

void LogisticRand::FitMM(const std::vector<double> &sample)
{
    FitLocationMM(sample);
    FitScaleMM(sample);
}

/// Maximum-likelihood
void LogisticRand::FitLocationMLE(const std::vector<double> &sample)
{
    double nHalf = 0.5 * sample.size();
    double root = 0;
    if (!RandMath::findRoot([this, sample, nHalf](double m)
    {
        double f1 = 0, f2 = 0;
        for (const double & x : sample)
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
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetLocation(root);
}
