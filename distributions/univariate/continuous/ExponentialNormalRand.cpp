#include "ExponentialNormalRand.h"

ExponentialNormalRand::ExponentialNormalRand(double location, double variance, double rate)
{
    SetParameters(location, variance, rate);
}

std::string ExponentialNormalRand::Name() const
{
    return "Exponential-Normal(" + toStringWithPrecision(GetLocation()) + ", "
                                 + toStringWithPrecision(X.Variance()) + ", "
                                 + toStringWithPrecision(GetRate()) + ")";
}

void ExponentialNormalRand::SetParameters(double location, double variance, double rate)
{
    X.SetLocation(location);
    X.SetVariance(variance);
    Y.SetRate(rate);

    double mu = X.GetLocation();
    double sigma = X.GetScale();
    double lambda = Y.GetRate();
    double var = sigma * sigma;
    a = 0.5 * lambda * var;
    c = mu + a;
    a += c;
    b = M_SQRT1_2 / sigma;
    v = lambda * sigma;
}

double ExponentialNormalRand::f(double x) const
{
    double lambda = Y.GetRate();
    double y = a - x;
    y *= b;
    y = std::erfc(y);
    y *= 0.5 * lambda;
    double exponent = c - x;
    exponent *= lambda;
    exponent = std::exp(exponent);
    return y * exponent;
}

double ExponentialNormalRand::F(double x) const
{
    double u = Y.GetRate() * (x - X.GetLocation());
    double y = X.F(x);
    double exponent = -u + 0.5 * v * v;
    exponent = std::exp(exponent);
    exponent *= X.F(x - v * X.GetScale());
    return y - exponent;
}

double ExponentialNormalRand::Variate() const
{
    return X.Variate() + Y.Variate();
}

double ExponentialNormalRand::Variate(double location, double rootVar, double rate)
{
    return NormalRand::Variate(location, rootVar) + ExponentialRand::Variate(rate);
}

double ExponentialNormalRand::Mean() const
{
    return X.Mean() + Y.Mean();
}

double ExponentialNormalRand::Variance() const
{
    return X.Variance() + Y.Variance();
}

std::complex<double> ExponentialNormalRand::CF(double t) const
{
    return X.CF(t) * Y.CF(t);
}

double ExponentialNormalRand::Skewness() const
{
    double sigma = X.GetScale();
    double lambda = Y.GetRate();
    double tmp = 1.0 / (sigma * lambda);
    double tmpSq = tmp * tmp;
    double y = 1.0 + tmpSq;
    y = y * y * y;
    y = std::sqrt(y);
    y = tmpSq * tmp / y;
    return y + y;
}

double ExponentialNormalRand::ExcessKurtosis() const
{
    double sigma = X.GetScale();
    double lambda = Y.GetRate();
    double tmp = 1.0 / (sigma * lambda);
    tmp *= tmp;
    double numerator = 1.0 + 2.0 * tmp + 3.0 * tmp * tmp;
    double denominator = 1.0 + tmp;
    denominator *= denominator;
    double y = numerator / denominator - 1.0;
    return 3.0 * y;
}

