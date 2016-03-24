#include "ExponentialNormalRand.h"

ExponentialNormalRand::ExponentialNormalRand(double location, double variance, double rate)
{
    setParameters(location, variance, rate);
}

std::string ExponentialNormalRand::name()
{
    return "Exponential-Normal(" + toStringWithPrecision(getLocation()) + ", "
                                 + toStringWithPrecision(X.getVariance()) + ", "
                                 + toStringWithPrecision(getRate()) + ")";
}

void ExponentialNormalRand::setParameters(double location, double variance, double rate)
{
    X.setLocation(location);
    X.setVariance(variance);
    Y.setRate(rate);

    double mu = X.getLocation();
    double sigma = X.getScale();
    double lambda = Y.getRate();
    double var = sigma * sigma;
    a = 0.5 * lambda * var;
    c = mu + a;
    a += c;
    b = M_SQRT1_2 / sigma;
    v = lambda * sigma;
}

double ExponentialNormalRand::f(double x) const
{
    double lambda = Y.getRate();
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
    double u = Y.getRate() * (x - X.getLocation());
    double y = X.F(x);
    double exponent = -u + 0.5 * v * v;
    exponent = std::exp(exponent);
    exponent *= X.F(x - v * X.getScale());
    return y - exponent;
}

double ExponentialNormalRand::variate() const
{
    return X.variate() + Y.variate();
}

double ExponentialNormalRand::variate(double location, double rootVar, double rate)
{
    return NormalRand::variate(location, rootVar) + ExponentialRand::variate(rate);
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
    double sigma = X.getScale();
    double lambda = Y.getRate();
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
    double sigma = X.getScale();
    double lambda = Y.getRate();
    double tmp = 1.0 / (sigma * lambda);
    tmp *= tmp;
    double numerator = 1.0 + 2.0 * tmp + 3.0 * tmp * tmp;
    double denominator = 1.0 + tmp;
    denominator *= denominator;
    double y = numerator / denominator - 1.0;
    return 3.0 * y;
}

