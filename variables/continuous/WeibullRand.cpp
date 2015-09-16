#include "WeibullRand.h"
#include "ExponentialRand.h"

WeibullRand::WeibullRand(double scale, double shape)
{
    setParameters(scale, shape);
}

std::string WeibullRand::name()
{
    return "Weibull(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void WeibullRand::setParameters(double scale, double shape)
{
    lambda = scale;
    if (lambda <= 0)
        lambda = MIN_POSITIVE;
    lambdaInv = 1.0 / lambda;
    
    k = shape;
    if (k <= 0)
        k = MIN_POSITIVE;
    kInv = 1.0 / k;
}

double WeibullRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x * lambdaInv;
    double xAdjPow = std::pow(xAdj, k - 1);
    return k * lambdaInv * xAdjPow * std::exp(-xAdj * xAdjPow);
}

double WeibullRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return 1 - std::exp(-std::pow(x * lambdaInv, k));
}

double WeibullRand::variate() const
{
    return lambda * std::pow(ExponentialRand::standardVariate(), kInv);
}

double WeibullRand::Mean() const
{
    return lambda * std::tgamma(1 + kInv);
}

double WeibullRand::Variance() const
{
    double res = std::tgamma(1 + kInv);
    res *= res;
    res += std::tgamma(1 + kInv + kInv);
    return lambda * lambda * res;
}

double WeibullRand::Median() const
{
    return lambda * std::pow(M_LN2, kInv);
}

double WeibullRand::Mode() const
{
    if (k < 1)
        return INFINITY;
    return lambda * std::pow(1 - kInv, kInv);
}

double WeibullRand::Skewness() const
{
    double mu = Mean();
    double var = Variance();
    double sigma = std::sqrt(var);
    double numerator = std::tgamma(1 + 3.0 * kInv);
    numerator *= lambda * lambda * lambda;
    numerator -= 3 * mu * var;
    numerator -= mu * mu * mu;
    double denominator = var * sigma;
    return numerator / denominator;
}

double WeibullRand::ExcessKurtosis() const
{
    double mu = Mean();
    double var = Variance();
    double sigma = std::sqrt(var);
    double skewness = Skewness();
    double numerator = lambda * lambda;
    numerator *= numerator;
    numerator *= std::tgamma(1 + 4.0 * kInv);
    numerator -= 4 * skewness * var * sigma * mu;
    double mu2 = mu * mu;
    numerator -= 6 * mu2 * var;
    numerator -= mu2 * mu2;
    double kurtosis = numerator / (var * var);
    return kurtosis - 3;
}
