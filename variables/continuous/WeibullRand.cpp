#include "WeibullRand.h"

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
    l = std::max(scale, MIN_POSITIVE);
    k = std::max(shape, MIN_POSITIVE);
    lInv = 1.0 / l;
    kInv = 1.0 / k;
}

double WeibullRand::f(double x) const
{
    if (x <= 0)
        return 0;
    double xAdj = x * lInv;
    double xAdjPow = std::pow(xAdj, k - 1);
    return k * lInv * xAdjPow * std::exp(-xAdj * xAdjPow);
}

double WeibullRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return 1 - std::exp(-std::pow(x * lInv, k));
}

double WeibullRand::variate() const
{
    return l * std::pow(ExponentialRand::standardVariate(), kInv);
}

double WeibullRand::E() const
{
    return l * std::tgamma(1 + kInv);
}

double WeibullRand::Var() const
{
    double res = std::tgamma(1 + kInv);
    res *= res;
    res += std::tgamma(1 + kInv + kInv);
    return l * l * res;
}

double WeibullRand::Median() const
{
    return l * std::pow(M_LN2, kInv);
}

double WeibullRand::Mode() const
{
    if (k < 1)
        return INFINITY;
    return l * std::pow(1 - kInv, kInv);
}

double WeibullRand::Skewness() const
{
    double mu = E();
    double var = Var();
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
    double mu = E();
    double var = Var();
    double sigma = std::sqrt(var);
    double skewness = Skewness();
    double numerator = l * l;
    numerator *= numerator;
    numerator *= std::tgamma(1 + 4.0 * kInv);
    numerator -= 4 * skewness * var * sigma * mu;
    double mu2 = mu * mu;
    numerator -= 6 * mu2 * var;
    numerator -= mu2 * mu2;
    double kurtosis = numerator / (var * var);
    return kurtosis - 3;
}
