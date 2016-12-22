#include "WeibullRand.h"
#include "ExponentialRand.h"

WeibullRand::WeibullRand(double scale, double shape)
{
    SetParameters(scale, shape);
}

std::string WeibullRand::Name() const
{
    return "Weibull(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void WeibullRand::SetParameters(double scale, double shape)
{
    lambda = (scale <= 0.0) ? 1.0 : scale;
    k = (shape <= 0.0) ? 1.0 : shape;
    kInv = 1.0 / k;
}

double WeibullRand::f(double x) const
{
    if (x < 0)
        return 0;
    if (x == 0) {
        if (k == 1)
            return 1.0 / lambda;
        return (k > 1) ? 0.0 : INFINITY;
    }
    double xAdj = x / lambda;
    double xAdjPow = std::pow(xAdj, k - 1);
    return k / lambda * xAdjPow * std::exp(-xAdj * xAdjPow);
}

double WeibullRand::F(double x) const
{
    if (x <= 0)
        return 0;
    return 1 - std::exp(-std::pow(x / lambda, k));
}

double WeibullRand::Variate() const
{
    return lambda * std::pow(ExponentialRand::StandardVariate(), kInv);
}

double WeibullRand::Mean() const
{
    return lambda * std::tgamma(1 + kInv);
}

double WeibullRand::Variance() const
{
    double res = std::tgamma(1 + kInv);
    res *= -res;
    res += std::tgamma(1 + kInv + kInv);
    return lambda * lambda * res;
}

double WeibullRand::quantileImpl(double p) const
{
    double x = -std::log1p(-p);
    x = std::pow(x, kInv);
    return lambda * x;
}

double WeibullRand::quantileImpl1m(double p) const
{
    double x = -std::log(p);
    x = std::pow(x, kInv);
    return lambda * x;
}

double WeibullRand::Median() const
{
    return lambda * std::pow(M_LN2, kInv);
}

double WeibullRand::Mode() const
{
    if (k <= 1)
        return 0;
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

double WeibullRand::Entropy() const
{
    return M_EULER * (1.0 - kInv) + std::log(lambda * kInv) + 1.0;
}
