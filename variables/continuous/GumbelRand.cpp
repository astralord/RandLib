#include "GumbelRand.h"

GumbelRand::GumbelRand(double location, double scale)
{
    setLocation(location);
    setScale(scale);
}

std::string GumbelRand::name()
{
    return "Gumbel(" + toStringWithPrecision(getLocation()) + ", " + toStringWithPrecision(getScale()) + ")";
}

void GumbelRand::setLocation(double location)
{
    mu = location;
}

void GumbelRand::setScale(double scale)
{
    beta = scale;
    if (beta <= 0)
        beta = 1.0;
    betaInv = 1.0 / beta;
}

double GumbelRand::f(double x) const
{
    double z = betaInv * (x - mu);
    double y = std::exp(-z);
    y = std::exp(-z - y);
    return betaInv * y;
}

double GumbelRand::F(double x) const
{
    double y = betaInv * (x - mu);
    y = std::exp(-y);
    return std::exp(-y);
}

double GumbelRand::variate() const
{
    return GumbelRand::variate(mu, beta);
}

double GumbelRand::variate(double location, double scale)
{
    double w = ExponentialRand::standardVariate();
    return location - scale * std::log(w);
}

double GumbelRand::Mean() const
{
    return mu + beta * M_EULER;
}

double GumbelRand::Variance() const
{
    double v = M_PI * beta;
    return v * v / 6;
}

double GumbelRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return -INFINITY;
    if (p == 1)
        return INFINITY;
    return mu - beta * std::log(-std::log(p));
}

double GumbelRand::Median() const
{
    static constexpr double LN_M_LN2 = std::log(M_LN2);
    return mu - beta * LN_M_LN2;
}

double GumbelRand::Mode() const
{
    return mu;
}

double GumbelRand::Skewness() const
{
    static constexpr double skew = 12 * M_SQRT2 * M_SQRT3 * M_APERY / (M_PI * M_PI * M_PI);
    return skew;
}

double GumbelRand::ExcessKurtosis() const
{
    return 2.4;
}

double GumbelRand::Entropy() const
{
    return std::log(beta) + M_EULER + 1.0;
}
