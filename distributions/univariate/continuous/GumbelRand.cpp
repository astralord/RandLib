#include "GumbelRand.h"
#include "ExponentialRand.h"

GumbelRand::GumbelRand(double location, double scale)
{
    SetLocation(location);
    SetScale(scale);
}

std::string GumbelRand::Name() const
{
    return "Gumbel(" + toStringWithPrecision(GetLocation()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void GumbelRand::SetLocation(double location)
{
    mu = location;
}

void GumbelRand::SetScale(double scale)
{
    beta = scale > 0.0 ? scale : 1.0;
    logBeta = std::log(beta);
}

double GumbelRand::f(const double & x) const
{
    return std::exp(logf(x));
}

double GumbelRand::logf(const double & x) const
{
    double z = (mu - x) / beta;
    double y = std::exp(z);
    return z - y - logBeta;
}

double GumbelRand::F(const double & x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
    return std::exp(-y);
}

double GumbelRand::S(const double & x) const
{
    double y = (mu - x) / beta;
    y = std::exp(y);
    return -std::expm1(-y);
}

double GumbelRand::Variate() const
{
    return GumbelRand::Variate(mu, beta);
}

double GumbelRand::Variate(double location, double scale)
{
    double w = ExponentialRand::StandardVariate();
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

double GumbelRand::quantileImpl(double p) const
{
    return mu - beta * std::log(-std::log(p));
}

double GumbelRand::quantileImpl1m(double p) const
{
    return mu - beta * std::log(-std::log1p(-p));
}

double GumbelRand::Median() const
{
    static constexpr double M_LN_LN2 = std::log(M_LN2);
    return mu - beta * M_LN_LN2;
}

double GumbelRand::Mode() const
{
    return mu;
}

double GumbelRand::Skewness() const
{
    static constexpr double skew = 12 * M_SQRT2 * M_SQRT3 * M_APERY / (M_PI_SQ * M_PI);
    return skew;
}

double GumbelRand::ExcessKurtosis() const
{
    return 2.4;
}

double GumbelRand::Entropy() const
{
    return logBeta + M_EULER + 1.0;
}
