#include "BetaPrimeRand.h"

BetaPrimeRand::BetaPrimeRand(double shape1, double shape2)
{
    SetParameters(shape1, shape2);
}

String BetaPrimeRand::Name() const
{
    return "Beta Prime(" + toStringWithPrecision(GetAlpha()) + ", " + toStringWithPrecision(GetBeta()) + ")";
}

void BetaPrimeRand::SetParameters(double shape1, double shape2)
{
    if (shape1 <= 0 || shape2 <= 0)
        throw std::invalid_argument("Beta-prime distribution: shapes should be positive");
    B.SetShapes(shape1, shape2);
    alpha = B.GetAlpha();
    beta = B.GetBeta();
}

double BetaPrimeRand::f(const double & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0) {
        if (alpha == 1.0)
            return 1.0 / GetBetaFunction();
        return (alpha > 1) ? 0.0 : INFINITY;
    }
    return std::exp(logf(x));
}

double BetaPrimeRand::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0) {
        if (alpha == 1.0)
            return -GetLogBetaFunction();
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    double y = (alpha - 1) * std::log(x);
    y -= (alpha + beta) * std::log1p(x);
    return y - GetLogBetaFunction();
}

double BetaPrimeRand::F(const double & x) const
{
    return (x > 0) ? B.F(x / (1.0 + x)) : 0;
}

double BetaPrimeRand::S(const double & x) const
{
    return (x > 0) ? B.S(x / (1.0 + x)) : 1;
}

double BetaPrimeRand::Variate() const
{
    double x = B.Variate();
    return x / (1.0 - x);
}

void BetaPrimeRand::Sample(std::vector<double> &outputData) const
{
    B.Sample(outputData);
    for (double &var : outputData)
        var = var / (1.0 - var);
}

double BetaPrimeRand::Mean() const
{
    return (beta > 1) ? alpha / (beta - 1) : INFINITY;
}

double BetaPrimeRand::Variance() const
{
    if (beta <= 2)
        return INFINITY;
    double betam1 = beta - 1;
    double numerator = alpha * (alpha + betam1);
    double denominator = (betam1 - 1) * betam1 * betam1;
    return numerator / denominator;
}

double BetaPrimeRand::Median() const
{
    return (alpha == beta) ? 1.0 : quantileImpl(0.5);
}

double BetaPrimeRand::Mode() const
{
    return (alpha < 1) ? 0 : (alpha - 1) / (beta + 1);
}

double BetaPrimeRand::Skewness() const
{
    if (beta <= 3)
        return INFINITY;
    double aux = alpha + beta - 1;
    double skewness = (beta - 2) / (alpha * aux);
    skewness = std::sqrt(skewness);
    aux += alpha;
    aux += aux;
    return aux * skewness / (beta - 3);
}

double BetaPrimeRand::ExcessKurtosis() const
{
    if (beta <= 4)
        return INFINITY;
    double betam1 = beta - 1;
    double numerator = betam1 * betam1 * (beta - 2) / (alpha * (alpha + betam1));
    numerator += 5 * beta - 11;
    double denominator = (beta - 3) * (beta - 4);
    return 6 * numerator / denominator;
}

double BetaPrimeRand::quantileImpl(double p) const
{
    double x = B.Quantile(p);
    return x / (1.0 - x);
}

double BetaPrimeRand::quantileImpl1m(double p) const
{
    double x = B.Quantile1m(p);
    return x / (1.0 - x);
}

std::complex<double> BetaPrimeRand::CFImpl(double t) const
{
    /// if no singularity - simple numeric integration
    if (alpha >= 1)
        return UnivariateDistribution::CFImpl(t);

    double re = ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x) - 1.0;
    }, 0.0, 1.0);
    re += F(1.0);
    re += ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    }, 1.0, INFINITY);

    double im = ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, 0.0, INFINITY);
    return std::complex<double>(re, im);
}
