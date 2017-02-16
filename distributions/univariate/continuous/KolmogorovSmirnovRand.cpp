#include "KolmogorovSmirnovRand.h"
#include "UniformRand.h"

KolmogorovSmirnovRand::KolmogorovSmirnovRand()
{
}

std::string KolmogorovSmirnovRand::Name() const
{
    return "Kolmogorov-Smirnov";
}

double KolmogorovSmirnovRand::PDF(double x)
{
    if (x <= 0.0)
        return 0.0;
    double sum = 0.0, addon = 0.0;
    int k = 1;
    double xSq = x * x;
    if (x < 1.0) {
        double aux = 0.125 / xSq;
        do {
            double temp = M_PI * (2 * k - 1);
            temp *= temp;
            addon = temp - 4 * xSq;
            addon *= std::exp(-temp * aux);
            sum += addon;
            ++k;
        } while (addon > MIN_POSITIVE * sum);
        return M_SQRT2PI * sum * 0.25 / std::pow(x, 4);
    }
    /// x > 1.0
    do {
        int temp = k * k;
        addon = std::exp(-2 * temp * xSq);
        addon *= temp;
        sum += (k & 1) ? addon : -addon;
        ++k;
    } while (addon > MIN_POSITIVE * sum);
    return 8 * sum * x;
}

double KolmogorovSmirnovRand::L(double x)
{
    if (x <= 0.0)
        return 0.0;
    double sum = 0.0, addon = 0.0;
    int k = 1;
    double aux = M_PI_SQ * 0.125 / (x * x);
    do {
        int temp = (2 * k - 1);
        temp *= temp;
        addon = std::exp(-temp * aux);
        sum += addon;
        ++k;
    } while (addon > MIN_POSITIVE * sum);
    return M_SQRT2PI * sum / x;
}

double KolmogorovSmirnovRand::K(double x)
{
    if (x <= 0.0)
        return 1.0;
    double sum = 0.0, addon = 0.0;
    int k = 1;
    double xSq = x * x;
    do {
        int temp = 2 * k * k;
        addon = std::exp(-temp * xSq);
        sum += (k & 1) ? addon : -addon;
        ++k;
    } while (addon > MIN_POSITIVE * sum);
    return 2 * sum;
}

double KolmogorovSmirnovRand::CDF(double x)
{
    return (x < 1.0) ? 1.0 - K(x) : L(x);
}

double KolmogorovSmirnovRand::CDFCompl(double x)
{
    return (x < 1.0) ? K(x) : 1.0 - L(x);
}

double KolmogorovSmirnovRand::f(double x) const
{
    return PDF(x);
}

double KolmogorovSmirnovRand::F(double x) const
{
    return CDF(x);
}

double KolmogorovSmirnovRand::S(double x) const
{
    return CDFCompl(x);
}

double KolmogorovSmirnovRand::Variate() const
{
    return KolmogorovSmirnovRand::Quantile(UniformRand::StandardVariate());
}

double KolmogorovSmirnovRand::Mean() const
{
    return M_SQRTPI * M_SQRT1_2 * M_LN2;
}

double KolmogorovSmirnovRand::Variance() const
{
    double mean = Mean();
    return M_PI_SQ / 12 - mean * mean;
}

double KolmogorovSmirnovRand::Mode() const
{
    return 0.735467812776958;
}

double KolmogorovSmirnovRand::Median() const
{
    return 0.82757355518990761;
}

double KolmogorovSmirnovRand::Quantile(double p)
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return 0.0;
    if (p == 1)
        return INFINITY;
    double guess = std::sqrt(-0.5 * std::log(0.5 - 0.5 * p));
    if (RandMath::findRoot([p] (double x)
    {
        double first = CDF(x) - p;
        double second = PDF(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double KolmogorovSmirnovRand::Quantile1m(double p)
{
    if (p < 0 || p > 1)
        return NAN;
    if (p == 0)
        return INFINITY;
    if (p == 1)
        return 0.0;
    double guess = std::sqrt(-0.5 * std::log(0.5 * p));
    if (RandMath::findRoot([p] (double x)
    {
        double first = p - CDFCompl(x);
        double second = PDF(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double KolmogorovSmirnovRand::quantileImpl(double p) const
{
    return KolmogorovSmirnovRand::Quantile(p);
}

double KolmogorovSmirnovRand::quantileImpl1m(double p) const
{
    return KolmogorovSmirnovRand::Quantile1m(p);
}


