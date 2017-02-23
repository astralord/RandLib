#include "KolmogorovSmirnovRand.h"
#include "ExponentialRand.h"
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
    return (x > 1.0) ? 1.0 - K(x) : L(x);
}

double KolmogorovSmirnovRand::CDFCompl(double x)
{
    return (x > 1.0) ? K(x) : 1.0 - L(x);
}

double KolmogorovSmirnovRand::f(double x) const
{
    return PDF(x);
}

double KolmogorovSmirnovRand::logf(double x) const
{
    return std::log(f(x));
}

double KolmogorovSmirnovRand::F(double x) const
{
    return CDF(x);
}

double KolmogorovSmirnovRand::S(double x) const
{
    return CDFCompl(x);
}

double KolmogorovSmirnovRand::truncatedGammaVariate() const
{
    /// Generator for truncated gamma distribution with shape = 1.5
    static constexpr long double tp = 2.193245422464302l; /// π^2 / (8 * 0.75^2)
    int iter = 0;
    do {
        double E0 = 1.2952909208355123l * ExponentialRand::StandardVariate();
        double E1 = 2 * ExponentialRand::StandardVariate();
        double G = tp + E0;
        if (E0 * E0 <= tp * E1 * (G + tp))
            return G;
        double Wp = E0 / tp;
        if (Wp - std::log1p(Wp) <= E1)
            return G;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN;
}

double KolmogorovSmirnovRand::variateForTheLeftMostInterval() const
{
    int iter1 = 0;
    do {
        double G = truncatedGammaVariate();
        double X = M_PI / std::sqrt(8 * G);
        double W = 0.0;
        double Z = 0.5 / G;
        int n = 1, iter2 = 0;
        double Q = 1.0;
        double U = UniformRand::StandardVariate();
        while (U >= W && ++iter2 <= MAX_ITER_REJECTION) {
            W += Z * Q;
            if (U >= W)
                return X;
            n += 2;
            int nSq = n * n;
            Q = std::exp(G - G * nSq);
            W -= nSq * Q;
        }
    } while (++iter1 <= MAX_ITER_REJECTION);
    return NAN;
}

double KolmogorovSmirnovRand::variateForTheRightMostInterval() const
{
    static constexpr double tSq = 0.5625; /// square of parameter t suggested in the book
    int iter1 = 0;
    do {
        double E = ExponentialRand::StandardVariate();
        double U = UniformRand::StandardVariate();
        double X = std::sqrt(tSq + 0.5 * E);
        double W = 0.0;
        int n = 1, iter2 = 0;
        double Z = -2 * X * X;
        while (U > W && ++iter2 < MAX_ITER_REJECTION) {
            ++n;
            int nSq = n * n;
            W += nSq * std::exp(Z * (nSq - 1));
            if (U >= W)
                return X;
            ++n;
            nSq = n * n;
            W -= nSq * std::exp(Z * (nSq - 1));
        }
    } while (++iter1 <= MAX_ITER_REJECTION);
    return NAN;
}

double KolmogorovSmirnovRand::Variate() const
{
    /// Luc Devroye, pp. 163-165
    /// alternating series method
    bool isLeft = UniformRand::StandardVariate() < 0.3728329582237386; /// F(0.75)
    return isLeft ? variateForTheLeftMostInterval() : variateForTheRightMostInterval();
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


