#include "KolmogorovSmirnovRand.h"
#include "ExponentialRand.h"
#include "UniformRand.h"

KolmogorovSmirnovRand::KolmogorovSmirnovRand()
{
}

String KolmogorovSmirnovRand::Name() const
{
    return "Kolmogorov-Smirnov";
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

double KolmogorovSmirnovRand::f(const double & x) const
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

double KolmogorovSmirnovRand::logf(const double & x) const
{
    return std::log(f(x));
}

double KolmogorovSmirnovRand::F(const double & x) const
{
    return (x > 1.0) ? 1.0 - K(x) : L(x);
}

double KolmogorovSmirnovRand::S(const double & x) const
{
    return (x > 1.0) ? K(x) : 1.0 - L(x);
}

double KolmogorovSmirnovRand::logF(const double & x) const
{
    return (x > 1.0) ? std::log1p(-K(x)) : std::log(L(x));
}

double KolmogorovSmirnovRand::logS(const double & x) const
{
    return (x > 1.0) ? std::log(K(x)) : std::log1p(-L(x));
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

double KolmogorovSmirnovRand::quantileImpl(double p) const
{
    double guess = std::sqrt(-0.5 * (std::log1p(-p) - M_LN2));
    if (p < 1e-5) {
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
            double logCdf = logF(x), logPdf = logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            return DoublePair(first, second);
        }, guess))
            return guess;
        return NAN;
    }
    if (RandMath::findRoot([p, this] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}

double KolmogorovSmirnovRand::quantileImpl1m(double p) const
{
    double guess = std::sqrt(-0.5 * std::log(0.5 * p));
    if (p < 1e-5) {
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
            double logCcdf = logS(x), logPdf = logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            return DoublePair(first, second);
        }, guess))
            return guess;
        return NAN;
    }
    if (RandMath::findRoot([p, this] (double x)
    {
        double first = p - S(x);
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        return guess;
    return NAN;
}


