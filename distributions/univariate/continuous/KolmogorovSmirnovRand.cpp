#include "KolmogorovSmirnovRand.h"
#include "ExponentialRand.h"
#include "UniformRand.h"

template < typename RealType >
KolmogorovSmirnovRand<RealType>::KolmogorovSmirnovRand()
{
}

template < typename RealType >
String KolmogorovSmirnovRand<RealType>::Name() const
{
    return "Kolmogorov-Smirnov";
}

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::L(RealType x)
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

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::K(RealType x)
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

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::f(const RealType & x) const
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

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::logf(const RealType & x) const
{
    return std::log(f(x));
}

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::F(const RealType & x) const
{
    return (x > 1.0) ? 1.0 - K(x) : L(x);
}

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::S(const RealType & x) const
{
    return (x > 1.0) ? K(x) : 1.0 - L(x);
}

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::logF(const RealType & x) const
{
    return (x > 1.0) ? std::log1pl(-K(x)) : std::log(L(x));
}

template < typename RealType >
double KolmogorovSmirnovRand<RealType>::logS(const RealType & x) const
{
    return (x > 1.0) ? std::log(K(x)) : std::log1pl(-L(x));
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::truncatedGammaVariate() const
{
    /// Generator for truncated gamma distribution with shape = 1.5
    static constexpr long double tp = 2.193245422464302l; ///< π^2 / (8 * 0.75^2)
    static constexpr long double rate = 1.2952909208355123l;
    size_t iter = 0;
    do {
        RealType E0 = rate * ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType E1 = 2 * ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType G = tp + E0;
        if (E0 * E0 <= tp * E1 * (G + tp))
            return G;
        RealType Wp = E0 / tp;
        if (Wp - std::log1pl(Wp) <= E1)
            return G;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Kolmogorov-Smirnov distribution: sampling failed");
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::variateForTheLeftMostInterval() const
{
    size_t iter1 = 0;
    do {
        RealType G = truncatedGammaVariate();
        RealType X = M_PI / std::sqrt(8 * G);
        RealType W = 0.0;
        RealType Z = 0.5 / G;
        size_t n = 1, iter2 = 0;
        RealType Q = 1.0;
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        while (U >= W && ++iter2 <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION) {
            W += Z * Q;
            if (U >= W)
                return X;
            n += 2;
            int nSq = n * n;
            Q = std::exp(G - G * nSq);
            W -= nSq * Q;
        }
    } while (++iter1 <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Kolmogorov-Smirnov distribution: sampling failed");
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::variateForTheRightMostInterval() const
{
    static constexpr double tSq = 0.5625; /// square of parameter t suggested in the book
    size_t iter1 = 0;
    do {
        RealType E = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType X = std::sqrt(tSq + 0.5 * E);
        RealType W = 0.0;
        size_t n = 1, iter2 = 0;
        RealType Z = -2 * X * X;
        while (U > W && ++iter2 <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION) {
            ++n;
            int nSq = n * n;
            W += nSq * std::exp(Z * (nSq - 1));
            if (U >= W)
                return X;
            ++n;
            nSq = n * n;
            W -= nSq * std::exp(Z * (nSq - 1));
        }
    } while (++iter1 <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    throw std::runtime_error("Kolmogorov-Smirnov distribution: sampling failed");
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::Variate() const
{
    /// Luc Devroye, pp. 163-165
    /// alternating series method
    bool isLeft = UniformRand<RealType>::StandardVariate(this->localRandGenerator) < 0.3728329582237386; /// F(0.75)
    return isLeft ? variateForTheLeftMostInterval() : variateForTheRightMostInterval();
}

template < typename RealType >
long double KolmogorovSmirnovRand<RealType>::Mean() const
{
    return M_SQRTPI * M_SQRT1_2 * M_LN2;
}

template < typename RealType >
long double KolmogorovSmirnovRand<RealType>::Variance() const
{
    long double mean = Mean();
    return M_PI_SQ / 12 - mean * mean;
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::Mode() const
{
    return 0.735467812776958l;
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::Median() const
{
    return 0.82757355518990761l;
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::quantileImpl(double p) const
{
    RealType guess = std::sqrt(-0.5 * (std::log1pl(-p) - M_LN2));
    if (p < 1e-5) {
        double logP = std::log(p);
        if (!RandMath::findRootNewtonFirstOrder<RealType>([this, logP] (RealType x)
        {
            double logCdf = logF(x), logPdf = logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            return DoublePair(first, second);
        }, guess))
            throw std::runtime_error("Kolmogorov-Smirnov distribution: failure in numerical procedure");
        return guess;
    }
    if (!RandMath::findRootNewtonFirstOrder<RealType>([p, this] (RealType x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        throw std::runtime_error("Kolmogorov-Smirnov distribution: failure in numerical procedure");
    return guess;
}

template < typename RealType >
RealType KolmogorovSmirnovRand<RealType>::quantileImpl1m(double p) const
{
    RealType guess = std::sqrt(-0.5 * std::log(0.5 * p));
    if (p < 1e-5) {
        double logP = std::log(p);
        if (!RandMath::findRootNewtonFirstOrder<RealType>([this, logP] (RealType x)
        {
            double logCcdf = logS(x), logPdf = logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            return DoublePair(first, second);
        }, guess))
            throw std::runtime_error("Kolmogorov-Smirnov distribution: failure in numerical procedure");
        return guess;
    }
    if (!RandMath::findRootNewtonFirstOrder<RealType>([p, this] (RealType x)
    {
        double first = p - S(x);
        double second = f(x);
        return DoublePair(first, second);
    }, guess))
        throw std::runtime_error("Kolmogorov-Smirnov distribution: failure in numerical procedure");
    return guess;
}



template class KolmogorovSmirnovRand<float>;
template class KolmogorovSmirnovRand<double>;
template class KolmogorovSmirnovRand<long double>;
