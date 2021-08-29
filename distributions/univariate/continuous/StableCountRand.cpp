#include "StableCountRand.h"


template < typename RealType >
StableCountRand<RealType>::StableCountRand(double exponent, double scale, double location)
{
    SetExponent(exponent);
    SetScale(scale);
    SetLocation(location);
    X.SetSkewness(1);
}

template < typename RealType >
String StableCountRand<RealType>::Name() const
{
    return "Stable count("
            + this->toStringWithPrecision(GetExponent()) + ", "
            + this->toStringWithPrecision(GetScale()) + ", "
            + this->toStringWithPrecision(GetLocation()) + ")";
}

template < typename RealType >
void StableCountRand<RealType>::SetExponent(double exponent)
{
    if (exponent < 0.1 || exponent >= 1.0)
        throw std::invalid_argument("Stable count distribution: exponent should be in the interval [0.1, 1), but it's equal to "
                                    + std::to_string(exponent));

    /// the following error should be removed in the future
    if (exponent != 1.0 && std::fabs(exponent - 1.0) < 0.01)
        throw std::invalid_argument("Stable count distribution: exponent close to 1 with non-zero skewness is not yet supported");

    alpha = exponent;
    logAlpha = std::log(alpha);
    gammaAlphaInv = std::tgammal(1 / alpha);
    lgammaAlphaInv = std::lgammal(1 / alpha);
    lgamma2AlphaInv = std::lgammal(2 / alpha);

    X.SetExponent(exponent);
    double scale = std::pow(std::cos(0.5 * alpha * M_PI), 1 / alpha);
    X.SetScale(scale);
}

template < typename RealType >
void StableCountRand<RealType>::SetLocation(double location)
{
    nu = location;
}

template < typename RealType >
void StableCountRand<RealType>::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Stable count distribution: scale should be positive, but it's equal to "
                                    + std::to_string(scale));
    theta = scale;
    logTheta = std::log(scale);
}

template < typename RealType >
long double StableCountRand<RealType>::Mean() const
{
    return nu + theta * std::exp(lgamma2AlphaInv - lgammaAlphaInv);
}

template < typename RealType >
long double StableCountRand<RealType>::Variance() const
{
    double lgamma3AlphaInv = std::lgamma(3 / alpha);
    double frac1 = 0.5 * std::exp(lgamma3AlphaInv - lgammaAlphaInv);
    double frac2 = std::exp(lgamma2AlphaInv - lgammaAlphaInv);
    double var = frac1 - frac2 * frac2;
    return theta * theta * var;
}

template < typename RealType >
double StableCountRand<RealType>::logf_adj(const RealType &x) const
{
    if (x <= 0)
        return -INFINITY;
    if (alpha == 0.5) {
        // Levy case -> Gamma distribution
        double y = -0.25 * x;
        y += 0.5 * (std::log(x) - M_LNPI);
        return y - 2 * M_LN2;
    }
    double y = logAlpha - lgammaAlphaInv;
    y -= std::log(x);
    y += X.logf(1 / x);
    return y;
}

template < typename RealType >
double StableCountRand<RealType>::logf(const RealType &x) const
{
    double x_adj = (x - nu) / theta;
    double y = logf_adj(x_adj);
    return y - logTheta;
}

template < typename RealType >
double StableCountRand<RealType>::f_adj(const RealType &x) const
{
    if (x <= 0)
        return 0.0;
    if (alpha == 0.5) // Gamma distribution
        return std::exp(logf_adj(x));
    return alpha * X.f(1. / x) / (gammaAlphaInv * x);
}

template < typename RealType >
double StableCountRand<RealType>::f(const RealType &x) const
{
    double x_adj = (x - nu) / theta;
    return f_adj(x_adj) / theta;
}

template < typename RealType >
double StableCountRand<RealType>::F_adj(const RealType &x) const
{
    if (x <= 0)
        return 0.0;
    if (alpha == 0.5) {
        // Levy case -> Gamma distribution
        static constexpr double logA = M_LN3 - M_LN2;
        static constexpr double lgammaA = 0.5 * M_LNPI - M_LN2;
        double x_adj = 0.25 * x;
        return RandMath::pgamma(1.5, x_adj, logA, lgammaA);
    }
    if (x <= 1) { // direct integration
        return RandMath::integral([this] (double t)
        {
            return f_adj(t);
        },
        0, x);
    }
    return 1 - S_adj(x);
}

template < typename RealType >
double StableCountRand<RealType>::F(const RealType &x) const
{
    double x_adj = (x - nu) / theta;
    return F_adj(x_adj);
}

template < typename RealType >
double StableCountRand<RealType>::S_adj(const RealType &x) const
{
    if (x <= 0)
        return 1.0;
    if (alpha == 0.5) {
        // Levy case -> Gamma distribution
        static constexpr double logA = M_LN3 - M_LN2;
        static constexpr double lgammaA = 0.5 * M_LNPI - M_LN2;
        double x_adj = 0.25 * x;
        return RandMath::qgamma(1.5, x_adj, logA, lgammaA);
    }
    if (x >= 1) { // direct integration
        return RandMath::integral([this] (double t)
        {
            if (t == 0)
                return 0.0;
            return alpha * X.f(t) / (gammaAlphaInv * t);
        },
        0, 1.0 / x);
    }
    return 1 - F_adj(x);
}


template < typename RealType >
double StableCountRand<RealType>::S(const RealType &x) const
{
    double x_adj = (x - nu) / theta;
    return S_adj(x_adj);
}

template < typename RealType >
RealType StableCountRand<RealType>::logsumexp(std::vector<RealType> &sample) const
{
    auto max_elem = *std::max_element(sample.begin(), sample.end());
    auto sum = std::accumulate(sample.begin(), sample.end(), 0,
       [max_elem](RealType a, RealType b) { return a + exp(b - max_elem); });
    return max_elem + log(sum);
}

template < typename RealType >
RealType StableCountRand<RealType>::Variate() const
{
    // TODO: firure out how to sample small numbers properly
    // rejection sampling?
    static RealType T_inf = 100000;
    static constexpr RealType cutoff = 100000;
    RealType x = 0, T_n = 0, T_np1 = 0;
    while (T_np1 < T_inf) {
        RealType tau = X.Variate();
        while (tau > cutoff) {
            // workaround for proper simulation
            // of small random values
            tau = X.Variate();
        }
        T_n = T_np1;
        T_np1 += tau;
        ++x;
    }

    return std::exp(std::log(x) / alpha - std::log(T_n));
}

template class StableCountRand<float>;
template class StableCountRand<double>;
template class StableCountRand<long double>;
