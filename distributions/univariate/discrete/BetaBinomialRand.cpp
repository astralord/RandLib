#include "BetaBinomialRand.h"
#include "BinomialRand.h"
#include "UniformDiscreteRand.h"

template< typename IntType >
BetaBinomialRand<IntType>::BetaBinomialRand(IntType number, double shape1, double shape2)
{
    SetParameters(number, shape1, shape2);
}

template< typename IntType >
String BetaBinomialRand<IntType>::Name() const
{
    return "Beta-Binomial(" + this->toStringWithPrecision(GetNumber()) + ", "
                            + this->toStringWithPrecision(GetAlpha()) + ", "
                            + this->toStringWithPrecision(GetBeta()) + ")";
}

template< typename IntType >
void BetaBinomialRand<IntType>::SetParameters(IntType number, double shape1, double shape2)
{
    if (shape1 <= 0.0 || shape2 <= 0.0)
        throw std::invalid_argument("Beta-Binomial distribution: shape parameters should be positive");
    if (number <= 0)
        throw std::invalid_argument("Beta-Binomial distribution: number should be positive");
    n = number;
    B.SetShapes(shape1, shape2);
    pmfCoef = RandMath::lfact(n);
    pmfCoef -= std::lgammal(B.GetAlpha() + B.GetBeta() + n);
    pmfCoef -= B.GetLogBetaFunction();
}

template< typename IntType >
double BetaBinomialRand<IntType>::P(const IntType &k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

template< typename IntType >
double BetaBinomialRand<IntType>::logP(const IntType &k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    double y = std::lgammal(k + B.GetAlpha());
    y += std::lgammal(n - k + B.GetBeta());
    y -= RandMath::lfact(k);
    y -= RandMath::lfact(n - k);
    return pmfCoef + y;
}

template< typename IntType >
double BetaBinomialRand<IntType>::F(const IntType & k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    double sum = 0.0;
    int i = 0;
    do {
        sum += P(i);
    } while (++i <= k);
    return sum;
}

template< typename IntType >
IntType BetaBinomialRand<IntType>::VariateUniform() const
{
    return UniformDiscreteRand<IntType>::StandardVariate(0, n, this->localRandGenerator);
}

template< typename IntType >
IntType BetaBinomialRand<IntType>::VariateBeta() const
{
    double p = B.Variate();
    return BinomialDistribution<IntType>::Variate(n, p, this->localRandGenerator);
}

template< typename IntType >
IntType BetaBinomialRand<IntType>::Variate() const
{
    return (B.GetAlpha() == 1 && B.GetBeta() == 1) ? VariateUniform() : VariateBeta();
}

template< typename IntType >
void BetaBinomialRand<IntType>::Sample(std::vector<IntType> &outputData) const
{
    if (B.GetAlpha() == 1 && B.GetBeta() == 1) {
        for (IntType & var : outputData)
            var = VariateUniform();
    }
    else {
        for (IntType & var : outputData)
            var = VariateBeta();
    }
}

template< typename IntType >
void BetaBinomialRand<IntType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    B.Reseed(seed + 1);
}

template< typename IntType >
long double BetaBinomialRand<IntType>::Mean() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    return n * alpha / (alpha + beta);
}

template< typename IntType >
long double BetaBinomialRand<IntType>::Variance() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    double alphaPBeta = alpha + beta;
    double numerator = n * alpha * beta * (alphaPBeta + n);
    double denominator = alphaPBeta * alphaPBeta;
    denominator *= (alphaPBeta + 1);
    return numerator / denominator;
}

template< typename IntType >
IntType BetaBinomialRand<IntType>::Mode() const
{
    /// for small n we use direct comparison of probabilities
    static constexpr int SMALL_N = 32;
    if (n < SMALL_N) {
        std::vector<double> logProbs(n + 1);
        for (int i = 0; i <= n; ++i)
            logProbs[i] = logP(i);
        std::vector<double>::iterator maxVar = std::max_element(logProbs.begin(), logProbs.end());
        return std::distance(logProbs.begin(), maxVar);
    }
    /// otherwise use numerical procedure to solve the equation f'(x) = 0
    double guess = n * B.Mean();
    double alpha = B.GetAlpha(), beta = B.GetBeta();
    if (!RandMath::findRoot<double>([this, alpha, beta] (double x)
    {
        double y = RandMath::digamma(x + alpha);
        y -= RandMath::digamma(n - x + beta);
        y -= RandMath::digamma(x + 1);
        y += RandMath::digamma(n - x + 1);
        return y;
    }, 0, n, guess))
        throw std::runtime_error("Beta-Binomial distribution: failure in numerical procedure");
    return std::round(guess);
}

template< typename IntType >
long double BetaBinomialRand<IntType>::Skewness() const
{
    long double alpha = B.GetAlpha();
    long double beta = B.GetBeta();
    long double alphaPBeta = alpha + beta;
    long double res = (1 + alphaPBeta) / (n * alpha * beta * (alphaPBeta + n));
    res = std::sqrt(res);
    res *= (alphaPBeta + 2 * n) * (beta - alpha);
    res /= alphaPBeta + 2;
    return res;
}

template< typename IntType >
long double BetaBinomialRand<IntType>::ExcessKurtosis() const
{
    long double alpha = B.GetAlpha();
    long double beta = B.GetBeta();
    long double alphaPBeta = alpha + beta;
    long double alphaBetaN = alpha * beta * n;
    long double res = alpha * beta * (n - 2);
    res += 2 * (double)n * n;
    res -= alphaBetaN * (6 - n) / alphaPBeta;
    res -= 6 * alphaBetaN * n / (alphaPBeta * alphaPBeta);
    res *= 3;
    res += alphaPBeta * (alphaPBeta - 1 + 6 * n);
    res *= alphaPBeta * alphaPBeta * (1 + alphaPBeta);
    res /= (alphaBetaN * (alphaPBeta + 2) * (alphaPBeta + 3) * (alphaPBeta + n));
    return res - 3.0;
}

template class BetaBinomialRand<int>;
template class BetaBinomialRand<long int>;
template class BetaBinomialRand<long long int>;
