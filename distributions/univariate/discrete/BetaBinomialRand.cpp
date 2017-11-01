#include "BetaBinomialRand.h"
#include "BinomialRand.h"
#include <thread>

BetaBinomialRand::BetaBinomialRand(int number, double shape1, double shape2)
{
    SetParameters(number, shape1, shape2);
}

String BetaBinomialRand::Name() const
{
    return "Beta-Binomial(" + toStringWithPrecision(GetNumber()) + ", "
                            + toStringWithPrecision(GetAlpha()) + ", "
                            + toStringWithPrecision(GetBeta()) + ")";
}

void BetaBinomialRand::SetParameters(int number, double shape1, double shape2)
{
    if (shape1 <= 0.0 || shape2 <= 0.0)
        throw std::invalid_argument("Beta-Binomial distribution: shape parameters should be positive");
    if (number <= 0)
        throw std::invalid_argument("Beta-Binomial distribution: number should be positive");
    n = number;
    B.SetShapes(shape1, shape2);
    pmfCoef = RandMath::lfact(n);
    pmfCoef -= std::lgamma(B.GetAlpha() + B.GetBeta() + n);
    pmfCoef -= B.GetLogBetaFunction();
}

double BetaBinomialRand::P(const int & k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

double BetaBinomialRand::logP(const int & k) const
{
    if (k < 0 || k > n)
        return 0.0;
    double y = std::lgamma(k + B.GetAlpha());
    y += std::lgamma(n - k + B.GetBeta());
    y -= RandMath::lfact(k);
    y -= RandMath::lfact(n - k);
    return pmfCoef + y;
}

double BetaBinomialRand::F(const int & k) const
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

int BetaBinomialRand::Variate() const
{
    double p = B.Variate();
    return BinomialDistribution::Variate(n, p);
}

double BetaBinomialRand::Mean() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    return n * alpha / (alpha + beta);
}

double BetaBinomialRand::Variance() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    double alphaPBeta = alpha + beta;
    double numerator = n * alpha * beta * (alphaPBeta + n);
    double denominator = alphaPBeta * alphaPBeta;
    denominator *= (alphaPBeta + 1);
    return numerator / denominator;
}

int BetaBinomialRand::Mode() const
{
    /// for small n we use direct comparison of probabilities
    if (n < 30) {
        double maxValue = logP(0);
        int index = 0;
        for (int i = 1; i <= n; ++i)
        {
            double value = logP(i);
            if (maxValue < value)
            {
                maxValue = value;
                index = i;
            }
        }
        return index;
    }
    /// otherwise use numerical procedure to solve the equation f'(x) = 0
    double guess = n * B.Mean();
    double alpha = B.GetAlpha(), beta = B.GetBeta();
    if (RandMath::findRoot([this, alpha, beta] (double x)
    {
        double y = RandMath::digamma(x + alpha);
        y -= RandMath::digamma(n - x + beta);
        y -= RandMath::digamma(x + 1);
        y += RandMath::digamma(n - x + 1);
        return y;
    }, 0, n, guess))
        return std::round(guess);
    /// if we can't find quantile, then probably something bad has happened
    return -1;
}

double BetaBinomialRand::Skewness() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    double alphaPBeta = alpha + beta;
    double res = (1 + alphaPBeta) / (n * alpha * beta * (alphaPBeta + n));
    res = std::sqrt(res);
    res *= (alphaPBeta + 2 * n) * (beta - alpha);
    res /= alphaPBeta + 2;
    return res;
}

double BetaBinomialRand::ExcessKurtosis() const
{
    double alpha = B.GetAlpha();
    double beta = B.GetBeta();
    double alphaPBeta = alpha + beta;
    double alphaBetaN = alpha * beta * n;
    double res = alpha * beta * (n - 2);
    res += 2 * (double)n * n;
    res -= alphaBetaN * (6 - n) / alphaPBeta;
    res -= 6 * alphaBetaN * n / (alphaPBeta * alphaPBeta);
    res *= 3;
    res += alphaPBeta * (alphaPBeta - 1 + 6 * n);
    res *= alphaPBeta * alphaPBeta * (1 + alphaPBeta);
    res /= (alphaBetaN * (alphaPBeta + 2) * (alphaPBeta + 3) * (alphaPBeta + n));
    return res - 3.0;
}

