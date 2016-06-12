#include "BetaBinomialRand.h"
#include "BinomialRand.h"
#include <thread>

BetaBinomialRand::BetaBinomialRand(int number, double shape1, double shape2)
{
    setParameters(number, shape1, shape2);
}

std::string BetaBinomialRand::name()
{
    return "Beta-Binomial(" + toStringWithPrecision(getNumber()) + ", "
                            + toStringWithPrecision(getAlpha()) + ", "
                            + toStringWithPrecision(getBeta()) + ")";
}

void BetaBinomialRand::setParameters(int number, double shape1, double shape2)
{
    n = std::max(number, 1);
    B.setParameters(shape1, shape2);
}

double BetaBinomialRand::P(int k) const
{
    if (k < 0 || k > n)
        return 0.0;
    return RandMath::binomialCoef(n, k) * RandMath::betaFun(k + B.getAlpha(), n - k + B.getBeta()) * B.getInverseBetaFunction();
}

double BetaBinomialRand::F(int k) const
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

int BetaBinomialRand::variate() const
{
    double p = B.variate();
    return BinomialRand::variate(n, p);
}

double BetaBinomialRand::Mean() const
{
    double alpha = B.getAlpha();
    double beta = B.getBeta();
    return n * alpha / (alpha + beta);
}

double BetaBinomialRand::Variance() const
{
    double alpha = B.getAlpha();
    double beta = B.getBeta();
    double alphaPBeta = alpha + beta;
    double numerator = n * alpha * beta * (alphaPBeta + n);
    double denominator = alphaPBeta * alphaPBeta;
    denominator *= (alphaPBeta + 1);
    return numerator / denominator;
}

int BetaBinomialRand::Mode() const
{
    double maxValue = 0.0;
    int index = 0;
    for (int i = 0; i <= n; ++i)
    {
        double value = P(i);
        if (maxValue < value)
        {
            maxValue = value;
            index = i;
        }
    }
    return index;
}

double BetaBinomialRand::Skewness() const
{
    double alpha = B.getAlpha();
    double beta = B.getBeta();
    double alphaPBeta = alpha + beta;
    double res = (1 + alphaPBeta) / (n * alpha * beta * (alphaPBeta + n));
    res = std::sqrt(res);
    res *= (alphaPBeta + n + n) * (beta - alpha);
    res /= (alphaPBeta + 2);
    return res;
}

double BetaBinomialRand::ExcessKurtosis() const
{
    double alpha = B.getAlpha();
    double beta = B.getBeta();
    double alphaPBeta = alpha + beta;
    double alphaBetaN = alpha * beta * n;
    double res = alpha * beta * (n - 2);
    res += 2 * n * n;
    res -= alphaBetaN * (6 - n) / alphaPBeta;
    res -= 6 * alphaBetaN * n / (alphaPBeta * alphaPBeta);
    res *= 3;
    res += alphaPBeta * (alphaPBeta - 1 + 6 * n);
    res *= alphaPBeta * alphaPBeta * (1 + alphaPBeta);
    res /= (alphaBetaN * (alphaPBeta + 2) * (alphaPBeta + 3) * (alphaPBeta + n));
    return res - 3.0;
}

