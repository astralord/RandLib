#include "BinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/NormalRand.h"
#include "../continuous/ExponentialRand.h"
#include "BernoulliRand.h"

BinomialRand::BinomialRand(int number, double probability)
{
    SetParameters(number, probability);
}

std::string BinomialRand::Name() const
{
    return "Binomial(" + toStringWithPrecision(GetNumber()) + ", " + toStringWithPrecision(GetProbability()) + ")";
}

void BinomialRand::SetGeneratorConstants()
{
    minpq = std::min(p, q);
    npFloor = std::floor(n * minpq);
    pFloor = npFloor / n;
    pRes = (RandMath::areClose(npFloor, n * minpq) ? 0.0 : minpq - pFloor);

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    if (genId == BERNOULLI_SUM)
        return;
    else if (genId == WAITING) {
        G.SetProbability(minpq);
        return;
    }

    nqFloor = n - npFloor;
    double qFloor = 1.0 - pFloor;
    if (pRes > 0)
        G.SetProbability(pRes / qFloor);

    /// Set deltas
    double npq = npFloor * qFloor;
    double coef = 128.0 * n / M_PI;
    delta1 = coef * pFloor / (81.0 * qFloor);
    delta1 = npq * std::log(delta1);
    if (delta1 > 1.0)
        delta1 = std::sqrt(delta1);
    else
        delta1 = 1.0;
    delta2 = coef * qFloor / pFloor;
    delta2 = npq * std::log(delta2);
    if (delta2 > 1.0)
        delta2 = std::sqrt(delta2);
    else
        delta2 = 1.0;

    /// Set sigmas and c
    double npqSqrt = std::sqrt(npq);
    sigma1 = npqSqrt * (1.0 + 0.25 * delta1 / npFloor);
    sigma2 = npqSqrt * (1.0 + 0.25 * delta2 / nqFloor);
    c = 2.0 * delta1 / npFloor;

    /// Set a's
    a1 = 0.5 * std::exp(c) * sigma1 * M_SQRT2PI;

    a2 = 0.5 * sigma2 * M_SQRT2PI;
    a2 += a1;

    coefa3 = 0.5 * delta1 / (sigma1 * sigma1);
    a3 = 1.0 / nqFloor - coefa3;
    a3 *= delta1;
    a3 = std::exp(a3);
    a3 /= coefa3;
    a3 += a2;

    coefa4 = 0.5 * delta2 / (sigma2 * sigma2);
    a4 = std::exp(-delta2 * coefa4) / coefa4;
    a4 += a3;

    logPFloor = std::log(pFloor);
    logQFloor = (pFloor == qFloor) ? logPFloor : std::log1p(-pFloor);

    logPnpInv = logProbFloor(npFloor);
}

void BinomialRand::SetParameters(int number, double probability)
{
    n = std::max(number, 1);
    p = std::min(probability, 1.0);
    p = std::max(p, 0.0);
    q = 1.0 - p;
    np = n * p;
    lgammaNp1 = std::lgamma(n + 1);
    logProb = std::log(p);
    log1mProb = std::log1p(-p);
    SetGeneratorConstants();
}

double BinomialRand::logProbFloor(int k) const
{
    double y = lgammaNp1;
    y -= std::lgamma(n - k);
    y -= std::lgamma(k);
    y += k * logPFloor;
    y += (n - k) * logQFloor;
    return y;
}

double BinomialRand::P(int k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

double BinomialRand::logP(int k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    double y = lgammaNp1;
    y -= std::lgamma(n - k + 1);
    y -= std::lgamma(k + 1);
    y += k * logProb;
    y += (n - k) * log1mProb;
    return y;
}

double BinomialRand::F(int k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    return RandMath::regularizedBetaFun(q, n - k, 1 + k);
}

double BinomialRand::S(int k) const
{
    if (k < 0)
        return 1.0;
    if (k >= n)
        return 0.0;
    return RandMath::regularizedBetaFun(p, 1 + k, n - k);
}

BinomialRand::GENERATOR_ID BinomialRand::GetIdOfUsedGenerator() const
{
    /// if (n is tiny and minpq is big) or p = 0.5 and n is not so large,
    /// we just sum Bernoulli random variables
    if ((n <= 3) || (n <= 13 && minpq > 0.025 * (n + 6)) || (n <= 200 && RandMath::areClose(p, 0.5)))
        return BERNOULLI_SUM;

    /// for small [np] we use simple waiting algorithm
    if ((npFloor <= 12) ||
        (pRes > 0 && npFloor <= 16))
        return WAITING;

    /// otherwise
    return REJECTION;
}

int BinomialRand::variateRejection() const
{
    /// a rejection algorithm by Devroye and Naderlsamanl (1980)
    /// p.533. Non-Uniform Random Variate Generation. Luc Devroye
    /// it can be used only when n * p is integer and p < 0.5
    bool reject = true;
    int iter = 0;
    double X, Y, V;
    do {
        double U = UniformRand::Variate(0, a4);
        if (U <= a1)
        {
            double N = NormalRand::StandardVariate();
            Y = sigma1 * std::fabs(N);
            reject = (Y >= delta1);
            if (!reject)
            {
                double W = ExponentialRand::StandardVariate();
                X = std::floor(Y);
                V = -W - 0.5 * N * N + c;
            }
        }
        else if (U <= a2)
        {
            double N = NormalRand::StandardVariate();
            Y = sigma2 * std::fabs(N);
            reject = (Y >= delta2);
            if (!reject)
            {
                double W = ExponentialRand::StandardVariate();
                X = std::floor(-Y);
                V = -W - 0.5 * N * N;
            }
        }
        else if (U <= a3)
        {
            double W1 = ExponentialRand::StandardVariate();
            double W2 = ExponentialRand::StandardVariate();
            Y = delta1 + W1 / coefa3;
            X = std::floor(Y);
            V = -W2 - coefa3 * Y + delta1 / nqFloor;
            reject = false;
        }
        else
        {
            double W1 = ExponentialRand::StandardVariate();
            double W2 = ExponentialRand::StandardVariate();
            Y = delta2 + W1 / coefa4;
            X = std::floor(-Y);
            V = -W2 - coefa4 * Y;
            reject = false;
        }

        X += npFloor;
        if (!reject && X >= 0 && X <= n && V <= logProbFloor(X) - logPnpInv)
            return X;
    } while (++iter <= MAX_ITER_REJECTION);
    return -1;
}

int BinomialRand::variateWaiting(int number) const
{
    /// waiting algorithm, using
    /// sum of geometrically distributed variables
    int X = -1, sum = 0;
    do {
        sum += G.Variate() + 1.0;
        ++X;
    } while (sum <= number);
    return X;
}

int BinomialRand::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING:
    {
        double var = variateWaiting(n);
        return (p <= 0.5) ? var : n - var;
    }
    case REJECTION:
    {
        /// if X ~ Bin(n, p') and Y ~ Bin(n - Y, (p - p') / (1 - p'))
        /// then Z = X + Y ~ Bin(n, p)
        int Y = variateRejection();
        if (pRes > 0)
            Y += variateWaiting(n - Y);
        return (p > 0.5) ? n - Y : Y;
    }
    case BERNOULLI_SUM:
    default:
        return variateBernoulliSum(n, p);
    }
    return -1; /// unexpected return
}

int BinomialRand::variateBernoulliSum(int number, double probability)
{
    int var = 0;
    if (RandMath::areClose(probability, 0.5)) {
        for (int i = 0; i != number; ++i)
            var += BernoulliRand::StandardVariate();
    }
    else {
        for (int i = 0; i != number; ++i)
            var += BernoulliRand::Variate(probability);
    }
    return var;
}

int BinomialRand::Variate(int number, double probability)
{
    return variateBernoulliSum(number, probability);
}

void BinomialRand::Sample(std::vector<int> &outputData) const
{
    if (p == 0.0) {
        std::fill(outputData.begin(), outputData.end(), 0);
        return;
    }
    if (RandMath::areClose(p, 1.0)) {
        std::fill(outputData.begin(), outputData.end(), n);
        return;
    }

    GENERATOR_ID genId = GetIdOfUsedGenerator();
    switch (genId) {
    case WAITING:
    {
        if (p <= 0.5) {
            for (int &var : outputData)
               var = variateWaiting(n);
        }
        else {
            for (int &var : outputData)
               var = n - variateWaiting(n);
        }
        return;
    }
    case REJECTION:
    {
        for (int &var : outputData)
            var = variateRejection();
        if (pRes > 0) {
            for (int &var : outputData)
                var += variateWaiting(n - var);
        }
        if (p > 0.5) {
            for (int &var : outputData)
               var = n - var;
        }
        return;
    }
    case BERNOULLI_SUM:
    default:
    {
        for (int &var : outputData)
           var = variateBernoulliSum(n, p);
        return;
    }
    }
}

double BinomialRand::Mean() const
{
    return np;
}

double BinomialRand::Variance() const
{
    return np * q;
}

std::complex<double> BinomialRand::CFImpl(double t) const
{
    std::complex<double> y(q + p * std::cos(t), p * std::sin(t));
    return std::pow(y, n);
}

double BinomialRand::Median() const
{
    return std::round(np);
}

int BinomialRand::Mode() const
{
    return std::round(np + p);
}

double BinomialRand::Skewness() const
{
    return (q - p) / std::sqrt(np * q);
}

double BinomialRand::ExcessKurtosis() const
{
    double y = 1.0 / (p * q);
    y -= 6.0;
    return y / n;
}

bool BinomialRand::FitProbabilityMLE(const std::vector<int> &sample)
{
    if (!checkValidity(sample))
        return false;
    SetParameters(n, sampleMean(sample) / n);
    return true;
}

bool BinomialRand::FitProbabilityMM(const std::vector<int> &sample)
{
    return FitProbabilityMLE(sample);
}

bool BinomialRand::FitProbabilityBayes(const std::vector<int> &sample, BetaRand &priorDistribution)
{
    int N = sample.size();
    double sum = sampleSum(sample);
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    priorDistribution.SetParameters(sum + alpha, N * n - sum + beta);
    SetParameters(n, priorDistribution.Mean());
    return true;
}

bool BinomialRand::FitProbabilityMinimax(const std::vector<int> &sample)
{
    double shape = 0.5 * std::sqrt(n);
    BetaRand B(shape, shape);
    return FitProbabilityBayes(sample, B);
}
