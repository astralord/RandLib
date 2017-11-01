#include "BinomialRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/NormalRand.h"
#include "../continuous/ExponentialRand.h"
#include "BernoulliRand.h"

BinomialDistribution::BinomialDistribution(int number, double probability)
{
    SetParameters(number, probability);
}

void BinomialDistribution::SetGeneratorConstants()
{
    minpq = std::min(p, q);
    npFloor = std::floor(n * minpq);
    pFloor = npFloor / n;
    pRes = RandMath::areClose(npFloor, n * minpq) ? 0.0 : minpq - pFloor;

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

void BinomialDistribution::SetParameters(int number, double probability)
{
    if (probability < 0.0 || probability > 1.0)
        throw std::invalid_argument("Binomial distribution: probability parameter should in interval [0, 1]");
    if (number <= 0)
        throw std::invalid_argument("Binomial distribution: number should be positive");
    n = number;
    p = probability;
    q = 1.0 - p;
    np = n * p;
    lfactn = RandMath::lfact(n);
    logProb = std::log(p);
    log1mProb = std::log1p(-p);
    SetGeneratorConstants();
}

double BinomialDistribution::logProbFloor(int k) const
{
    double y = lfactn;
    y -= RandMath::lfact(n - k);
    y -= RandMath::lfact(k);
    y += k * logPFloor;
    y += (n - k) * logQFloor;
    return y;
}

double BinomialDistribution::P(const int & k) const
{
    return (k < 0 || k > n) ? 0.0 : std::exp(logP(k));
}

double BinomialDistribution::logP(const int & k) const
{
    if (k < 0 || k > n)
        return -INFINITY;
    double y = lfactn;
    y -= RandMath::lfact(n - k);
    y -= RandMath::lfact(k);
    y += k * logProb;
    y += (n - k) * log1mProb;
    return y;
}

double BinomialDistribution::F(const int & k) const
{
    if (k < 0)
        return 0.0;
    if (k >= n)
        return 1.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = RandMath::lfact(n - kp1);
    logBetaFun += RandMath::lfact(k);
    logBetaFun -= lfactn;
    return RandMath::ibeta(q, nmk, kp1, logBetaFun, log1mProb, logProb);
}

double BinomialDistribution::S(const int & k) const
{
    if (k < 0)
        return 1.0;
    if (k >= n)
        return 0.0;
    int nmk = n - k, kp1 = k + 1;
    double logBetaFun = RandMath::logBeta(kp1, nmk);
    return RandMath::ibeta(p, kp1, nmk, logBetaFun, logProb, log1mProb);
}

BinomialDistribution::GENERATOR_ID BinomialDistribution::GetIdOfUsedGenerator() const
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

int BinomialDistribution::variateRejection() const
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

        if (!reject) {
            X += npFloor;
            if (X >= 0 && X <= n && V <= logProbFloor(X) - logPnpInv)
                return X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return -1;
}

int BinomialDistribution::variateWaiting(int number) const
{
    /// waiting algorithm, using
    /// sum of geometrically distributed variables
    int X = -1, sum = 0;
    do {
        sum += G.Variate() + 1;
        ++X;
    } while (sum <= number);
    return X;
}

int BinomialDistribution::variateWaiting(int number, double probability)
{
    int X = -1, sum = 0;
    do {
        sum += GeometricRand::Variate(probability) + 1;
        ++X;
    } while (sum <= number);
    return X;
}

int BinomialDistribution::Variate() const
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
        /// if X ~ Bin(n, p') and Y ~ Bin(n - X, (p - p') / (1 - p'))
        /// then Z = X + Y ~ Bin(n, p)
        int Z = variateRejection();
        if (pRes > 0)
            Z += variateWaiting(n - Z);
        return (p > 0.5) ? n - Z : Z;
    }
    case BERNOULLI_SUM:
    default:
        return variateBernoulliSum(n, p);
    }
    return -1; /// unexpected return
}

int BinomialDistribution::variateBernoulliSum(int number, double probability)
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

int BinomialDistribution::Variate(int number, double probability)
{
    /// sanity check
    if (number < 0 || probability < 0.0 || probability > 1.0)
        return -1;
    if (probability == 0.0)
        return 0;
    if (probability == 1.0)
        return number;

    if (number < 10)
        return variateBernoulliSum(number, probability);
    if (probability < 0.5)
        return variateWaiting(number, probability);
    return number - variateWaiting(number, 1.0 - probability);
}

void BinomialDistribution::Sample(std::vector<int> &outputData) const
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

double BinomialDistribution::Mean() const
{
    return np;
}

double BinomialDistribution::Variance() const
{
    return np * q;
}

std::complex<double> BinomialDistribution::CFImpl(double t) const
{
    std::complex<double> y(q + p * std::cos(t), p * std::sin(t));
    return std::pow(y, n);
}

int BinomialDistribution::Median() const
{
    return std::round(np);
}

int BinomialDistribution::Mode() const
{
    return std::floor(np + p);
}

double BinomialDistribution::Skewness() const
{
    return (q - p) / std::sqrt(np * q);
}

double BinomialDistribution::ExcessKurtosis() const
{
    double y = 1.0 / (p * q);
    y -= 6.0;
    return y / n;
}

void BinomialDistribution::FitProbability(const std::vector<int> &sample)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    if (!allElementsAreNotBiggerThan(n, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(n)));
    SetParameters(n, GetSampleMean(sample) / n);
}

BetaRand BinomialDistribution::FitProbabilityBayes(const std::vector<int> &sample, const BetaDistribution &priorDistribution)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    if (!allElementsAreNotBiggerThan(n, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(n)));
    int N = sample.size();
    double sum = GetSampleSum(sample);
    double alpha = priorDistribution.GetAlpha();
    double beta = priorDistribution.GetBeta();
    BetaRand posteriorDistribution(sum + alpha, N * n - sum + beta);
    SetParameters(n, posteriorDistribution.Mean());
    return posteriorDistribution;
}

BetaRand BinomialDistribution::FitProbabilityMinimax(const std::vector<int> &sample)
{
    double shape = 0.5 * std::sqrt(n);
    BetaRand B(shape, shape);
    return FitProbabilityBayes(sample, B);
}

String BinomialRand::Name() const
{
    return "Binomial(" + toStringWithPrecision(GetNumber()) + ", " + toStringWithPrecision(GetProbability()) + ")";
}
