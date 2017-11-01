#include "PoissonRand.h"
#include "../continuous/UniformRand.h"
#include "../continuous/ExponentialRand.h"

PoissonRand::PoissonRand(double rate)
{
    SetRate(rate);
}

String PoissonRand::Name() const
{
    return "Poisson(" + toStringWithPrecision(GetRate()) + ")";
}

void PoissonRand::SetGeneratorConstants()
{
    delta = std::round(std::sqrt(2 * mu * std::log(M_1_PI * 128 * mu)));
    delta = std::max(6.0, std::min(mu, delta));

    c1 = std::sqrt(0.5 * M_PI * mu * M_E);

    c2 = 1.0 / (2 * mu + delta);
    c2 += 0.5 * (M_LNPI - M_LN2 + std::log(mu + 0.5 * delta));
    c2 = c1 + std::exp(c2);

    c3 = c2 + 1.0;
    c4 = c3 + 1.0;

    zeta = (4.0 * mu) / delta + 2;
    c = std::exp(-(2.0 + delta) / zeta);
    c *= zeta;
    c += c4;

    sqrtMu = std::sqrt(mu);
    sqrtMupHalfDelta = std::sqrt(mu + 0.5 * delta);
    lfactMu = RandMath::lfact(mu);
}

void PoissonRand::SetRate(double rate)
{
    if (rate <= 0.0)
        throw std::invalid_argument("Poisson distribution: rate should be positive");
    lambda = rate;

    logLambda = std::log(lambda);
    mu = std::floor(lambda);
    Fmu = F(mu);
    Pmu = P(mu);

    if (!generateByInversion())
        SetGeneratorConstants();
}

double PoissonRand::P(const int & k) const
{
    return (k < 0) ? 0.0 : std::exp(logP(k));
}

double PoissonRand::logP(const int & k) const
{
    if (k < 0)
        return -INFINITY;
    double y = k * logLambda - lambda;
    return y - RandMath::lfact(k);
}

double PoissonRand::F(const int & k) const
{
    return (k >= 0.0) ? RandMath::qgamma(k + 1, lambda, logLambda) : 0.0;
}

double PoissonRand::S(const int & k) const
{
    return (k >= 0.0) ? RandMath::pgamma(k + 1, lambda, logLambda) : 1.0;
}

double PoissonRand::acceptanceFunction(int X) const
{
    if (X == 0)
        return 0.0;
    double q = X * logLambda;
    q += lfactMu;
    q -= RandMath::lfact(X + mu);
    return q;
}

bool PoissonRand::generateByInversion() const
{
    /// the inversion generator is much faster than rejection one,
    /// however the lost of precision for large rate affects it drastically
    return lambda < 10;
}

int PoissonRand::variateRejection() const
{
    int iter = 0;
    int X = 0;
    do {
        bool reject = false;
        double W = 0.0;
        double U = c * UniformRand::StandardVariate();
        if (U <= c1) {
            double N = NormalRand::StandardVariate();
            double Y = -std::fabs(N) * sqrtMu;
            X = std::floor(Y);
            if (X < -mu) {
                reject = true;
            }
            else {
                W = -0.5 * (N * N - 1.0);
            }
        }
        else if (U <= c2) {
            double N = NormalRand::StandardVariate();
            double Y = 1.0 + std::fabs(N) * sqrtMupHalfDelta;
            X = std::ceil(Y);
            if (X > delta) {
                reject = true;
            }
            else {
                W = Y * (2.0 - Y) / (2.0 * mu + delta);
            }
        }
        else if (U <= c3) {
            return mu;
        }
        else if (U <= c4) {
            X = 1;
        }
        else {
            double V = ExponentialRand::StandardVariate();
            double Y = delta + V * zeta;
            X = std::ceil(Y);
            W = -(2.0 + Y) / zeta;
        }

        if (!reject && W - ExponentialRand::StandardVariate() <= acceptanceFunction(X)) {
            return X + mu;
        }

    } while (++iter < MAX_ITER_REJECTION);
    return -1;
}

int PoissonRand::variateInversion() const
{
    double U = UniformRand::StandardVariate();
    int k = mu;
    double s = Fmu, p = Pmu;
    if (s < U)
    {
        do {
            ++k;
            p *= lambda / k;
            s += p;
        } while (s < U && p > 0);
    }
    else
    {
        s -= p;
        while (k > 0 && s > U) {
            p /= lambda / k;
            --k;
            s -= p;
        }
    }
    return k;
}

int PoissonRand::Variate() const
{
    return generateByInversion() ? variateInversion() : variateRejection();
}

int PoissonRand::Variate(double rate)
{
    /// check validness of parameter
    if (rate <= 0.0)
        return -1;
    if (rate > 1000) {
        /// approximate with normal distribution
        double X = NormalRand::StandardVariate();
        return std::floor(rate + std::sqrt(rate) * X);
    }
    int k = -1;
    double s = 0;
    do {
        s += ExponentialRand::StandardVariate();
        ++k;
    } while (s < rate);
    return k;
}

void PoissonRand::Sample(std::vector<int> &outputData) const
{
    if (generateByInversion()) {
        for (int & var : outputData)
            var = variateInversion();
    }
    else {
        for (int & var : outputData)
            var = variateRejection();
    }
}

double PoissonRand::Mean() const
{
    return lambda;
}

double PoissonRand::Variance() const
{
    return lambda;
}

std::complex<double> PoissonRand::CFImpl(double t) const
{
    std::complex<double> y(std::cos(t) - 1.0, std::sin(t));
    return std::exp(lambda * y);
}

int PoissonRand::Median() const
{
    /// this value is approximate
    return std::max(std::floor(lambda + 1.0 / 3 - 0.02 / lambda), 0.0);
}

int PoissonRand::Mode() const
{
    return mu;
}

double PoissonRand::Skewness() const
{
    return 1.0 / std::sqrt(lambda);
}

double PoissonRand::ExcessKurtosis() const
{
    return 1.0 / lambda;
}

void PoissonRand::Fit(const std::vector<int> &sample)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    SetRate(GetSampleMean(sample));
}

void PoissonRand::Fit(const std::vector<int> &sample, DoublePair &confidenceInterval, double significanceLevel)
{
    size_t n = sample.size();

    if (significanceLevel <= 0 || significanceLevel > 1)
        throw std::invalid_argument(fitErrorDescription(WRONG_LEVEL, "Alpha is equal to " + toStringWithPrecision(significanceLevel)));

    Fit(sample);

    double halfAlpha = 0.5 * significanceLevel;
    ErlangRand ErlangRV(n);
    confidenceInterval.first = ErlangRV.Quantile(halfAlpha);
    ErlangRV.SetShape(n + 1);
    confidenceInterval.second = ErlangRV.Quantile1m(halfAlpha);
}

GammaRand PoissonRand::FitBayes(const std::vector<int> &sample, const GammaDistribution &priorDistribution)
{
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, NON_NEGATIVITY_VIOLATION));
    double alpha = priorDistribution.GetShape();
    double beta = priorDistribution.GetRate();
    GammaRand posteriorDistribution(alpha + GetSampleSum(sample), beta + sample.size());
    SetRate(posteriorDistribution.Mean());
    return posteriorDistribution;
}
