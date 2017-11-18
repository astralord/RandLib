#include "ParetoRand.h"
#include "UniformRand.h"

ParetoRand::ParetoRand(double shape, double scale)
{
    SetShape(shape);
    SetScale(scale);
}

String ParetoRand::Name() const
{
    return "Pareto(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetScale()) + ")";
}

void ParetoRand::SetShape(double shape)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Pareto distribution: shape should be positive");
    alpha = shape;
    logAlpha = std::log(alpha);
}

void ParetoRand::SetScale(double scale)
{
    if (scale <= 0.0)
        throw std::invalid_argument("Pareto distribution: scale should be positive");
    sigma = scale;
    logSigma = std::log(sigma);
}

double ParetoRand::f(const double & x) const
{
    return (x < sigma) ? 0.0 : std::exp(logf(x));
}

double ParetoRand::logf(const double & x) const
{
    if (x < sigma)
        return -INFINITY;
    double logX = std::log(x);
    double y = logSigma - logX;
    y *= alpha;
    y -= logX;
    y += logAlpha;
    return y;
}

double ParetoRand::F(const double & x) const
{
    return (x > sigma) ? -std::expm1l(alpha * std::log(sigma / x)) : 0.0;
}

double ParetoRand::S(const double & x) const
{
    return (x > sigma) ? std::pow(sigma / x, alpha) : 1.0;
}

double ParetoRand::variateForAlphaOne(RandGenerator &randGenerator)
{
    return 1.0 / UniformRand::StandardVariate(randGenerator);
}

double ParetoRand::variateForAlphaTwo(RandGenerator &randGenerator)
{
    return 1.0 / std::sqrt(UniformRand::StandardVariate(randGenerator));
}

double ParetoRand::variateForGeneralAlpha(double shape, RandGenerator &randGenerator)
{
    return std::exp(ExponentialRand::StandardVariate(randGenerator) / shape);
}

double ParetoRand::Variate() const
{
    return sigma * StandardVariate(alpha, localRandGenerator);
}

double ParetoRand::StandardVariate(double shape, RandGenerator &randGenerator)
{
    if (RandMath::areClose(shape, 1.0))
        return variateForAlphaOne(randGenerator);
    if (RandMath::areClose(shape, 2.0))
        return variateForAlphaTwo(randGenerator);
    return variateForGeneralAlpha(shape, randGenerator);
}

void ParetoRand::Sample(std::vector<double> &outputData) const
{
    if (RandMath::areClose(alpha, 1.0)) {
        for (double &var : outputData)
            var = sigma * variateForAlphaOne(localRandGenerator);
    }
    else if (RandMath::areClose(alpha, 2.0)) {
        for (double &var : outputData)
            var = sigma * variateForAlphaTwo(localRandGenerator);
    }
    else {
        for (double &var : outputData)
            var = sigma * variateForGeneralAlpha(alpha, localRandGenerator);
    }
}

long double ParetoRand::Mean() const
{
    return (alpha > 1) ? alpha * sigma / (alpha - 1) : INFINITY;
}

long double ParetoRand::Variance() const
{
    if (alpha > 2)
    {
        double var = sigma / (alpha - 1);
        var *= var;
        return alpha * var / (alpha - 2);
    }
    return (alpha > 1) ? INFINITY : NAN;
}

double ParetoRand::Median() const
{
    return std::exp(logSigma + M_LN2 / alpha);
}

double ParetoRand::Mode() const
{
    return sigma;
}

long double ParetoRand::Skewness() const
{
    if (alpha <= 3)
        return INFINITY;
    double skewness = (alpha - 2.0) / alpha;
    skewness = std::sqrt(skewness);
    skewness *= (1 + alpha) / (alpha - 3);
    return skewness + skewness;
}

long double ParetoRand::ExcessKurtosis() const
{
    if (alpha <= 4)
        return INFINITY;
    double numerator = alpha + 1;
    numerator *= alpha;
    numerator -= 6;
    numerator *= alpha;
    numerator -= 2;
    double denominator = alpha * (alpha - 3) * (alpha - 4);
    return 6.0 * numerator / denominator;
}

double ParetoRand::quantileImpl(double p) const
{
    double y = logSigma - std::log1pl(-p) / alpha;
    return std::exp(y);
}

double ParetoRand::quantileImpl1m(double p) const
{
    double y = logSigma - std::log(p) / alpha;
    return std::exp(y);
}

double ParetoRand::Entropy() const
{
    return logSigma - logAlpha + 1.0 / alpha + 1;
}

void ParetoRand::FitShape(const std::vector<double> &sample, bool unbiased)
{
    if (!allElementsAreNotSmallerThan(sigma, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(sigma)));
    double invShape = GetSampleLogMean(sample) - logSigma;
    if (invShape == 0.0)
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, "Possibly all the elements of the sample coincide with the lower boundary σ."));
    double shape = 1.0 / invShape;
    if (unbiased)
        shape *= (1.0 - 1.0 / sample.size());
    SetShape(shape);
}

void ParetoRand::FitScale(const std::vector<double> &sample, bool unbiased)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double scale = unbiased ? (1.0 - 1.0 / (sample.size() * alpha)) * minVar : minVar;
    SetScale(scale);
}

void ParetoRand::Fit(const std::vector<double> &sample, bool unbiased)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    if (minVar <= 0)
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));

    double logBiasedSigma = std::log(minVar);
    double invShape = GetSampleLogMean(sample) - logBiasedSigma;
    if (invShape == 0.0)
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, "Possibly all the elements of the sample are the same."));
    double shape = 1.0 / invShape;

    double scale = minVar;
    if (unbiased) {
        int n = sample.size();
        shape *= 1.0 - 2.0 / n;
        scale *= 1.0 - invShape / (n - 1);
    }
    SetScale(scale);
    SetShape(shape);
}

GammaRand ParetoRand::FitShapeBayes(const std::vector<double> &sample, const GammaDistribution &priorDistribution, bool MAP)
{
    if (!allElementsAreNotSmallerThan(sigma, sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(sigma)));
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newRate = priorDistribution.GetRate() + n * (GetSampleLogMean(sample) - logSigma);
    GammaRand posteriorDistribution(newShape, newRate);
    SetShape(MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}
