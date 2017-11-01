#include "GammaRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

GammaDistribution::GammaDistribution(double shape, double rate)
{
    SetParameters(shape, rate);
}

void GammaDistribution::SetParameters(double shape, double rate)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Gamma distribution: shape should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Gamma distribution: rate should be positive");

    alpha = shape > 0 ? shape : 1.0;
    
    beta = (rate > 0.0) ? rate : 1.0;
    theta = 1.0 / beta;

    lgammaAlpha = std::lgamma(alpha);
    logAlpha = std::log(alpha);
    logBeta = std::log(beta);
    pdfCoef = -lgammaAlpha + alpha * logBeta;

    if (getIdOfUsedGenerator(alpha) == SMALL_SHAPE) {
        /// set constants for generator
        genCoef.t = 0.5 * std::log1p(-alpha);
        genCoef.t = 0.07 + 0.75 * std::exp(genCoef.t);
        genCoef.b = 1.0 + std::exp(-genCoef.t) * alpha / genCoef.t;
    }
}

void GammaDistribution::SetShape(double shape)
{
    SetParameters(shape, beta);
}

double GammaDistribution::f(const double & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0)
    {
        if (alpha > 1.0)
            return 0.0;
        return (alpha == 1.0) ? beta : INFINITY;
    }
    return std::exp(logf(x));
}

double GammaDistribution::logf(const double & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0)
    {
        if (alpha > 1.0)
            return -INFINITY;
        return (alpha == 1.0) ? logBeta : INFINITY;
    }
    double y = (alpha - 1.0) * std::log(x);
    y -= x * beta;
    y += pdfCoef;
    return y;
}

double GammaDistribution::F(const double & x) const
{
    return (x > 0.0) ? RandMath::pgamma(alpha, x * beta, logAlpha, lgammaAlpha) : 0.0;
}

double GammaDistribution::logF(const double &x) const
{
    return (x > 0.0) ? RandMath::lpgamma(alpha, x * beta, logAlpha, lgammaAlpha) : -INFINITY;
}

double GammaDistribution::S(const double & x) const
{
    return (x > 0.0) ? RandMath::qgamma(alpha, x * beta, logAlpha, lgammaAlpha) : 1.0;
}

double GammaDistribution::logS(const double &x) const
{
    return (x > 0.0) ? RandMath::lqgamma(alpha, x * beta, logAlpha, lgammaAlpha) : 0.0;
}

GammaDistribution::GENERATOR_ID GammaDistribution::getIdOfUsedGenerator(double shape)
{
    if (shape < 0.34)
        return SMALL_SHAPE;
    if (shape <= 3.0 && RandMath::areClose(shape, std::round(shape)))
        return INTEGER_SHAPE;
    if (RandMath::areClose(shape, 1.5))
        return ONE_AND_A_HALF_SHAPE;
    if (shape > 1.0 && shape < 1.2)
        return FISHMAN;
    return MARSAGLIA_TSANG;
}

double GammaDistribution::variateThroughExponentialSum(int shape)
{
    double X = 0.0;
    for (int i = 0; i < shape; ++i)
        X += ExponentialRand::StandardVariate();
    return X;
}

double GammaDistribution::variateForShapeOneAndAHalf()
{
    double W = ExponentialRand::StandardVariate();
    double N = NormalRand::StandardVariate();
    return W + 0.5 * N * N;
}

double GammaDistribution::variateBest() const
{
    /// Algorithm RGS for gamma variates (Best, 1983)
    double X = 0;
    int iter = 0;
    do {
        double V = genCoef.b * UniformRand::StandardVariate();
        double W = UniformRand::StandardVariate();
        if (V <= 1) {
            X = genCoef.t * std::pow(V, 1.0 / alpha);
            if (W <= (2.0 - X) / (2.0 + X) || W <= std::exp(-X))
                return X;
        }
        else {
            X = -std::log(genCoef.t * (genCoef.b - V) / alpha);
            double Y = X / genCoef.t;
            if (W * (alpha + Y - alpha * Y) <= 1 || W <= std::pow(Y, alpha - 1))
                return X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}


double GammaDistribution::variateAhrensDieter(double shape)
{
    /// Rejection algorithm GS for gamma variates (Ahrens and Dieter, 1974)
    double X = 0;
    int iter = 0;
    double shapeInv = 1.0 / shape;
    double t = shapeInv + M_1_E;
    do {
        double U = UniformRand::StandardVariate();
        double p = shape * t * U;
        double W = ExponentialRand::StandardVariate();
        if (p <= 1)
        {
            X = std::pow(p, shapeInv);
            if (X <= W)
                return X;
        }
        else
        {
            X = -std::log(t * (1 - U));
            if ((1 - shape) * std::log(X) <= W)
                return X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

double GammaDistribution::variateFishman(double shape)
{
    /// G. Fishman algorithm (shape > 1)
    double W1, W2;
    double shapem1 = shape - 1;
    do {
        W1 = ExponentialRand::StandardVariate();
        W2 = ExponentialRand::StandardVariate();
    } while (W2 < shapem1 * (W1 - std::log(W1) - 1));
    return shape * W1;
}

double GammaDistribution::variateMarsagliaTsang(double shape)
{
    /// Marsaglia and Tsang’s Method (shape > 1/3)
    double d = shape - 1.0 / 3;
    double c = 3 * std::sqrt(d);
    int iter = 0;
    do {
        double N;
        do {
            N = NormalRand::StandardVariate();
        } while (N <= -c);
        double v = 1 + N / c;
        v = v * v * v;
        N *= N;
        double U = UniformRand::StandardVariate();
        if (U < 1.0 - 0.331 * N * N || std::log(U) < 0.5 * N + d * (1.0 - v + std::log(v))) {
            return d * v;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

double GammaDistribution::StandardVariate(double shape)
{
    if (shape <= 0)
        return NAN;

    GENERATOR_ID genId = getIdOfUsedGenerator(shape);

    switch(genId) {
    case INTEGER_SHAPE:
        return variateThroughExponentialSum(std::round(shape));
    case ONE_AND_A_HALF_SHAPE:
        return variateForShapeOneAndAHalf();
    case SMALL_SHAPE:
        return variateAhrensDieter(shape);
    case FISHMAN:
        return variateFishman(shape);
    case MARSAGLIA_TSANG:
        return variateMarsagliaTsang(shape);
    default:
        return NAN;
    }
}

double GammaDistribution::Variate(double shape, double rate)
{
    return (shape <= 0.0 || rate <= 0.0) ? NAN : StandardVariate(shape) / rate;
}

double GammaDistribution::Variate() const
{
    GENERATOR_ID genId = getIdOfUsedGenerator(alpha);

    switch(genId) {
    case INTEGER_SHAPE:
        return theta * variateThroughExponentialSum(alpha);
    case ONE_AND_A_HALF_SHAPE:
        return theta * variateForShapeOneAndAHalf();
    case SMALL_SHAPE:
        return theta * variateBest();
    case FISHMAN:
        return theta * variateFishman(alpha);
    case MARSAGLIA_TSANG:
        return theta * variateMarsagliaTsang(alpha);
    default:
        return NAN;
    }
}

void GammaDistribution::Sample(std::vector<double> &outputData) const
{
    GENERATOR_ID genId = getIdOfUsedGenerator(alpha);

    switch(genId) {
    case INTEGER_SHAPE:
        for (double &var : outputData)
            var = theta * variateThroughExponentialSum(alpha);
        break;
    case ONE_AND_A_HALF_SHAPE:
        for (double &var : outputData)
            var = theta * variateForShapeOneAndAHalf();
        break;
    case SMALL_SHAPE:
        for (double &var : outputData)
            var = theta * variateBest();
        break;
    case FISHMAN:
        for (double &var : outputData)
            var = theta * variateFishman(alpha);
        break;
    case MARSAGLIA_TSANG:
        for (double &var : outputData)
            var = theta * variateMarsagliaTsang(alpha);
        break;
    default:
        return;
    }
}

double GammaDistribution::Mean() const
{
    return alpha * theta;
}

double GammaDistribution::GeometricMean() const
{
    return RandMath::digamma(alpha) - logBeta;
}

double GammaDistribution::Variance() const
{
    return alpha * theta * theta;
}

double GammaDistribution::GeometricVariance() const
{
    return RandMath::trigamma(alpha);
}

double GammaDistribution::Mode() const
{
    return (alpha <= 1) ? 0 : (alpha - 1) * theta;
}

double GammaDistribution::Skewness() const
{
    return 2.0 / std::sqrt(alpha);
}

double GammaDistribution::ExcessKurtosis() const
{
    return 6.0 / alpha;
}

double GammaDistribution::initRootForSmallP(double r) const
{
    double root = 0;
    double c[5];
    c[4] = 1;
    /// first coefficient
    double denominator = alpha + 1;
    c[3] = 1.0 / denominator;
    /// second coefficient
    denominator *= denominator;
    denominator *= alpha + 2;
    c[2] = 0.5 * (3 * alpha + 5) / denominator;
    /// third coefficient
    denominator *= 3 * (alpha + 1) * (alpha + 3);
    c[1] = 8 * alpha + 33;
    c[1] *= alpha;
    c[1] += 31;
    c[1] /= denominator;
    /// fourth coefficient
    denominator *= 8 * (alpha + 1) * (alpha + 2) * (alpha + 4);
    c[0] = 125 * alpha + 1179;
    c[0] *= alpha;
    c[0] += 3971;
    c[0] *= alpha;
    c[0] += 5661;
    c[0] *= alpha;
    c[0] += 2888;
    c[0] /= denominator;
    /// now calculate root
    for (int i = 0; i != 5; ++i) {
        root += c[i];
        root *= r;
    }
    return root;
}

double GammaDistribution::initRootForLargeP(double logQ) const
{
    /// look for approximate value of x -> INFINITY
    double x = (logQ + lgammaAlpha) / alpha;
    x = -std::exp(x) / alpha;
    return -alpha * RandMath::Wm1Lambert(x);
}

double GammaDistribution::initRootForLargeShape(double p) const
{
    if (p == 0.5)
        return alpha;
    double x = RandMath::erfcinv(2 * p);
    double lambda = x * x / alpha + 1;
    lambda = -std::exp(-lambda);
    lambda = (x < 0) ? -RandMath::Wm1Lambert(lambda) : -RandMath::W0Lambert(lambda);
    return lambda * alpha;
}

double GammaDistribution::initRootForLargeShape1m(double p) const
{
    if (p == 0.5)
        return alpha;
    double x = -RandMath::erfcinv(2 * p);
    double lambda = x * x / alpha + 1;
    lambda = -std::exp(-lambda);
    lambda = (x < 0) ? -RandMath::Wm1Lambert(lambda) : -RandMath::W0Lambert(lambda);
    return lambda * alpha;
}

double GammaDistribution::quantileInitialGuess(double p) const
{
    /// Method is taken from
    /// "Efficient and accurate algorithms
    /// for the computation and inversion
    /// of the incomplete gamma function ratios"
    /// (Amparo Gil, Javier Segura and Nico M. Temme)
    double guess = 0;
    if (alpha < 10) {
        double r = std::log(p * alpha) + lgammaAlpha;
        r = std::exp(r / alpha);
        /// if p -> 0
        if (r < 0.2 * (alpha + 1)) {
            guess = initRootForSmallP(r);
        }
        else {
            double logQ = std::log1p(-p);
            /// boundary adviced in a paper
            double maxBoundary1 = -0.5 * alpha - logAlpha - lgammaAlpha;
            /// the maximum possible value to have a solution
            double maxBoundary2 = alpha * (logAlpha - 1) - lgammaAlpha;
            /// if p -> 1
            if (logQ < std::min(maxBoundary1, maxBoundary2))
                guess = initRootForLargeP(logQ);
            else if (alpha < 1)
                guess = r;
            else
                guess = initRootForLargeShape(p);
        }
    }
    else
        guess = initRootForLargeShape(p);
    return guess / beta;
}

double GammaDistribution::quantileInitialGuess1m(double p) const
{
    if (alpha < 10) {
        double logQ = std::log(p);
        /// boundary adviced in a paper
        double maxBoundary1 = -0.5 * alpha - logAlpha - lgammaAlpha;
        /// the maximum possible value to have a solution
        double maxBoundary2 = alpha * (logAlpha - 1) - lgammaAlpha;
        /// if p -> 0
        if (logQ < std::min(maxBoundary1, maxBoundary2))
            return initRootForLargeP(logQ) / beta;
    }
    else {
        return initRootForLargeShape1m(p) / beta;
    }
    return quantileInitialGuess(1.0 - p);
}

double GammaDistribution::df(double x) const
{
    double z = (alpha - 1) - beta * x;
    double y = (alpha - 2) * std::log(x);
    y -= beta * x;
    y += pdfCoef;
    return z * std::exp(y);
}

double GammaDistribution::dfDivf(double x) const
{
    return x / (alpha - 1 - beta * x);
}

double GammaDistribution::quantileImpl(double p) const
{
    double guess = quantileInitialGuess(p);
    if (p < 1e-5) { /// too small p
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
            if (x <= 0)
               return DoubleTriplet(-INFINITY, 0, 0);
            double logCdf = logF(x), logPdf = logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            double third = second * (dfDivf(x) - second);
            return DoubleTriplet(first, second, third);
        }, guess))
            return guess;
        /// if we can't find quantile, then probably something bad has happened
        return NAN;
    }
    if (RandMath::findRoot([this, p] (double x)
    {
        if (x <= 0)
            return DoubleTriplet(-p, 0, 0);
        double first = F(x) - p;
        double second = f(x);
        double third = df(x);
        return DoubleTriplet(first, second, third);
    }, guess))
        return guess;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

double GammaDistribution::quantileImpl1m(double p) const
{
    double guess = quantileInitialGuess1m(p);
    if (p < 1e-5) { /// too small p
        double logP = std::log(p);
        if (RandMath::findRoot([this, logP] (double x)
        {
           if (x <= 0)
               return DoubleTriplet(logP, 0, 0);
            double logCcdf = logS(x), logPdf = logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            double third = second * (dfDivf(x) + second);
            return DoubleTriplet(first, second, third);
        }, guess))
            return guess;
        /// if we can't find quantile, then probably something bad has happened
        return NAN;
    }
    if (RandMath::findRoot([this, p] (double x)
    {
        if (x <= 0)
            return DoubleTriplet(p - 1.0, 0, 0);
        double first = p - S(x);
        double second = f(x);
        double third = df(x);
        return DoubleTriplet(first, second, third);
    }, guess))
        return guess;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

std::complex<double> GammaDistribution::CFImpl(double t) const
{
    return std::pow(std::complex<double>(1.0, -theta * t), -alpha);
}

void FreeScaleGammaDistribution::SetRate(double rate)
{
    SetParameters(alpha, rate);
}

void FreeScaleGammaDistribution::SetScale(double scale)
{
    SetRate(1.0 / scale);
}

void FreeScaleGammaDistribution::FitRate(const std::vector<double> &sample, bool unbiased)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double mean = GetSampleMean(sample);
    double coef = alpha - (unbiased ? 1.0 / sample.size() : 0.0);
    SetParameters(alpha, coef / mean);
}

GammaRand FreeScaleGammaDistribution::FitRateBayes(const std::vector<double> &sample, const GammaDistribution & priorDistribution)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double alpha0 = priorDistribution.GetShape();
    double beta0 = priorDistribution.GetRate();
    double newAlpha = alpha * sample.size() + alpha0;
    double newBeta = GetSampleSum(sample) + beta0;
    GammaRand posteriorDistribution(newAlpha, newBeta);
    SetParameters(alpha, posteriorDistribution.Mean());
    return posteriorDistribution;
}

String GammaRand::Name() const
{
    return "Gamma(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetRate()) + ")";
}

void GammaRand::FitShape(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));

    /// Calculate log-average
    long double logAverage = 0.0L;
    for (double var : sample)
        logAverage += std::log(var);
    logAverage /= sample.size();

    /// Calculate initial guess via method of moments
    double shape = GetSampleMean(sample) * beta;
    /// Run root-finding procedure
    double s = logAverage + logBeta;
    if (!RandMath::findRoot([s] (double x)
    {
        double first = RandMath::digamma(x) - s;
        double second = RandMath::trigamma(x);
        return DoublePair(first, second);
    }, shape))
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetParameters(shape, beta);
}

void GammaRand::Fit(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsArePositive(sample))
        throw std::invalid_argument(fitErrorDescription(WRONG_SAMPLE, POSITIVITY_VIOLATION));

    /// Calculate average and log-average
    double average = GetSampleMean(sample);
    long double logAverage = 0.0L;
    for (double var : sample)
        logAverage += std::log(var);
    logAverage /= sample.size();

    /// Calculate initial guess for shape
    double s = std::log(average) - logAverage;
    double sm3 = s - 3.0, sp12 = 12.0 * s;
    double shape = sm3 * sm3 + 2 * sp12;
    shape = std::sqrt(shape);
    shape -= sm3;
    shape /= sp12;

    if (!RandMath::findRoot([s] (double x)
    {
        double first = RandMath::digammamLog(x) + s;
        double second = RandMath::trigamma(x) - 1.0 / x;
        return DoublePair(first, second);
    }, shape))
        throw std::runtime_error(fitErrorDescription(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetParameters(shape, shape / average);
}

String ChiSquaredRand::Name() const
{
    return "Chi-squared(" + toStringWithPrecision(GetDegree()) + ")";
}

void ChiSquaredRand::SetDegree(size_t degree)
{
    GammaDistribution::SetParameters(0.5 * degree, 0.5);
}


String ErlangRand::Name() const
{
    return "Erlang(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetRate()) + ")";
}

void ErlangRand::SetParameters(size_t shape, double rate)
{
    GammaDistribution::SetParameters(shape, rate);
}

void ErlangRand::SetShape(size_t shape)
{
    SetParameters(shape, beta);
}
