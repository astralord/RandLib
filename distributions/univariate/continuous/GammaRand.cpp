#include "GammaRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

GammaRand::GammaRand(double shape, double rate)
{
    SetParameters(shape, rate);
}

std::string GammaRand::Name() const
{
    return "Gamma(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetRate()) + ")";
}

void GammaRand::SetConstantsForGenerator()
{
    t = 0.07 + 0.75 * std::sqrt(1.0 - alpha);
    b = 1.0 + std::exp(-t) * alpha / t;
}

void GammaRand::SetParameters(double shape, double rate)
{
    alpha = shape > 0 ? shape : 1.0;
    alphaInv = 1.0 / alpha;
    
    beta = (rate > 0.0) ? rate : 1.0;
    theta = 1.0 / beta;

    mLgammaShape = -std::lgamma(alpha);
    logBeta = std::log(beta);
    pdfCoef = mLgammaShape + alpha * logBeta;

    if (GetIdOfUsedGenerator(alpha) == SMALL_SHAPE)
        SetConstantsForGenerator();
}

double GammaRand::f(const double & x) const
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

double GammaRand::logf(const double & x) const
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

double GammaRand::F(const double & x) const
{
    return (x > 0.0) ? RandMath::pgamma(alpha, x * beta) : 0.0;
}

double GammaRand::S(const double & x) const
{
    return (x > 0.0) ? RandMath::qgamma(alpha, x * beta) : 1.0;
}

GammaRand::GENERATOR_ID GammaRand::GetIdOfUsedGenerator(double shape)
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

double GammaRand::variateThroughExponentialSum(int shape)
{
    double X = 0.0;
    for (int i = 0; i < shape; ++i)
        X += ExponentialRand::StandardVariate();
    return X;
}

double GammaRand::variateForShapeOneAndAHalf()
{
    double W = ExponentialRand::StandardVariate();
    double N = NormalRand::StandardVariate();
    return W + 0.5 * N * N;
}

double GammaRand::variateBest() const
{
    /// Algorithm RGS for gamma variates (Best, 1983)
    double X = 0;
    int iter = 0;
    do {
        double V = b * UniformRand::StandardVariate();
        double W = UniformRand::StandardVariate();
        if (V <= 1) {
            X = t * std::pow(V, alphaInv);
            if (W <= (2.0 - X) / (2.0 + X) || W <= std::exp(-X))
                return X;
        }
        else {
            X = -std::log(alphaInv * t * (b - V));
            double Y = X / t;
            if (W * (alpha + Y - alpha * Y) <= 1 || W <= std::pow(Y, alpha - 1))
                return X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}


double GammaRand::variateAhrensDieter(double shape)
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

double GammaRand::variateFishman(double shape)
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

double GammaRand::variateMarsagliaTsang(double shape)
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

double GammaRand::StandardVariate(double shape)
{
    if (shape <= 0)
        return NAN;

    GENERATOR_ID genId = GetIdOfUsedGenerator(shape);

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
    return NAN;
}

double GammaRand::Variate(double shape, double rate)
{
    return StandardVariate(shape) / rate;
}

double GammaRand::Variate() const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator(alpha);

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
    return NAN;
}

void GammaRand::Sample(std::vector<double> &outputData) const
{
    GENERATOR_ID genId = GetIdOfUsedGenerator(alpha);

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

double GammaRand::Mean() const
{
    return alpha * theta;
}

double GammaRand::Variance() const
{
    return alpha * theta * theta;
}

double GammaRand::Mode() const
{
    return (alpha <= 1) ? 0 : (alpha - 1) * theta;
}

double GammaRand::Skewness() const
{
    return 2.0 / std::sqrt(alpha);
}

double GammaRand::ExcessKurtosis() const
{
    return 6.0 / alpha;
}

double GammaRand::initRootForSmallP(double r) const
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

double GammaRand::initRootForLargeP(double logQ) const
{
    /// look for approximate value of x -> INFINITY
    double x = (logQ - mLgammaShape) / alpha;
    x = -std::exp(x) / alpha;
    return -alpha * RandMath::Wm1Lambert(x);
}

double GammaRand::initRootForLargeShape(double p) const
{
    if (p == 0.5)
        return alpha;
    double x = RandMath::erfcinv(2 * p);
    double lambda = x * x / alpha + 1;
    lambda = -std::exp(-lambda);
    if (x < 0)
        lambda = -RandMath::Wm1Lambert(lambda);
    else
        lambda = -RandMath::W0Lambert(lambda);
    return lambda * alpha;
}

double GammaRand::df(double x) const
{
    double z = (alpha - 1) - beta * x;
    double y = (alpha - 2) * std::log(x);
    y -= beta * x;
    y += pdfCoef;
    return z * std::exp(y);
}

double GammaRand::quantileInitialGuess(double p) const
{
    /// Method is taken from
    /// "Efficient and accurate algorithms
    /// for the computation and inversion
    /// of the incomplete gamma function ratios"
    /// (Amparo Gil, Javier Segura and Nico M. Temme)
    double guess = 0;
    if (alpha < 10) {
        double r = std::log(p * alpha) - mLgammaShape;
        r = std::exp(r * alphaInv);
        /// if p -> 0
        if (r < 0.2 * (alpha + 1)) {
            guess = initRootForSmallP(r);
        }
        else {
            double logQ = std::log1p(-p);
            double logAlpha = std::log(alpha); // can be hashed
            double maxBoundary1 = -0.5 * alpha - logAlpha + mLgammaShape; /// boundary adviced in a paper
            double maxBoundary2 = alpha * (logAlpha - 1) + mLgammaShape; /// the maximum possible value to have a solution
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

double GammaRand::quantileInitialGuess1m(double p) const
{
    if (alpha < 10) {
        double logQ = std::log(p);
        double logAlpha = std::log(alpha); // can be hashed
        double maxBoundary1 = -0.5 * alpha - logAlpha + mLgammaShape; /// boundary adviced in a paper
        double maxBoundary2 = alpha * (logAlpha - 1) + mLgammaShape; /// the maximum possible value to have a solution
        /// if p -> 0
        if (logQ < std::min(maxBoundary1, maxBoundary2))
            return initRootForLargeP(logQ) / beta;
    }
    return quantileInitialGuess(1.0 - p);
}

double GammaRand::quantileImpl(double p) const
{
    double guess = quantileInitialGuess(p);
    if (p < 1e-5) { /// too small p
        if (RandMath::findRoot([this, p] (double x)
        {
            double cdf = F(x), pdf = f(x);
            double first = std::log(cdf / p);
            double second = pdf / cdf;
            double third = (df(x) * cdf - pdf * pdf) / (cdf * cdf);
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

double GammaRand::quantileImpl1m(double p) const
{
    double guess = quantileInitialGuess1m(p);
    if (p < 1e-5) { /// too small p
        if (RandMath::findRoot([this, p] (double x)
        {
            double ccdf = S(x), pdf = f(x);
            double first = std::log(ccdf / p);
            double second = -pdf / ccdf;
            double third = (-df(x) * ccdf - pdf * pdf) / (ccdf * ccdf);
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

std::complex<double> GammaRand::CFImpl(double t) const
{
    return std::pow(std::complex<double>(1.0, -theta * t), -alpha);
}

void GammaRand::FitScaleMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    SetParameters(alpha, alpha / sampleMean(sample));
}

void GammaRand::FitMLE(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));

    /// Calculate average and log-average
    double average = sampleMean(sample);
    long double logAverage = 0.0L;
    for (double var : sample) {
        logAverage += std::log(var);
    }
    logAverage /= sample.size();

    /// Calculate initial guess for shape
    double s = std::log(average) - logAverage;
    double sm3 = s - 3.0, sp12 = 12.0 * s;
    double shape = sm3 * sm3 + sp12 + sp12;
    shape = std::sqrt(shape);
    shape -= sm3;
    shape /= sp12;

    if (!RandMath::findRoot([s] (double x)
    {
        double first = std::log(x) - RandMath::digamma(x) - s;
        double second = RandMath::trigamma(x) - 1.0 / x;
        return DoublePair(first, second);
    }, shape))
        throw std::runtime_error(fitError(UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetParameters(shape, shape / average);
}

void GammaRand::FitShapeMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    SetParameters(sampleMean(sample) * beta, beta);
}

void GammaRand::FitScaleMM(const std::vector<double> &sample)
{
    FitScaleMLE(sample);
}

void GammaRand::FitMM(const std::vector<double> &sample)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double mu1 = sampleMean(sample);
    double var = sampleVariance(sample, mu1);
    double shape = mu1 * mu1 / var;
    SetParameters(shape, shape / mu1);
}

GammaRand GammaRand::FitRateBayes(const std::vector<double> &sample, const GammaRand & priorDistribution)
{
    /// Sanity check
    if (!allElementsAreNonNegative(sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, POSITIVITY_VIOLATION));
    double alpha0 = priorDistribution.GetShape();
    double beta0 = priorDistribution.GetRate();
    double newAlpha = alpha * sample.size() + alpha0;
    double newBeta = sampleSum(sample) + beta0;
    GammaRand posteriorDistribution(newAlpha, newBeta);
    SetParameters(alpha, posteriorDistribution.Mean());
    return posteriorDistribution;
}


/// CHI-SQUARED
ChiSquaredRand::ChiSquaredRand(int degree)
{
    SetDegree(degree);
}

std::string ChiSquaredRand::Name() const
{
    return "Chi-squared(" + toStringWithPrecision(GetDegree()) + ")";
}

void ChiSquaredRand::SetDegree(int degree)
{
    GammaRand::SetParameters((degree < 1) ? 0.5 : 0.5 * degree, 0.5);
}


/// ERLANG
ErlangRand::ErlangRand(int shape, double rate)
{
    GammaRand::SetParameters(std::max(shape, 1), rate);
}

std::string ErlangRand::Name() const
{
    return "Erlang(" + toStringWithPrecision(GetShape()) + ", " + toStringWithPrecision(GetRate()) + ")";
}

int ErlangRand::GetShape() const
{
    return static_cast<int>(GammaRand::GetShape());
}

double ErlangRand::GetRate() const
{
    return 1.0 / GammaRand::GetScale();
}
