#include "GammaRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

template < typename RealType >
GammaDistribution<RealType>::GammaDistribution(double shape, double rate)
{
    SetParameters(shape, rate);
}

template < typename RealType >
void GammaDistribution<RealType>::SetParameters(double shape, double rate)
{
    if (shape <= 0.0)
        throw std::invalid_argument("Gamma distribution: shape should be positive");
    if (rate <= 0.0)
        throw std::invalid_argument("Gamma distribution: rate should be positive");

    this->alpha = shape > 0 ? shape : 1.0;
    
    this->beta = (rate > 0.0) ? rate : 1.0;
    theta = 1.0 / this->beta;

    lgammaAlpha = std::lgammal(this->alpha);
    logAlpha = std::log(this->alpha);
    this->logBeta = std::log(this->beta);
    pdfCoef = -lgammaAlpha + this->alpha * this->logBeta;

    if (getIdOfUsedGenerator(this->alpha) == SMALL_SHAPE) {
        /// set constants for generator
        genCoef.t = 0.5 * std::log1pl(-alpha);
        genCoef.t = 0.07 + 0.75 * std::exp(genCoef.t);
        genCoef.b = 1.0 + std::exp(-genCoef.t) * this->alpha / genCoef.t;
    }
}

template < typename RealType >
void GammaDistribution<RealType>::SetShape(double shape)
{
    SetParameters(shape, this->beta);
}

template < typename RealType >
double GammaDistribution<RealType>::f(const RealType & x) const
{
    if (x < 0.0)
        return 0.0;
    if (x == 0.0)
    {
        if (this->alpha > 1.0)
            return 0.0;
        return (this->alpha == 1.0) ? this->beta : INFINITY;
    }
    return std::exp(logf(x));
}

template < typename RealType >
double GammaDistribution<RealType>::logf(const RealType & x) const
{
    if (x < 0.0)
        return -INFINITY;
    if (x == 0.0)
    {
        if (this->alpha > 1.0)
            return -INFINITY;
        return (this->alpha == 1.0) ? this->logBeta : INFINITY;
    }
    double y = (this->alpha - 1.0) * std::log(x);
    y -= x * this->beta;
    y += pdfCoef;
    return y;
}

template < typename RealType >
double GammaDistribution<RealType>::F(const RealType & x) const
{
    return (x > 0.0) ? RandMath::pgamma(this->alpha, x * this->beta, logAlpha, lgammaAlpha) : 0.0;
}

template < typename RealType >
double GammaDistribution<RealType>::logF(const RealType &x) const
{
    return (x > 0.0) ? RandMath::lpgamma(this->alpha, x * this->beta, logAlpha, lgammaAlpha) : -INFINITY;
}

template < typename RealType >
double GammaDistribution<RealType>::S(const RealType & x) const
{
    return (x > 0.0) ? RandMath::qgamma(this->alpha, x * this->beta, logAlpha, lgammaAlpha) : 1.0;
}

template < typename RealType >
double GammaDistribution<RealType>::logS(const RealType &x) const
{
    return (x > 0.0) ? RandMath::lqgamma(this->alpha, x * this->beta, logAlpha, lgammaAlpha) : 0.0;
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateThroughExponentialSum(int shape, RandGenerator& randGenerator)
{
    double X = 0.0;
    for (int i = 0; i < shape; ++i)
        X += ExponentialRand<RealType>::StandardVariate(randGenerator);
    return X;
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateForShapeOneAndAHalf(RandGenerator& randGenerator)
{
    RealType W = ExponentialRand<RealType>::StandardVariate(randGenerator);
    RealType N = NormalRand<RealType>::StandardVariate(randGenerator);
    return W + 0.5 * N * N;
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateBest(RandGenerator &randGenerator) const
{
    /// Algorithm RGS for gamma variates (Best, 1983)
    double X = 0;
    size_t iter = 0;
    do {
        double V = genCoef.b * UniformRand<RealType>::StandardVariate(randGenerator);
        double W = UniformRand<RealType>::StandardVariate(randGenerator);
        if (V <= 1) {
            X = genCoef.t * std::pow(V, 1.0 / this->alpha);
            if (W <= (2.0 - X) / (2.0 + X) || W <= std::exp(-X))
                return X;
        }
        else {
            X = -std::log(genCoef.t * (genCoef.b - V) / this->alpha);
            double Y = X / genCoef.t;
            if (W * (this->alpha + Y - this->alpha * Y) <= 1 || W <= std::pow(Y, this->alpha - 1))
                return X;
        }
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateAhrensDieter(double shape, RandGenerator &randGenerator)
{
    /// Rejection algorithm GS for gamma variates (Ahrens and Dieter, 1974)
    double X = 0;
    size_t iter = 0;
    double shapeInv = 1.0 / shape;
    double t = shapeInv + M_1_E;
    do {
        double U = UniformRand<RealType>::StandardVariate(randGenerator);
        double p = shape * t * U;
        double W = ExponentialRand<RealType>::StandardVariate(randGenerator);
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
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateFishman(double shape, RandGenerator& randGenerator)
{
    /// G. Fishman algorithm (shape > 1)
    double W1, W2;
    double shapem1 = shape - 1;
    do {
        W1 = ExponentialRand<RealType>::StandardVariate(randGenerator);
        W2 = ExponentialRand<RealType>::StandardVariate(randGenerator);
    } while (W2 < shapem1 * (W1 - std::log(W1) - 1));
    return shape * W1;
}

template < typename RealType >
RealType GammaDistribution<RealType>::variateMarsagliaTsang(double shape, RandGenerator &randGenerator)
{
    /// Marsaglia and Tsang’s Method (shape > 1/3)
    RealType d = shape - 1.0 / 3;
    RealType c = 3 * std::sqrt(d);
    size_t iter = 0;
    do {
        RealType N;
        do {
            N = NormalRand<RealType>::StandardVariate(randGenerator);
        } while (N <= -c);
        RealType v = 1 + N / c;
        v = v * v * v;
        N *= N;
        RealType U = UniformRand<RealType>::StandardVariate(randGenerator);
        if (U < 1.0 - 0.331 * N * N || std::log(U) < 0.5 * N + d * (1.0 - v + std::log(v))) {
            return d * v;
        }
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

template < typename RealType >
RealType GammaDistribution<RealType>::StandardVariate(double shape, RandGenerator& randGenerator)
{
    if (shape <= 0)
        return NAN;

    GENERATOR_ID genId = getIdOfUsedGenerator(shape);

    switch(genId) {
    case INTEGER_SHAPE:
        return variateThroughExponentialSum(std::round(shape), randGenerator);
    case ONE_AND_A_HALF_SHAPE:
        return variateForShapeOneAndAHalf(randGenerator);
    case SMALL_SHAPE:
        return variateAhrensDieter(shape, randGenerator);
    case FISHMAN:
        return variateFishman(shape, randGenerator);
    case MARSAGLIA_TSANG:
        return variateMarsagliaTsang(shape, randGenerator);
    default:
        return NAN;
    }
}

template < typename RealType >
RealType GammaDistribution<RealType>::Variate(double shape, double rate, RandGenerator& randGenerator)
{
    return (shape <= 0.0 || rate <= 0.0) ? NAN : StandardVariate(shape, randGenerator) / rate;
}

template < typename RealType >
RealType GammaDistribution<RealType>::Variate() const
{
    GENERATOR_ID genId = getIdOfUsedGenerator(this->alpha);

    switch(genId) {
    case INTEGER_SHAPE:
        return theta * variateThroughExponentialSum(this->alpha, this->localRandGenerator);
    case ONE_AND_A_HALF_SHAPE:
        return theta * variateForShapeOneAndAHalf(this->localRandGenerator);
    case SMALL_SHAPE:
        return theta * variateBest(this->localRandGenerator);
    case FISHMAN:
        return theta * variateFishman(this->alpha, this->localRandGenerator);
    case MARSAGLIA_TSANG:
        return theta * variateMarsagliaTsang(this->alpha, this->localRandGenerator);
    default:
        return NAN;
    }
}

template < typename RealType >
void GammaDistribution<RealType>::Sample(std::vector<RealType> &outputData) const
{
    GENERATOR_ID genId = getIdOfUsedGenerator(this->alpha);

    switch(genId) {
    case INTEGER_SHAPE:
        for (RealType &var : outputData)
            var = theta * variateThroughExponentialSum(this->alpha, this->localRandGenerator);
        break;
    case ONE_AND_A_HALF_SHAPE:
        for (RealType &var : outputData)
            var = theta * variateForShapeOneAndAHalf(this->localRandGenerator);
        break;
    case SMALL_SHAPE:
        for (RealType &var : outputData)
            var = theta * variateBest(this->localRandGenerator);
        break;
    case FISHMAN:
        for (RealType &var : outputData)
            var = theta * variateFishman(this->alpha, this->localRandGenerator);
        break;
    case MARSAGLIA_TSANG:
        for (RealType &var : outputData)
            var = theta * variateMarsagliaTsang(this->alpha, this->localRandGenerator);
        break;
    default:
        return;
    }
}

template < typename RealType >
long double GammaDistribution<RealType>::Mean() const
{
    return this->alpha * theta;
}

template < typename RealType >
long double GammaDistribution<RealType>::GeometricMean() const
{
    return RandMath::digamma(this->alpha) - this->logBeta;
}

template < typename RealType >
long double GammaDistribution<RealType>::Variance() const
{
    return this->alpha * theta * theta;
}

template < typename RealType >
long double GammaDistribution<RealType>::GeometricVariance() const
{
    return RandMath::trigamma(this->alpha);
}

template < typename RealType >
RealType GammaDistribution<RealType>::Mode() const
{
    return (this->alpha <= 1) ? 0 : (this->alpha - 1) * theta;
}

template < typename RealType >
RealType GammaDistribution<RealType>::Median() const
{
    return (this->alpha == 1.0) ? theta * M_LN2 : quantileImpl(0.5);
}

template < typename RealType >
long double GammaDistribution<RealType>::Skewness() const
{
    return 2.0l / std::sqrt(this->alpha);
}

template < typename RealType >
long double GammaDistribution<RealType>::ExcessKurtosis() const
{
    return 6.0l / this->alpha;
}

template < typename RealType >
RealType GammaDistribution<RealType>::initRootForSmallP(double r) const
{
    double root = 0;
    double c[5];
    c[4] = 1;
    /// first coefficient
    double denominator = this->alpha + 1;
    c[3] = 1.0 / denominator;
    /// second coefficient
    denominator *= denominator;
    denominator *= this->alpha + 2;
    c[2] = 0.5 * (3 * this->alpha + 5) / denominator;
    /// third coefficient
    denominator *= 3 * (this->alpha + 1) * (this->alpha + 3);
    c[1] = 8 * this->alpha + 33;
    c[1] *= this->alpha;
    c[1] += 31;
    c[1] /= denominator;
    /// fourth coefficient
    denominator *= 8 * (this->alpha + 1) * (this->alpha + 2) * (this->alpha + 4);
    c[0] = 125 * this->alpha + 1179;
    c[0] *= this->alpha;
    c[0] += 3971;
    c[0] *= this->alpha;
    c[0] += 5661;
    c[0] *= this->alpha;
    c[0] += 2888;
    c[0] /= denominator;
    /// now calculate root
    for (int i = 0; i != 5; ++i) {
        root += c[i];
        root *= r;
    }
    return root;
}

template < typename RealType >
RealType GammaDistribution<RealType>::initRootForLargeP(double logQ) const
{
    /// look for approximate value of x -> INFINITY
    double x = (logQ + lgammaAlpha) / this->alpha;
    x = -std::exp(x) / this->alpha;
    return -alpha * RandMath::Wm1Lambert(x);
}

template < typename RealType >
RealType GammaDistribution<RealType>::initRootForLargeShape(double p) const
{
    if (p == 0.5)
        return this->alpha;
    double x = RandMath::erfcinv(2 * p);
    double lambda = x * x / this->alpha + 1;
    lambda = -std::exp(-lambda);
    lambda = (x < 0) ? -RandMath::Wm1Lambert(lambda) : -RandMath::W0Lambert(lambda);
    return lambda * this->alpha;
}

template < typename RealType >
RealType GammaDistribution<RealType>::initRootForLargeShape1m(double p) const
{
    if (p == 0.5)
        return this->alpha;
    double x = -RandMath::erfcinv(2 * p);
    double lambda = x * x / this->alpha + 1;
    lambda = -std::exp(-lambda);
    lambda = (x < 0) ? -RandMath::Wm1Lambert(lambda) : -RandMath::W0Lambert(lambda);
    return lambda * this->alpha;
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileInitialGuess(double p) const
{
    /// Method is taken from
    /// "Efficient and accurate algorithms
    /// for the computation and inversion
    /// of the incomplete gamma function ratios"
    /// (Amparo Gil, Javier Segura and Nico M. Temme)
    double guess = 0;
    if (this->alpha < 10) {
        double r = std::log(p * this->alpha) + lgammaAlpha;
        r = std::exp(r / this->alpha);
        /// if p -> 0
        if (r < 0.2 * (this->alpha + 1)) {
            guess = initRootForSmallP(r);
        }
        else {
            double logQ = std::log1pl(-p);
            /// boundary adviced in a paper
            double maxBoundary1 = -0.5 * this->alpha - logAlpha - lgammaAlpha;
            /// the maximum possible value to have a solution
            double maxBoundary2 = this->alpha * (logAlpha - 1) - lgammaAlpha;
            /// if p -> 1
            if (logQ < std::min(maxBoundary1, maxBoundary2))
                guess = initRootForLargeP(logQ);
            else if (this->alpha < 1)
                guess = r;
            else
                guess = initRootForLargeShape(p);
        }
    }
    else
        guess = initRootForLargeShape(p);
    return guess / this->beta;
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileInitialGuess1m(double p) const
{
    if (this->alpha < 10) {
        double logQ = std::log(p);
        /// boundary adviced in a paper
        double maxBoundary1 = -0.5 * this->alpha - logAlpha - lgammaAlpha;
        /// the maximum possible value to have a solution
        double maxBoundary2 = this->alpha * (logAlpha - 1) - lgammaAlpha;
        /// if p -> 0
        if (logQ < std::min(maxBoundary1, maxBoundary2))
            return initRootForLargeP(logQ) / this->beta;
    }
    else {
        return initRootForLargeShape1m(p) / this->beta;
    }
    return quantileInitialGuess(1.0 - p);
}

template < typename RealType >
double GammaDistribution<RealType>::df(RealType x) const
{
    double z = (this->alpha - 1) - this->beta * x;
    double y = (this->alpha - 2) * std::log(x);
    y -= this->beta * x;
    y += pdfCoef;
    return z * std::exp(y);
}

template < typename RealType >
double GammaDistribution<RealType>::dfDivf(RealType x) const
{
    return x / (this->alpha - 1 - this->beta * x);
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileImpl(double p, RealType initValue) const
{
    if (p < 1e-5) { /// too small p
        double logP = std::log(p);
        if (RandMath::findRoot<RealType>([this, logP] (double x)
        {
            if (x <= 0)
               return DoubleTriplet(-INFINITY, 0, 0);
            double logCdf = logF(x), logPdf = logf(x);
            double first = logCdf - logP;
            double second = std::exp(logPdf - logCdf);
            double third = second * (dfDivf(x) - second);
            return DoubleTriplet(first, second, third);
        }, initValue))
            return initValue;
        /// if we can't find quantile, then probably something bad has happened
        return NAN;
    }
    if (RandMath::findRoot<RealType>([this, p] (double x)
    {
        if (x <= 0)
            return DoubleTriplet(-p, 0, 0);
        double first = F(x) - p;
        double second = f(x);
        double third = df(x);
        return DoubleTriplet(first, second, third);
    }, initValue))
        return initValue;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileImpl(double p) const
{
    return (this->alpha == 1.0) ? -theta * std::log1pl(-p) : quantileImpl(p, quantileInitialGuess(p));
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileImpl1m(double p, RealType initValue) const
{
    if (p < 1e-5) { /// too small p
        double logP = std::log(p);
        if (RandMath::findRoot<RealType>([this, logP] (double x)
        {
           if (x <= 0)
               return DoubleTriplet(logP, 0, 0);
            double logCcdf = logS(x), logPdf = logf(x);
            double first = logP - logCcdf;
            double second = std::exp(logPdf - logCcdf);
            double third = second * (dfDivf(x) + second);
            return DoubleTriplet(first, second, third);
        }, initValue))
            return initValue;
        /// if we can't find quantile, then probably something bad has happened
        return NAN;
    }
    if (RandMath::findRoot<RealType>([this, p] (double x)
    {
        if (x <= 0)
            return DoubleTriplet(p - 1.0, 0, 0);
        double first = p - S(x);
        double second = f(x);
        double third = df(x);
        return DoubleTriplet(first, second, third);
    }, initValue))
        return initValue;
    /// if we can't find quantile, then probably something bad has happened
    return NAN;
}

template < typename RealType >
RealType GammaDistribution<RealType>::quantileImpl1m(double p) const
{
    return (this->alpha == 1.0) ? -theta * std::log(p) : quantileImpl1m(p, quantileInitialGuess1m(p));
}

template < typename RealType >
std::complex<double> GammaDistribution<RealType>::CFImpl(double t) const
{
    return std::pow(std::complex<double>(1.0, -theta * t), -alpha);
}

template class GammaDistribution<float>;
template class GammaDistribution<double>;
template class GammaDistribution<long double>;


// FREE-SCALE-GAMMA-DISTRIBUTION

template < typename RealType >
void FreeScaleGammaDistribution<RealType>::SetRate(double rate)
{
    this->SetParameters(this->alpha, rate);
}

template < typename RealType >
void FreeScaleGammaDistribution<RealType>::SetScale(double scale)
{
    SetRate(1.0 / scale);
}

template < typename RealType >
void FreeScaleGammaDistribution<RealType>::FitRate(const std::vector<RealType> &sample, bool unbiased)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    double mean = this->GetSampleMean(sample);
    double coef = this->alpha - (unbiased ? 1.0 / sample.size() : 0.0);
    this->SetParameters(this->alpha, coef / mean);
}

template < typename RealType >
GammaRand<RealType> FreeScaleGammaDistribution<RealType>::FitRateBayes(const std::vector<RealType> &sample, const GammaDistribution<RealType> &priorDistribution, bool MAP)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));
    double kappa = priorDistribution.GetShape();
    double gamma = priorDistribution.GetRate();
    double newShape = this->alpha * sample.size() + kappa;
    double newRate = this->GetSampleSum(sample) + gamma;
    GammaRand<RealType> posteriorDistribution(newShape, newRate);
    this->SetParameters(this->alpha, MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean());
    return posteriorDistribution;
}

template class FreeScaleGammaDistribution<float>;
template class FreeScaleGammaDistribution<double>;
template class FreeScaleGammaDistribution<long double>;


// GAMMARAND

template < typename RealType >
String GammaRand<RealType>::Name() const
{
    return "Gamma(" + this->toStringWithPrecision(this->GetShape()) + ", " + this->toStringWithPrecision(this->GetRate()) + ")";
}

template < typename RealType >
void GammaRand<RealType>::FitShape(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));

    /// Calculate initial guess via method of moments
    double shape = this->GetSampleMean(sample) * this->beta;
    /// Run root-finding procedure
    double s = this->GetSampleLogMean(sample) + this->logBeta;
    if (!RandMath::findRoot<double>([s] (double x)
    {
        double first = RandMath::digamma(x) - s;
        double second = RandMath::trigamma(x);
        return DoublePair(first, second);
    }, shape))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure"));
    SetParameters(shape, this->beta);
}

template < typename RealType >
void GammaRand<RealType>::Fit(const std::vector<RealType> &sample)
{
    /// Sanity check
    if (!this->allElementsArePositive(sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->POSITIVITY_VIOLATION));

    /// Calculate initial guess for shape
    double average = this->GetSampleMean(sample);
    double s = std::log(average) - this->GetSampleLogMean(sample);
    double sm3 = s - 3.0, sp12 = 12.0 * s;
    double shape = sm3 * sm3 + 2 * sp12;
    shape = std::sqrt(shape);
    shape -= sm3;
    shape /= sp12;

    if (!RandMath::findRoot<double>([s] (double x)
    {
        double first = RandMath::digammamLog(x) + s;
        double second = RandMath::trigamma(x) - 1.0 / x;
        return DoublePair(first, second);
    }, shape))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure"));

    SetParameters(shape, shape / average);
}

template class GammaRand<float>;
template class GammaRand<double>;
template class GammaRand<long double>;


// CHI-SQUAREDRAND

template < typename RealType >
String ChiSquaredRand<RealType>::Name() const
{
    return "Chi-squared(" + this->toStringWithPrecision(GetDegree()) + ")";
}

template < typename RealType >
void ChiSquaredRand<RealType>::SetDegree(size_t degree)
{
    GammaDistribution<RealType>::SetParameters(0.5 * degree, 0.5);
}

template class ChiSquaredRand<float>;
template class ChiSquaredRand<double>;
template class ChiSquaredRand<long double>;


// ERLANGRAND

template < typename RealType >
String ErlangRand<RealType>::Name() const
{
    return "Erlang(" + this->toStringWithPrecision(this->GetShape()) + ", " + this->toStringWithPrecision(this->GetRate()) + ")";
}

template < typename RealType >
void ErlangRand<RealType>::SetParameters(size_t shape, double rate)
{
    GammaDistribution<RealType>::SetParameters(shape, rate);
}

template < typename RealType >
void ErlangRand<RealType>::SetShape(size_t shape)
{
    SetParameters(shape, this->beta);
}

template class ErlangRand<float>;
template class ErlangRand<double>;
template class ErlangRand<long double>;
