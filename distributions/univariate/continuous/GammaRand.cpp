#include "GammaRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"
#include "NormalRand.h"

GammaRand::GammaRand(double shape, double rate)
{
    setParameters(shape, rate);
}

std::string GammaRand::name() const
{
    return "Gamma(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

void GammaRand::setConstantsForGenerator()
{
    t = 0.07 + 0.75 * std::sqrt(1.0 - alpha);
    b = 1.0 + std::exp(-t) * alpha / t;
    alphaInv = 1.0 / alpha;
}

void GammaRand::setParameters(double shape, double rate)
{
    alpha = shape;
    if (alpha <= 0)
        alpha = 1.0;
    
    beta = rate;
    if (beta <= 0)
        beta = 1.0;
    theta = 1.0 / beta;

    mLgammaShape = -std::lgamma(alpha);
    pdfCoef = mLgammaShape + alpha * std::log(beta);

    if (getIdOfUsedGenerator(alpha) == SMALL_SHAPE)
        setConstantsForGenerator();
}

double GammaRand::f(double x) const
{
    if (x < 0)
        return 0;
    if (x == 0)
    {
        if (alpha > 1)
            return 0.0;
        return (alpha == 1) ? beta : INFINITY;
    }
    double y = (alpha - 1) * std::log(x);
    y -= x * beta;
    y += pdfCoef;
    return std::exp(y);
}

double GammaRand::F(double x) const
{
    if (x <= 0)
        return 0;
    double y = mLgammaShape + RandMath::logLowerIncGamma(alpha, x * beta);
    return std::exp(y);
}

GammaRand::GENERATOR_ID GammaRand::getIdOfUsedGenerator(double shape)
{
    if (shape < 0.34)
        return SMALL_SHAPE;
    if (shape <= 3 && RandMath::areClose(shape, std::round(shape)))
        return INTEGER_SHAPE;
    if (RandMath::areClose(shape, 1.5))
        return ONE_AND_A_HALF_SHAPE;
    if (shape > 1 && shape < 1.2)
        return FISHMAN;
    return MARSAGLIA_TSANG;
}

double GammaRand::variateThroughExponentialSum(int shape)
{
    double X = 0;
    for (int i = 0; i < shape; ++i)
        X += ExponentialRand::standardVariate();
    return X;
}

double GammaRand::variateForShapeOneAndAHalf()
{
    double W = ExponentialRand::standardVariate();
    double N = NormalRand::standardVariate();
    return W + 0.5 * N * N;
}

double GammaRand::variateBest() const
{
    /// Algorithm RGS for gamma variates (Best, 1983)
    double X = 0;
    int iter = 0;
    do {
        double V = b * UniformRand::standardVariate();
        double W = UniformRand::standardVariate();
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
        double U = UniformRand::standardVariate();
        double p = shape * t * U;
        double W = ExponentialRand::standardVariate();
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
        W1 = ExponentialRand::standardVariate();
        W2 = ExponentialRand::standardVariate();
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
            N = NormalRand::standardVariate();
        } while (N <= -c);
        double v = 1 + N / c;
        v = v * v * v;
        N *= N;
        double U = UniformRand::standardVariate();
        if (U < 1.0 - 0.331 * N * N || std::log(U) < 0.5 * N + d * (1.0 - v + std::log(v))) {
            return d * v;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// shouldn't end up here
}

double GammaRand::standardVariate(double shape)
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
    return NAN;
}

double GammaRand::variate(double shape, double rate)
{
    return standardVariate(shape) / rate;
}

double GammaRand::variate() const
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
    return NAN;
}

void GammaRand::sample(std::vector<double> &outputData) const
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

double GammaRand::Mean() const
{
    return alpha * theta;
}

double GammaRand::Variance() const
{
    return alpha * theta * theta;
}

std::complex<double> GammaRand::CF(double t) const
{
    return std::pow(std::complex<double>(1.0, -theta * t), -alpha);
}

double GammaRand::QuantileImpl(double p) const
{
    double root = p * Mean() / (1 - p); /// good starting point
    if (RandMath::findRoot([this, p] (double x)
    {
        double first = F(x) - p;
        double second = f(x);
        return DoublePair(first, second);
    }, root))
        return root;
    /// if we can't find quantile, then probably p -> 1
    return INFINITY;
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

bool GammaRand::fitScaleMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setParameters(alpha, alpha / sampleMean(sample));
    return true;
}

bool GammaRand::fitMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0 || !checkValidity(sample))
        return false;

    /// Calculate averages
    double average = sampleMean(sample);
    long double logAverage = 0.0L;
    for (double var : sample) {
        logAverage += std::log(var);
    }
    logAverage /= n;

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
        double second = 1.0 / x - RandMath::trigamma(x);
        return DoublePair(first, second);
    }, shape))
        return false;

    setParameters(shape, shape / average);
    return true;
}

bool GammaRand::fitShapeMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setParameters(sampleMean(sample) * beta, beta);
    return true;
}

bool GammaRand::fitScaleMM(const std::vector<double> &sample)
{
    return fitScaleMLE(sample);
}

bool GammaRand::fitMM(const std::vector<double> &sample)
{  
    if (!checkValidity(sample))
        return false;
    double mu1 = sampleMean(sample);
    double var = sampleVariance(sample, mu1);
    double shape = mu1 * mu1 / var;
    
    setParameters(shape, shape / mu1);
    return true;
}

bool GammaRand::fitRateBayes(const std::vector<double> &sample, GammaRand &priorDistribution)
{
    int n = sample.size();
    if (n <= 0 || !checkValidity(sample))
        return false;
    double alpha0 = priorDistribution.getShape();
    double beta0 = priorDistribution.getRate();
    double newAlpha = alpha * n + alpha0;
    double newBeta = sampleSum(sample) + beta0;
    priorDistribution.setParameters(newAlpha, newBeta);
    setParameters(alpha, priorDistribution.Mean());
    return true;
}


/// CHI-SQUARED
ChiSquaredRand::ChiSquaredRand(int degree)
{
    setDegree(degree);
}

std::string ChiSquaredRand::name() const
{
    return "Chi-squared(" + toStringWithPrecision(getDegree()) + ")";
}

void ChiSquaredRand::setDegree(int degree)
{
    GammaRand::setParameters((degree < 1) ? 0.5 : 0.5 * degree, 0.5);
}


/// ERLANG
ErlangRand::ErlangRand(int shape, double rate)
{
    GammaRand::setParameters(std::max(shape, 1), rate);
}

std::string ErlangRand::name() const
{
    return "Erlang(" + toStringWithPrecision(getShape()) + ", " + toStringWithPrecision(getRate()) + ")";
}

int ErlangRand::getShape() const
{
    return static_cast<int>(GammaRand::getShape());
}

double ErlangRand::getRate() const
{
    return 1.0 / GammaRand::getScale();
}
