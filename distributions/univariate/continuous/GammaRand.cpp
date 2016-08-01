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
    m = alpha - 1;
    static constexpr double M_SQRT8_3 = 1.63299316185545206546;
    double sqrtAlpha = std::sqrt(alpha);
    s_2 = M_SQRT8_3 * sqrtAlpha + alpha;
    s = std::sqrt(s_2);
    d = M_SQRT2 * M_SQRT3 * s_2;
    b = d + m;
    w = s_2 / (m - 1);
    v = (s_2 + s_2) / (m * sqrtAlpha);
    c = b + std::log(s * d / b) - m - m - 3.7203285;
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

    if (getIdOfUsedGenerator(alpha) == LARGE_SHAPE)
        setConstantsForGenerator();
}

double GammaRand::f(double x) const
{
    if (x < 0)
        return 0;
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
    if (shape < 1)
        return SMALL_SHAPE;
    if (shape <= 5) {
        if (RandMath::areClose(shape, std::round(shape)))
            return INTEGER_SHAPE;
        if (shape <= 3.5) {
            double shapeMHalf = shape - 0.5;
            if (RandMath::areClose(shapeMHalf, std::round(shapeMHalf)))
                return HALF_INTEGER_SHAPE;
            if (shape < 3)
                return MEDIUM_SHAPE;
        }
    }
    return LARGE_SHAPE;
}

double GammaRand::variateForIntegerShape(int shape)
{
    double rv = 0;
    for (int i = 0; i < shape; ++i)
        rv += ExponentialRand::standardVariate();
    return rv;
}

double GammaRand::variateForHalfIntegerShape(int shape)
{
    double rv = 0;
    for (int i = 0; i < shape; ++i)
        rv += ExponentialRand::standardVariate();
    double N = NormalRand::standardVariate();
    return rv + .5 * N * N;
}

double GammaRand::variateForSmallShape(double shape)
{
    double rv = 0;
    int iter = 0;
    double shapeInv = 1.0 / shape;
    double t = shapeInv + M_1_E;
    do {
        double U = UniformRand::standardVariate();
        double p = shape * t * U;
        double W = ExponentialRand::standardVariate();
        if (p <= 1)
        {
            rv = std::pow(p, shapeInv);
            if (rv <= W)
                return rv;
        }
        else
        {
            rv = -std::log(t * (1 - U));
            if ((1 - shape) * std::log(rv) <= W)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return NAN; /// shouldn't end up here
}

double GammaRand::variateForMediumShape(double shape)
{
    double W1, W2;
    double shapem1 = shape - 1;
    do {
        W1 = ExponentialRand::standardVariate();
        W2 = ExponentialRand::standardVariate();
    } while (W2 < shapem1 * (W1 - std::log(W1) - 1));
    return shape * W1;
}

double GammaRand::variateForLargeShape() const
{
    double rv = 0;
    int iter = 0;
    do {
        double U = UniformRand::standardVariate();
        if (U <= 0.0095722652)
        {
            double W1 = ExponentialRand::standardVariate();
            double W2 = ExponentialRand::standardVariate();
            rv = b * (1 + W1 / d);
            if (m * (rv / b - std::log(rv / m)) + c <= W2)
                return rv;
        }
        else
        {
            double N;
            do {
                N = NormalRand::standardVariate();
                rv = s * N + m;
            } while (rv < 0 || rv > b);
            U = UniformRand::standardVariate();
            double S = .5 * N * N;
            if (N > 0) {
                if (U < 1 - w * S)
                    return rv;
            }
            else if (U < 1 + S * (v * N - w))
                return rv;

            if (std::log(U) < m * std::log(rv / m) + m - rv + S)
                return rv;
        }
    } while (++iter < 1e9); /// one billion should be enough
    return NAN; /// shouldn't end up here
}

double GammaRand::variateForLargeShape(double shape)
{
    /// Marsaglia and Tsang’s Method for shape > 1
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
    } while (++iter < 1e9); /// one billion should be enough
    return NAN; /// shouldn't end up here
}

double GammaRand::standardVariate(double shape)
{
    GENERATOR_ID genId = getIdOfUsedGenerator(shape);

    switch(genId) {
    case INTEGER_SHAPE:
        return variateForIntegerShape(std::round(shape));
    case HALF_INTEGER_SHAPE:
        return variateForHalfIntegerShape(std::round(shape));
    case SMALL_SHAPE:
        return variateForSmallShape(shape);
    case MEDIUM_SHAPE:
        return variateForMediumShape(shape);
    case LARGE_SHAPE:
        return variateForLargeShape(shape);
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
        return theta * variateForIntegerShape(alpha);
    case HALF_INTEGER_SHAPE:
        return theta * variateForHalfIntegerShape(alpha);
    case SMALL_SHAPE:
        return theta * variateForSmallShape(alpha);
    case MEDIUM_SHAPE:
        return theta * variateForMediumShape(alpha);
    case LARGE_SHAPE:
        return theta * variateForLargeShape();
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
            var = theta * variateForIntegerShape(alpha);
        break;
    case HALF_INTEGER_SHAPE:
        for (double &var : outputData)
            var = theta * variateForHalfIntegerShape(alpha);
        break;
    case SMALL_SHAPE:
        for (double &var : outputData)
            var = theta * variateForSmallShape(alpha);
        break;
    case MEDIUM_SHAPE:
        for (double &var : outputData)
            var = theta * variateForMediumShape(alpha);
        break;
    case LARGE_SHAPE:
        for (double &var : outputData)
            var = theta * variateForLargeShape();
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

bool GammaRand::checkValidity(const std::vector<double> &sample)
{
    for (double var : sample) {
        if (var < 0)
            return false;
    }
    return true;
}

bool GammaRand::fitScaleMLE(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setParameters(alpha, alpha / RandMath::sampleMean(sample));
    return true;
}

bool GammaRand::fitMLE(const std::vector<double> &sample)
{
    int n = sample.size();
    if (n <= 0 || !checkValidity(sample))
        return false;

    /// Calculate averages
    double average = RandMath::sampleMean(sample);
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
        return std::log(x) - RandMath::digamma(x) - s;
    },
    [] (double x)
    {
        return 1.0 / x - RandMath::trigamma(x);
    }, shape))
        return false;

    setParameters(shape, shape / average);
    return true;
}

bool GammaRand::fitShapeMM(const std::vector<double> &sample)
{
    if (!checkValidity(sample))
        return false;
    setParameters(RandMath::sampleMean(sample) * beta, beta);
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
    double mu1 = RandMath::sampleMean(sample);
    double var = RandMath::sampleVariance(sample, mu1);
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
    double newBeta = RandMath::sum(sample) + beta0;
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
