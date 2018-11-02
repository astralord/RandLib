#include "BetaRand.h"
#include "../discrete/BernoulliRand.h"
#include "UniformRand.h"
#include "ExponentialRand.h"

template < typename RealType >
BetaDistribution<RealType>::BetaDistribution(double shape1, double shape2, double minValue, double maxValue)
{
    SetShapes(shape1, shape2);
    SetSupport(minValue, maxValue);
}

template < typename RealType >
void BetaDistribution<RealType>::setCoefficientsForGenerator()
{
    GENERATOR_ID id = getIdOfUsedGenerator();
    if (id == REJECTION_NORMAL) {
        double alpham1 = alpha - 1;
        genCoef.s = alpham1 * std::log1pl(0.5 / alpham1) - 0.5;
        genCoef.t = 1.0 / std::sqrt(8 * alpha - 4);
    }
    else if (id == CHENG) {
        genCoef.s = alpha + beta;
        genCoef.t = std::min(alpha, beta);
        if (genCoef.t > 1)
            genCoef.t = std::sqrt((2 * alpha * beta - genCoef.s) / (genCoef.s - 2));
        genCoef.u = alpha + genCoef.t;
    }
    else if (id == ATKINSON_WHITTAKER) {
        genCoef.t = std::sqrt(alpha * (1 - alpha));
        genCoef.t /= (genCoef.t + std::sqrt(beta * (1 - beta)));
        genCoef.s = beta * genCoef.t;
        genCoef.s /= (genCoef.s + alpha * (1 - genCoef.t));
    }
}

template < typename RealType >
void BetaDistribution<RealType>::SetShapes(double shape1, double shape2)
{
    if (shape1 <= 0 || shape2 <= 0)
        throw std::invalid_argument("Beta distribution: shapes should be positive");
    GammaRV1.SetParameters(shape1, 1);
    GammaRV2.SetParameters(shape2, 1);
    alpha = GammaRV1.GetShape();
    beta = GammaRV2.GetShape();
    logBetaFun = -std::lgammal(alpha + beta) + GammaRV1.GetLogGammaShape() + GammaRV2.GetLogGammaShape();
    betaFun = std::exp(logBetaFun);
    setCoefficientsForGenerator();
}

template < typename RealType >
void BetaDistribution<RealType>::SetSupport(double minValue, double maxValue)
{
    if (minValue >= maxValue)
        throw std::invalid_argument("Beta distribution: minimum value should be smaller than maximum value");

    a = minValue;
    b = maxValue;
    bma = b - a;
    bmaInv = 1.0 / bma;
    logbma = std::log(bma);
}

template < typename RealType >
double BetaDistribution<RealType>::f(const RealType &x) const
{
    if (x < a || x > b)
        return 0.0;
    if (x == a) {
        if (alpha == 1)
            return beta / bma;
        return (alpha > 1) ? 0 : INFINITY;
    }
    if (x == b) {
        if (beta == 1)
            return alpha / bma;
        return (beta > 1) ? 0 : INFINITY;
    }
    return std::exp(logf(x));
}

template < typename RealType >
double BetaDistribution<RealType>::logf(const RealType &x) const
{
    /// Standardize
    double xSt = (x - a) / bma;
    if (xSt < 0.0 || xSt > 1.0)
        return -INFINITY;
    if (xSt == 0.0) {
        if (alpha == 1)
            return std::log(beta / bma);
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    if (xSt == 1.0) {
        if (beta == 1)
            return std::log(alpha / bma);
        return (beta > 1) ? -INFINITY : INFINITY;
    }
    double y = (alpha - 1) * std::log(xSt);
    y += (beta - 1) * std::log1pl(-xSt);
    return y - logBetaFun - logbma;
}

template < typename RealType >
double BetaDistribution<RealType>::F(const RealType &x) const
{
    if (x <= a)
        return 0.0;
    if (x >= b)
        return 1.0;
    /// Standardize
    double xSt = (x - a) / bma;
    /// Workaround known case
    if (alpha == beta && beta == 0.5)
        return M_2_PI * std::asin(std::sqrt(xSt));
    return RandMath::ibeta(xSt, alpha, beta, logBetaFun, std::log(xSt), std::log1pl(-xSt));
}

template < typename RealType >
double BetaDistribution<RealType>::S(const RealType &x) const
{
    if (x <= a)
        return 1.0;
    if (x >= b)
        return 0.0;
    /// Standardize
    double xSt = (x - a) / bma;
    /// Workaround known case
    if (alpha == beta && beta == 0.5)
        return M_2_PI * std::acos(std::sqrt(xSt));
    return RandMath::ibeta(1.0 - xSt, beta, alpha, logBetaFun, std::log1pl(-xSt), std::log(xSt));
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateArcsine() const
{
    RealType U = 2 * UniformRand<RealType>::StandardVariate(this->localRandGenerator) - 1;
    RealType X = std::sin(M_PI * U);
    return X * X;
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateRejectionUniform() const
{
    size_t iter = 0;
    do {
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType V = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        if (0.25 * V * V <= U - U * U)
            return U;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateRejectionUniformExtended() const
{
    size_t iter = 0;
    static constexpr double M_LN4 = M_LN2 + M_LN2;
    do {
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType edge = M_LN4 + std::log(U - U * U);
        if (W >= (1.0 - alpha) * edge)
            return U;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateRejectionNormal() const
{
    size_t iter = 0;
    RealType N = 0, Z = 0;
    RealType alpham1 = alpha - 1;
    RealType alpha2m1 = alpha + alpham1;
    do {
        do {
            N = NormalRand<RealType>::StandardVariate(this->localRandGenerator);
            Z = N * N;
        } while (Z >= alpha2m1);

        RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator) + genCoef.s;
        RealType aux = 0.5 - alpham1 / (alpha2m1 - Z);
        aux *= Z;
        if (W + aux >= 0)
            return 0.5 + N * genCoef.t;
        aux = std::log1pl(-Z / alpha2m1);
        aux *= alpham1;
        aux += W + 0.5 * Z;
        if (aux >= 0)
            return 0.5 + N * genCoef.t;
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateJohnk() const
{
    RealType X = 0, Z = 0;
    RealType W = 0, V = 0;
    do {
        W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator) / alpha;
        V = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator) / beta;
        X = std::exp(-W);
        Z = X + std::exp(-V);
    } while (Z > 1);
    return (Z > 0) ? (X / Z) : (W < V);;
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateCheng() const
{
    RealType R, T, Y;
    do {
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType V = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType X = std::log(U / (1 - U)) / genCoef.t;
        Y = alpha * std::exp(X);
        R = 1.0 / (beta + Y);
        T = 4 * U * U * V;
        T = std::log(T);
        T -= genCoef.u * X;
        T -= genCoef.s * std::log(genCoef.s * R);
    } while (T > 0);
    return Y * R;
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateAtkinsonWhittaker() const
{
    size_t iter = 0;
    do {
        RealType U = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        RealType W = ExponentialRand<RealType>::StandardVariate(this->localRandGenerator);
        if (U <= genCoef.s) {
            RealType X = genCoef.t * std::pow(U / genCoef.s, 1.0 / alpha);
            if (W >= (1.0 - beta) * std::log((1.0 - X) / (1.0 - genCoef.t)))
                return X;
        }
        else {
            RealType X = 1.0 - (1.0 - genCoef.t) * std::pow((1.0 - U) / (1.0 - genCoef.s), 1.0 / beta);
            if (W >= (1.0 - alpha) * std::log(X / genCoef.t))
                return X;
        }
    } while (++iter <= this->MAX_ITER_REJECTION);
    throw std::runtime_error("Beta distribution: sampling failed");
}

template < typename RealType >
RealType BetaDistribution<RealType>::variateGammaRatio() const
{
    RealType Y = GammaRV1.Variate();
    RealType Z = GammaRV2.Variate();
    return Y / (Y + Z);
}

template < typename RealType >
RealType BetaDistribution<RealType>::Variate() const
{
    double var = 0;
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM:
        var = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        break;
    case ARCSINE:
        var = variateArcsine();
        break;
    case CHENG:
        var = variateCheng();
        break;
    case REJECTION_UNIFORM:
        var = variateRejectionUniform();
        break;
    case REJECTION_UNIFORM_EXTENDED:
        var = variateRejectionUniformExtended();
        break;
    case REJECTION_NORMAL:
        var = variateRejectionNormal();
        break;
    case JOHNK:
        var = variateJohnk();
        break;
    case ATKINSON_WHITTAKER:
        var = variateAtkinsonWhittaker();
        break;
    case GAMMA_RATIO:
    default:
        var = variateGammaRatio();
        break;
    }

    return a + bma * var;
}

template < typename RealType >
void BetaDistribution<RealType>::Sample(std::vector<RealType> &outputData) const
{
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM: {
        for (RealType &var : outputData)
            var = UniformRand<RealType>::StandardVariate(this->localRandGenerator);
        }
        break;
    case ARCSINE: {
        for (RealType &var : outputData)
            var = variateArcsine();
        }
        break;
    case CHENG: {
        for (RealType &var : outputData)
            var = variateCheng();
        }
        break;
    case REJECTION_UNIFORM: {
        for (RealType &var : outputData)
            var = variateRejectionUniform();
        }
        break;
    case REJECTION_UNIFORM_EXTENDED: {
        for (RealType &var : outputData)
            var = variateRejectionUniformExtended();
        }
        break;
    case REJECTION_NORMAL: {
        for (RealType &var : outputData)
            var = variateRejectionNormal();
        }
        break;
    case JOHNK: {
        for (RealType &var : outputData)
            var = variateJohnk();
        }
        break;
    case ATKINSON_WHITTAKER: {
        for (RealType &var : outputData)
            var = variateAtkinsonWhittaker();
        }
        break;
    case GAMMA_RATIO:
    default: {
        GammaRV1.Sample(outputData);
        for (RealType &var : outputData)
            var /= (var + GammaRV2.Variate());
        }
        break;
    }

    /// Shift and scale
    for (RealType &var : outputData)
        var = a + bma * var;
}

template < typename RealType >
void BetaDistribution<RealType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    GammaRV1.Reseed(seed + 1);
    GammaRV2.Reseed(seed + 2);
}

template < typename RealType >
long double BetaDistribution<RealType>::Mean() const
{
    double mean = alpha / (alpha + beta);
    return a + bma * mean;
}

template < typename RealType >
long double BetaDistribution<RealType>::GeometricMean() const
{
    return RandMath::digamma(alpha) - RandMath::digamma(alpha + beta);
}

template < typename RealType >
long double BetaDistribution<RealType>::Variance() const
{
    double var = alpha + beta;
    var *= var * (var + 1);
    var = alpha * beta / var;
    return bma * bma * var;
}

template < typename RealType >
long double BetaDistribution<RealType>::GeometricVariance() const
{
    return RandMath::trigamma(alpha) - RandMath::trigamma(alpha + beta);
}

template < typename RealType >
RealType BetaDistribution<RealType>::Median() const
{
    if (alpha == beta)
        return a + bma * 0.5;
    if (alpha == 1.0)
        return a - bma * std::expm1l(-M_LN2 / beta);
    if (beta == 1.0)
        return a + bma * std::exp(-M_LN2 / alpha);
    if (alpha >= 1.0 && beta >= 1.0) {
        double initValue = 3 * alpha - 1.0;
        initValue /= 3 * (alpha + beta) - 2.0;
        initValue *= bma;
        initValue += a;
        return ContinuousDistribution<RealType>::quantileImpl(0.5, initValue);
    }
    return ContinuousDistribution<RealType>::quantileImpl(0.5);
}

template < typename RealType >
RealType BetaDistribution<RealType>::Mode() const
{
    double mode;
    if (alpha > 1)
        mode = (beta > 1) ? (alpha - 1) / (alpha + beta - 2) : 1.0;
    else
        mode = (beta > 1) ? 0.0 : (alpha > beta);
    return a + bma * mode;
}

template < typename RealType >
long double BetaDistribution<RealType>::Skewness() const
{
    long double skewness = (alpha + beta + 1) / (alpha * beta);
    skewness = std::sqrt(skewness);
    skewness *= beta - alpha;
    skewness /= alpha + beta + 2;
    return 2 * skewness;
}

template < typename RealType >
long double BetaDistribution<RealType>::ExcessKurtosis() const
{
    long double sum = alpha + beta;
    long double kurtosis = alpha - beta;
    kurtosis *= kurtosis;
    kurtosis *= (sum + 1);
    kurtosis /= (alpha * beta * (sum + 2));
    --kurtosis;
    kurtosis /= (sum + 3);
    return 6 * kurtosis;
}

template < typename RealType >
long double BetaDistribution<RealType>::MeanAbsoluteDeviation() const
{
    double y = M_LN2;
    y += alpha * std::log(alpha);
    y += beta * std::log(beta);
    y -= (alpha + beta + 1) * std::log(alpha + beta);
    y -= logBetaFun;
    y += logbma;
    return std::exp(y);
}

template < typename RealType >
RealType BetaDistribution<RealType>::quantileImpl(double p) const
{
    if (alpha == beta)
    {
        if (alpha == 0.5) {
            double x = std::sin(0.5 * M_PI * p);
            return a + bma * x * x;
        }
        if (alpha == 1.0)
            return a + bma * p;
    }
    if (alpha == 1.0)
        return a - bma * std::expm1l(std::log1pl(-p) / beta);
    if (beta == 1.0)
        return a + bma * std::pow(p, 1.0 / alpha);
    return ContinuousDistribution<RealType>::quantileImpl(p);
}

template < typename RealType >
RealType BetaDistribution<RealType>::quantileImpl1m(double p) const
{
    if (alpha == beta)
    {
        if (alpha == 0.5) {
            double x = std::cos(0.5 * M_PI * p);
            return a + bma * x * x;
        }
        if (alpha == 1.0)
            return b - bma * p;
    }
    if (alpha == 1.0)
        return a - bma * std::expm1l(std::log(p) / beta);
    if (beta == 1.0)
        return a + bma * std::exp(std::log1pl(-p) / alpha);
    return ContinuousDistribution<RealType>::quantileImpl1m(p);
}

template < typename RealType >
std::complex<double> BetaDistribution<RealType>::CFImpl(double t) const
{
    /// if we don't have singularity points, we can use direct integration
    if (alpha >= 1 && beta >= 1)
        return UnivariateDistribution<RealType>::CFImpl(t);

    double z = bma * t;
    double sinZ = std::sin(z);
    double cosZm1 = std::cos(z) - 1.0;

    double re = RandMath::integral([this, z, cosZm1](double x) {
        if (x >= 1)
            return 0.0;
        if (x <= 0)
            return -cosZm1;
        double y = std::cos(z * x) - 1;
        y *= std::pow(x, alpha - 1);
        y -= cosZm1;
        return std::pow(1.0 - x, beta - 1) * y;
    }, 0, 1);
    re += betaFun;
    re += cosZm1 / beta;

    double im = RandMath::integral([this, z, sinZ](double x) {
        if (x >= 1)
            return 0.0;
        if (x <= 0)
            return -sinZ;
        double y = std::sin(z * x);
        y *= std::pow(x, alpha - 1);
        y -= sinZ;
        return std::pow(1.0 - x, beta - 1) * y;
    }, 0, 1);
    im += sinZ / beta;

    std::complex<double> y(re, im);
    double cosTA = std::cos(t * a), sinTA = std::sin(t * a);
    return y * std::complex<double>(cosTA, sinTA) / betaFun;
}

template < typename RealType >
constexpr char BetaDistribution<RealType>::ALPHA_ZERO[];
template < typename RealType >
constexpr char BetaDistribution<RealType>::BETA_ZERO[];

template class BetaDistribution<float>;
template class BetaDistribution<double>;
template class BetaDistribution<long double>;

// BETARAND

template < typename RealType >
String BetaRand<RealType>::Name() const
{
    return "Beta(" + this->toStringWithPrecision(this->GetAlpha()) + ", "
                   + this->toStringWithPrecision(this->GetBeta()) + ", "
                   + this->toStringWithPrecision(this->MinValue()) + ", "
            + this->toStringWithPrecision(this->MaxValue()) + ")";
}

template < typename RealType >
DoublePair BetaRand<RealType>::SufficientStatistic(RealType x) const
{
    double y = (x - this->a) * this->bmaInv;
    return {std::log(y), std::log1p(-y)};
}

template < typename RealType >
DoublePair BetaRand<RealType>::SourceParameters() const
{
    return {this->alpha, this->beta};
}

template < typename RealType >
DoublePair BetaRand<RealType>::SourceToNatural(DoublePair sourceParameters) const
{
    return {sourceParameters.first - 1, sourceParameters.second - 1};
}

template < typename RealType >
double BetaRand<RealType>::LogNormalizer(DoublePair theta) const
{
    return this->logbma + RandMath::logBeta(theta.first + 1, theta.second + 1);
}

template < typename RealType >
DoublePair BetaRand<RealType>::LogNormalizerGradient(DoublePair theta) const
{
    double psi1 = RandMath::digamma(theta.first + 1);
    double psi2 = RandMath::digamma(theta.second + 1);
    double psisum = RandMath::digamma(theta.first + theta.second + 2);
    return {psi1 - psisum, psi2 - psisum};
}

template < typename RealType >
double BetaRand<RealType>::CarrierMeasure(RealType) const
{
    return 0;
}

template < typename RealType >
long double BetaRand<RealType>::GetSampleLogMeanNorm(const std::vector<RealType> &sample) const
{
    long double lnG = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG += std::log(x);
    }
    return lnG / sample.size();
}

template < typename RealType >
long double BetaRand<RealType>::GetSampleLog1pMeanNorm(const std::vector<RealType> &sample) const
{
    long double lnG1p = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG1p += std::log1pl(x);
    }
    return lnG1p / sample.size();
}

template < typename RealType >
long double BetaRand<RealType>::GetSampleLog1mMeanNorm(const std::vector<RealType> &sample) const
{
    long double lnG1m = 0;
    for (RealType var : sample) {
        RealType x = (var - this->a) * this->bmaInv;
        lnG1m += std::log1pl(-x);
    }
    return lnG1m / sample.size();
}

template < typename RealType >
void BetaRand<RealType>::FitAlpha(long double lnG, long double meanNorm)
{
    if (meanNorm <= 0 || meanNorm >= 1)
        throw std::invalid_argument(this->fitErrorDescription(this->NOT_APPLICABLE, "Normalized mean of the sample should be in interval of (0, 1)"));
    if (this->beta == 1.0) {
        /// for β = 1 we have explicit expression for estimator
        SetShapes(-1.0 / lnG, this->beta);
    }
    else {
        /// get initial value for shape by method of moments
        RealType shape = meanNorm;
        shape /= 1.0 - meanNorm;
        shape *= this->beta;
        /// run root-finding procedure
        if (!RandMath::findRootNewtonFirstOrder<RealType>([this, lnG] (double x)
        {
            double first = RandMath::digamma(x) - RandMath::digamma(x + this->beta) - lnG;
            double second = RandMath::trigamma(x);
            return DoublePair(first, second);
        }, shape))
            throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure."));
        SetShapes(shape, this->beta);
    }
}

template < typename RealType >
void BetaRand<RealType>::FitAlpha(const std::vector<RealType> &sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b)));

    long double lnG = this->GetSampleLogMeanNorm(sample);
    if (!std::isfinite(lnG))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->ALPHA_ZERO));
    long double mean = this->GetSampleMean(sample);
    mean -= this->a;
    mean *= this->bmaInv;
    FitAlpha(lnG, mean);
}

template < typename RealType >
void BetaRand<RealType>::FitBeta(long double lnG1m, long double meanNorm)
{
    if (meanNorm <= 0 || meanNorm >= 1)
        throw std::invalid_argument(this->fitErrorDescription(this->NOT_APPLICABLE, "Normalized mean of the sample should be in interval of (0, 1)"));
    if (this->alpha == 1.0) {
        /// for α = 1 we have explicit expression for estimator
        SetShapes(this->alpha, -1.0 / lnG1m);
    }
    else {
        /// get initial value for shape by method of moments
        RealType shape = this->alpha / meanNorm - this->alpha;
        /// run root-finding procedure
        if (!RandMath::findRootNewtonFirstOrder<RealType>([this, lnG1m] (double x)
        {
            double first = RandMath::digamma(x) - RandMath::digamma(x + this->alpha) - lnG1m;
            double second = RandMath::trigamma(x);
            return DoublePair(first, second);
        }, shape))
            throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure."));
        this->SetShapes(this->alpha, shape);
    }
}

template < typename RealType >
void BetaRand<RealType>::FitBeta(const std::vector<RealType> &sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b)));

    long double lnG1m = this->GetSampleLog1mMeanNorm(sample);
    if (!std::isfinite(lnG1m))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->BETA_ZERO));
    long double mean = this->GetSampleMean(sample);
    mean -= this->a;
    mean *= this->bmaInv;
    FitBeta(lnG1m, mean);
}

template < typename RealType >
void BetaRand<RealType>::FitShapes(long double lnG, long double lnG1m, long double mean, long double variance)
{
    /// get initial values for shapes by method of moments
    double scaledMean = (mean - this->a) * this->bmaInv;
    double scaledVar = variance * this->bmaInv * this->bmaInv;
    double temp = scaledMean * (1.0 - scaledMean) / scaledVar - 1.0;
    double shape1 = 0.001, shape2 = shape1;
    if (temp > 0) {
        shape1 = scaledMean * temp;
        shape2 = (1.0 - scaledMean) * temp;
    }
    DoublePair shapes = std::make_pair(shape1, shape2);

    /// run root-finding procedure
    if (!RandMath::findRootNewtonFirstOrder2d([lnG, lnG1m] (DoublePair x)
    {
        double digammaAlphapBeta = RandMath::digamma(x.first + x.second);
        double digammaAlpha = RandMath::digamma(x.first);
        double digammaBeta = RandMath::digamma(x.second);
        double first = lnG + digammaAlphapBeta - digammaAlpha;
        double second = lnG1m + digammaAlphapBeta - digammaBeta;
        return DoublePair(first, second);
    },
    [] (DoublePair x)
    {
        double trigammaAlphapBeta = RandMath::trigamma(x.first + x.second);
        double trigammaAlpha = RandMath::trigamma(x.first);
        double trigammaBeta = RandMath::trigamma(x.second);
        DoublePair first = std::make_pair(trigammaAlphapBeta - trigammaAlpha, trigammaAlphapBeta);
        DoublePair second = std::make_pair(trigammaAlphapBeta, trigammaAlphapBeta - trigammaBeta);
        return std::make_tuple(first, second);
    },
    shapes))
        throw std::runtime_error(this->fitErrorDescription(this->UNDEFINED_ERROR, "Error in root-finding procedure."));
    SetShapes(shapes.first, shapes.second);
}

template < typename RealType >
void BetaRand<RealType>::FitShapes(const std::vector<RealType> &sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b)));

    long double lnG = this->GetSampleLogMeanNorm(sample);
    if (!std::isfinite(lnG))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->ALPHA_ZERO));
    long double lnG1m = this->GetSampleLog1mMeanNorm(sample);
    if (!std::isfinite(lnG1m))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->BETA_ZERO));

    /// get initial values for shapes by method of moments
    DoublePair stats = this->GetSampleMeanAndVariance(sample);
    FitShapes(lnG, lnG1m, stats.first, stats.second);
}

template class BetaRand<float>;
template class BetaRand<double>;
template class BetaRand<long double>;

// ARCSINERAND

template < typename RealType >
String ArcsineRand<RealType>::Name() const
{
    return "Arcsine(" + this->toStringWithPrecision(GetShape()) + ", "
                      + this->toStringWithPrecision(this->MinValue()) + ", "
                      + this->toStringWithPrecision(this->MaxValue()) + ")";
}

template < typename RealType >
void ArcsineRand<RealType>::SetShape(double shape)
{
    BetaDistribution<RealType>::SetShapes(1.0 - shape, shape);
}

template < typename RealType >
void ArcsineRand<RealType>::FitShape(long double lnG, long double lnG1m)
{
    double shape = M_PI / (lnG1m - lnG);
    if (!std::isfinite(shape))
        SetShape(0.5);
    shape = -M_1_PI * RandMath::atan(shape);
    SetShape(shape > 0 ? shape : shape + 1);
}

template < typename RealType >
void ArcsineRand<RealType>::FitShape(const std::vector<RealType> &sample)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b)));

    int n = sample.size();
    long double lnG = 0, lnG1m = 0;
    for (double var : sample) {
        double x = (var - this->a) * this->bmaInv;
        lnG += std::log(x);
        lnG1m += std::log1pl(-x);
    }
    if (!std::isfinite(lnG))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->ALPHA_ZERO));
    if (!std::isfinite(lnG1m))
        throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, this->BETA_ZERO));
    lnG /= n;
    lnG1m /= n;
    FitShape(lnG, lnG1m);
}

template class ArcsineRand<float>;
template class ArcsineRand<double>;
template class ArcsineRand<long double>;

// BALDINGNICHOLSRAND

template < typename RealType >
BaldingNicholsRand<RealType>::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    SetFixatingIndexAndFrequency(fixatingIndex, frequency);
}

template < typename RealType >
String BaldingNicholsRand<RealType>::Name() const
{
    return "Balding-Nichols(" + this->toStringWithPrecision(GetFixatingIndex()) + ", " + this->toStringWithPrecision(GetFrequency()) + ")";
}

template < typename RealType >
void BaldingNicholsRand<RealType>::SetFixatingIndexAndFrequency(double fixatingIndex, double frequency)
{
    F = fixatingIndex;
    if (F <= 0 || F >= 1)
        F = 0.5;

    p = frequency;
    if (p <= 0 || p >= 1)
        p = 0.5;

    double frac = (1.0 - F) / F, fracP = frac * p;
    BetaDistribution<RealType>::SetShapes(fracP, frac - fracP);
}

template class BaldingNicholsRand<float>;
template class BaldingNicholsRand<double>;
template class BaldingNicholsRand<long double>;
