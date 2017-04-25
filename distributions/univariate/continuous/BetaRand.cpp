#include "BetaRand.h"
#include "../discrete/BernoulliRand.h"
#include "UniformRand.h"

BetaRand::BetaRand(double shape1, double shape2)
{
    SetShapes(shape1, shape2);
    SetSupport(0, 1);
}

BetaRand::BetaRand(double shape1, double shape2, double minValue, double maxValue)
{
    SetShapes(shape1, shape2);
    SetSupport(minValue, maxValue);
}

std::string BetaRand::Name() const
{
    return "Beta(" + toStringWithPrecision(GetAlpha()) + ", "
                   + toStringWithPrecision(GetBeta()) + ", "
                   + toStringWithPrecision(MinValue()) + ", "
                   + toStringWithPrecision(MaxValue()) + ")";
}

BetaRand::GENERATOR_ID BetaRand::getIdOfUsedGenerator() const
{
    if (alpha < 1 && beta < 1 && alpha + beta > 1)
        return ATKINSON_WHITTAKER;

    if (RandMath::areClose(alpha, beta)) {
        /// Generators for symmetric density
        if (RandMath::areClose(alpha, 1.0))
            return UNIFORM; /// Standard uniform variate
        else if (RandMath::areClose(alpha, 0.5))
            return ARCSINE; /// Arcsine method
        else if (RandMath::areClose(alpha, 1.5))
            return REJECTION_UNIFORM; /// Rejection method from uniform distribution
        else if (alpha > 1)
            /// Rejection from uniform or normal distribution
            return (alpha < 2) ? REJECTION_UNIFORM_EXTENDED : REJECTION_NORMAL;
    }
    if (std::min(alpha, beta) > 0.5 && std::max(alpha, beta) > 1)
        return CHENG;
    return (alpha + beta < 2) ? JOHNK : GAMMA_RATIO;
}

void BetaRand::setCoefficientsForGenerator()
{
    GENERATOR_ID id = getIdOfUsedGenerator();
    if (id == REJECTION_NORMAL) {
        double alpham1 = alpha - 1;
        s = alpham1 * std::log1p(0.5 / alpham1) - 0.5;
        t = 1.0 / std::sqrt(8 * alpha - 4);
    }
    else if (id == CHENG) {
        s = alpha + beta;
        t = std::min(alpha, beta);
        if (t > 1)
            t = std::sqrt((2 * alpha * beta - s) / (s - 2));
        u = alpha + t;
    }
    else if (id == ATKINSON_WHITTAKER) {
        t = std::sqrt(alpha * (1 - alpha));
        t /= (t + std::sqrt(beta * (1 - beta)));
        s = beta * t;
        s /= (s + alpha * (1 - t));
    }
}

void BetaRand::SetShapes(double shape1, double shape2)
{
    GammaRV1.SetParameters(shape1, 1);
    GammaRV2.SetParameters(shape2, 1);
    alpha = GammaRV1.GetShape();
    beta = GammaRV2.GetShape();
    /// we use log(Gamma(x)) in order to avoid too big numbers
    logBetaFun = -std::lgamma(alpha + beta) + GammaRV1.GetLogGammaFunction() + GammaRV2.GetLogGammaFunction();
    betaFun = std::exp(logBetaFun);
    setCoefficientsForGenerator();
}

void BetaRand::SetSupport(double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;
    /// Sanity check
    if (a >= b)
        b = a + 1.0;
    bma = b - a;
    bmaInv = 1.0 / bma;
    logBma = std::log(bma);
}

double BetaRand::f(const double & x) const
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

double BetaRand::logf(const double & x) const
{
    if (x < a || x > b)
        return -INFINITY;
    if (x == a) {
        if (alpha == 1)
            return std::log(beta / bma);
        return (alpha > 1) ? -INFINITY : INFINITY;
    }
    if (x == b) {
        if (beta == 1)
            return std::log(alpha / bma);
        return (beta > 1) ? -INFINITY : INFINITY;
    }
    /// Standardize
    double xSt = (x - a) / bma;
    double y = (alpha - 1) * std::log(xSt);
    y += (beta - 1) * std::log1p(-xSt);
    return y - logBetaFun - logBma;
}

double BetaRand::F(const double & x) const
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
    return RandMath::ibeta(xSt, alpha, beta, logBetaFun, std::log(xSt), std::log1p(-xSt));
}

double BetaRand::S(const double & x) const
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
    return RandMath::ibeta(1.0 - xSt, beta, alpha, logBetaFun, std::log1p(-xSt), std::log(xSt));
}

double BetaRand::variateArcsine() const
{
    double U = UniformRand::Variate(-M_PI, M_PI);
    double X = std::sin(U);
    return X * X;
}

double BetaRand::variateRejectionUniform() const
{
    /// Rejection method from uniform distribution
    /// boosted for specific value of shapes α = β = 1.5
    int iter = 0;
    do {
        double U = UniformRand::StandardVariate();
        double V = UniformRand::StandardVariate();
        if (0.25 * V * V <= U - U * U)
            return U;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail
}

double BetaRand::variateRejectionUniformExtended() const
{
    /// Rejection method from uniform distribution
    /// boosted by using exponential distribution
    int iter = 0;
    static constexpr double M_LN4 = M_LN2 + M_LN2;
    do {
        double U = UniformRand::StandardVariate();
        double W = ExponentialRand::StandardVariate();
        double edge = M_LN4 + std::log(U - U * U);
        if (W >= (1.0 - alpha) * edge)
            return U;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail
}

double BetaRand::variateCheng() const
{
    double R, T, Y;
    do {
        double U = UniformRand::StandardVariate();
        double V = UniformRand::StandardVariate();
        double X = std::log(U / (1 - U)) / t;
        Y = alpha * std::exp(X);
        R = 1.0 / (beta + Y);
        T = 4 * U * U * V;
        T = std::log(T);
        T -= u * X;
        T -= s * std::log(s * R);
    } while (T > 0);
    return Y * R;
}

double BetaRand::variateAtkinsonWhittaker() const
{
    int iter = 0;
    do {
        double U = UniformRand::StandardVariate();
        double W = ExponentialRand::StandardVariate();
        if (U <= s) {
            double X = t * std::pow(U / s, 1.0 / alpha);
            if (W >= (1.0 - beta) * std::log((1.0 - X) / (1.0 - t)))
                return X;
        }
        else {
            double X = 1.0 - (1.0 - t) * std::pow((1.0 - U) / (1.0 - s), 1.0 / beta);
            if (W >= (1.0 - alpha) * std::log(X / t))
                return X;
        }
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail
}

double BetaRand::variateGammaRatio() const
{
    double Y = GammaRV1.Variate();
    double Z = GammaRV2.Variate();
    return Y / (Y + Z);
}

double BetaRand::variateRejectionNormal() const
{
    int iter = 0;
    double N = 0, Z = 0;
    double alpham1 = alpha - 1;
    double alpha2m1 = alpha + alpham1;
    do {
        do {
            N = NormalRand::StandardVariate();
            Z = N * N;
        } while (Z >= alpha2m1);

        double W = ExponentialRand::StandardVariate() + s;
        double aux = 0.5 - alpham1 / (alpha2m1 - Z);
        aux *= Z;
        if (W + aux >= 0)
            return 0.5 + N * t;
        aux = std::log1p(-Z / alpha2m1);
        aux *= alpham1;
        aux += W + 0.5 * Z;
        if (aux >= 0)
            return 0.5 + N * t;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail
}

double BetaRand::variateJohnk() const
{
    double X = 0, Z = 0;
    do {
        double U = UniformRand::StandardVariate();
        double V = UniformRand::StandardVariate();
        X = std::pow(U, 1.0 / alpha);
        Z = X + std::pow(V, 1.0 / beta);
    } while (Z > 1);
    return X / Z;
}

double BetaRand::Variate() const
{
    double var = 0;
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM:
        var = UniformRand::StandardVariate();
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

void BetaRand::Sample(std::vector<double> &outputData) const
{
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM: {
        for (double &var : outputData)
            var = UniformRand::StandardVariate();
        }
        break;
    case ARCSINE: {
        for (double &var : outputData)
            var = variateArcsine();
        }
        break;
    case CHENG: {
        for (double &var : outputData)
            var = variateCheng();
        }
        break;
    case REJECTION_UNIFORM: {
        for (double &var : outputData)
            var = variateRejectionUniform();
        }
        break;
    case REJECTION_UNIFORM_EXTENDED: {
        for (double &var : outputData)
            var = variateRejectionUniformExtended();
        }
        break;
    case REJECTION_NORMAL: {
        for (double &var : outputData)
            var = variateRejectionNormal();
        }
        break;
    case JOHNK: {
        for (double &var : outputData)
            var = variateJohnk();
        }
        break;
    case ATKINSON_WHITTAKER: {
        for (double &var : outputData)
            var = variateAtkinsonWhittaker();
        }
        break;
    case GAMMA_RATIO:
    default: {
        GammaRV1.Sample(outputData);
        for (double &var : outputData)
            var /= (var + GammaRV2.Variate());
        }
        break;
    }

    /// Shift and scale
    for (double &var : outputData)
        var = a + bma * var;
}

double BetaRand::Mean() const
{
    double mean = alpha / (alpha + beta);
    return a + bma * mean;
}

double BetaRand::Variance() const
{
    double var = alpha + beta;
    var *= var * (var + 1);
    var = alpha * beta / var;
    return bma * bma * var;
}

double BetaRand::Median() const
{
    return (alpha == beta) ? 0.5 * (b + a) : quantileImpl(0.5);
}

double BetaRand::Mode() const
{
    double mode;
    if (alpha > 1)
        mode = (beta > 1) ? (alpha - 1) / (alpha + beta - 2) : 1.0;
    else
        mode = (beta > 1) ? 0.0 : (alpha > beta);
    return a + bma * mode;
}

double BetaRand::Skewness() const
{
    double skewness = (alpha + beta + 1) / (alpha * beta);
    skewness = std::sqrt(skewness);
    skewness *= beta - alpha;
    skewness /= alpha + beta + 2;
    return 2 * skewness;
}

double BetaRand::ExcessKurtosis() const
{
    double sum = alpha + beta;
    double kurtosis = alpha - beta;
    kurtosis *= kurtosis;
    kurtosis *= (sum + 1);
    kurtosis /= (alpha * beta * (sum + 2));
    --kurtosis;
    kurtosis /= (sum + 3);
    return 6 * kurtosis;
}

double BetaRand::quantileImpl(double p) const
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
        return a - bma * std::expm1(std::log1p(-p) / beta);
    if (beta == 1.0)
        return a + bma * std::pow(p, 1.0 / alpha);
    return ContinuousDistribution::quantileImpl(p);
}

double BetaRand::quantileImpl1m(double p) const
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
        return a - bma * std::expm1(std::log(p) / beta);
    if (beta == 1.0)
        return a + bma * std::exp(std::log1p(-p) / alpha);
    return ContinuousDistribution::quantileImpl1m(p);
}

std::complex<double> BetaRand::CFImpl(double t) const
{
    /// if we don't have singularity points, we can use direct integration
    if (alpha >= 1 && beta >= 1)
        return UnivariateProbabilityDistribution::CFImpl(t);

    double z = bma * t;
    double sinZ = std::sin(z);
    double cosZm1 = std::cos(z) - 1.0;

    double re = RandMath::integral([this, z, cosZm1](double x) {
        if (x >= 1)
            return 0.0;
        if (x <= 0)
            return -cosZm1;
        double f = std::cos(z * x) - 1;
        f *= std::pow(x, alpha - 1);
        f -= cosZm1;
        return std::pow(1.0 - x, beta - 1) * f;
    }, 0, 1);
    re += betaFun;
    re += cosZm1 / beta;

    double im = RandMath::integral([this, z, sinZ](double x) {
        if (x >= 1)
            return 0.0;
        if (x <= 0)
            return -sinZ;
        double f = std::sin(z * x);
        f *= std::pow(x, alpha - 1);
        f -= sinZ;
        return std::pow(1.0 - x, beta - 1) * f;
    }, 0, 1);
    im += sinZ / beta;

    std::complex<double> y(re, im);
    double cosTA = std::cos(t * a), sinTA = std::sin(t * a);
    return y * std::complex<double>(cosTA, sinTA) / betaFun;
}

void BetaRand::FitAlphaMM(const std::vector<double> &sample)
{
    if (!allElementsAreNotLessThen(a, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    if (!allElementsAreNotBiggerThen(b, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    double mean = sampleMean(sample);
    double shape = mean - a;
    shape /= b - mean;
    shape *= beta;
    SetShapes(shape, beta);
}

void BetaRand::FitBetaMM(const std::vector<double> &sample)
{
    if (!allElementsAreNotLessThen(a, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, LOWER_LIMIT_VIOLATION + toStringWithPrecision(a)));
    if (!allElementsAreNotBiggerThen(b, sample))
        throw std::invalid_argument(fitError(WRONG_SAMPLE, UPPER_LIMIT_VIOLATION + toStringWithPrecision(b)));
    double mean = sampleMean(sample);
    double shape = b - mean;
    shape /= mean - a;
    shape *= alpha;
    SetShapes(alpha, shape);
}

ArcsineRand::ArcsineRand(double shape, double minValue, double maxValue) :
    BetaRand(1.0 - shape, shape, minValue, maxValue)
{
}

std::string ArcsineRand::Name() const
{
    return "Arcsine(" + toStringWithPrecision(GetShape()) + ", "
                      + toStringWithPrecision(MinValue()) + ", "
                      + toStringWithPrecision(MaxValue()) + ")";
}

void ArcsineRand::SetShape(double shape)
{
    BetaRand::SetShapes(1.0 - shape, shape);
}


BaldingNicholsRand::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    SetFixatingIndexAndFrequency(fixatingIndex, frequency);
}

std::string BaldingNicholsRand::Name() const
{
    return "Balding-Nichols(" + toStringWithPrecision(GetFixatingIndex()) + ", " + toStringWithPrecision(GetFrequency()) + ")";
}

void BaldingNicholsRand::SetFixatingIndexAndFrequency(double fixatingIndex, double frequency)
{
    F = fixatingIndex;
    if (F <= 0 || F >= 1)
        F = 0.5;

    p = frequency;
    if (p <= 0 || p >= 1)
        p = 0.5;

    double frac = (1.0 - F) / F, fracP = frac * p;
    BetaRand::SetShapes(fracP, frac - fracP);
}
