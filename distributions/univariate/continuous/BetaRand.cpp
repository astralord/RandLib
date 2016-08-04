#include "BetaRand.h"
#include "../discrete/BernoulliRand.h"
#include "UniformRand.h"

BetaRand::BetaRand(double shape1, double shape2, double minValue, double maxValue)
{
    setParameters(shape1, shape2, minValue, maxValue);
}

std::string BetaRand::name() const
{
    return "Beta(" + toStringWithPrecision(getAlpha()) + ", "
                   + toStringWithPrecision(getBeta()) + ", "
                   + toStringWithPrecision(getMin()) + ", "
                   + toStringWithPrecision(getMax()) + ")";
}

void BetaRand::setParameters(double shape1, double shape2, double minValue, double maxValue)
{
    a = minValue;
    b = maxValue;

    if (a >= b)
        b = a + 1.0;

    bma = b - a;

    GammaRV1.setParameters(shape1, 1);
    GammaRV2.setParameters(shape2, 1);

    alpha = GammaRV1.getShape();
    beta = GammaRV2.getShape();

    /// we use log(Gamma(x)) in order to avoid too big numbers
    betaFunInv = std::lgamma(alpha + beta) - GammaRV1.getLogGammaFunction() - GammaRV2.getLogGammaFunction();
    mLogBetaFun = betaFunInv - std::log(bma);
    betaFunInv = std::exp(betaFunInv);

    setCoefficientsForGenerator();
}

double BetaRand::f(double x) const
{
    if (x < a || x > b)
        return 0;

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

    /// Standardize
    x -= a;
    x /= bma;

    double y = 0.0;
    if (alpha == beta)
        y = (alpha - 1) * std::log(x - x * x);
    else {
        y = (alpha - 1) * std::log(x);
        y += (beta - 1) * std::log(1 - x);
    }
    return std::exp(mLogBetaFun + y);
}

double BetaRand::F(double x) const
{
    if (x <= a)
        return 0;
    if (x >= b)
        return 1;

    /// Standardize
    x -= a;
    x /= bma;

    if (alpha == beta && beta == 0.5)
        return M_2_PI * std::asin(std::sqrt(x));
    return betaFunInv * RandMath::incompleteBetaFun(x, alpha, beta);
}

double BetaRand::variateArcsine() const
{
    double U = UniformRand::variate(-M_PI, M_PI);
    double X = std::sin(U);
    return X * X;
}

double BetaRand::variateRejectionUniform() const
{
    /// Rejection method from uniform distribution
    /// boosted for specific value of shapes α = β = 1.5
    int iter = 0;
    do {
        double U = UniformRand::standardVariate();
        double V = UniformRand::standardVariate();
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
        double U = UniformRand::standardVariate();
        double W = ExponentialRand::standardVariate();
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
        double U = UniformRand::standardVariate();
        double V = UniformRand::standardVariate();
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
        double U = UniformRand::standardVariate();
        double W = ExponentialRand::standardVariate();
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
    double Y = GammaRV1.variate();
    double Z = GammaRV2.variate();
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
            N = NormalRand::standardVariate();
            Z = N * N;
        } while (Z >= alpha2m1);

        double W = ExponentialRand::standardVariate() + s;
        double aux = 0.5 - alpham1 / (alpha2m1 - Z);
        aux *= Z;
        if (W + aux >= 0)
            return 0.5 + N * t;
        aux = std::log(1.0 - Z / alpha2m1);
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
        double U = UniformRand::standardVariate();
        double V = UniformRand::standardVariate();
        X = std::pow(U, 1.0 / alpha);
        Z = X + std::pow(V, 1.0 / beta);
    } while (Z > 1);
    return X / Z;
}

double BetaRand::variate() const
{
    double var = 0;
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM:
        var = UniformRand::standardVariate();
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

void BetaRand::sample(std::vector<double> &outputData) const
{
    GENERATOR_ID id = getIdOfUsedGenerator();

    switch (id) {
    case UNIFORM: {
        for (double &var : outputData)
            var = UniformRand::standardVariate();
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
        GammaRV1.sample(outputData);
        for (double &var : outputData)
            var /= (var + GammaRV2.variate());
        }
        break;
    }

    /// Shift and scale
    for (double &var : outputData)
        var = a + bma * var;
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

std::complex<double> BetaRand::CF(double t) const
{
    if (t == 0)
        return 1;
    /// if we don't have singularity points, we can use direct integration
    if (alpha >= 1 && beta >= 1)
        return UnivariateProbabilityDistribution::CF(t);

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
    re += 1.0 / betaFunInv;
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
    y *= std::exp(std::complex<double>(0, t * a));
    return betaFunInv * y;
}

double BetaRand::Quantile(double p) const
{
    if (p < 0 || p > 1)
        return NAN;
    double root = p;
    if (alpha == beta && beta == 0.5)
    {
        double x = std::sin(0.5 * M_PI * p);
        return a + bma * x * x;
    }
    if (RandMath::findRoot([this, p] (double x)
    {
        return BetaRand::F(x) - p;
    },
    a, b, root))
        return root;
    return NAN;
}

double BetaRand::Median() const
{
    return (alpha == beta) ? 0.5 * (b + a) : Quantile(0.5);
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
    skewness *= (alpha - beta);
    skewness /= (alpha + beta + 2);
    return skewness + skewness;
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


ArcsineRand::ArcsineRand(double shape, double minValue, double maxValue)
{
    BetaRand::setParameters(1.0 - shape, shape, minValue, maxValue);
}

std::string ArcsineRand::name() const
{
    return "Arcsine(" + toStringWithPrecision(getShape()) + ", "
                      + toStringWithPrecision(getMin()) + ", "
                      + toStringWithPrecision(getMax()) + ")";
}

void ArcsineRand::setShape(double shape)
{
    BetaRand::setParameters(1.0 - shape, shape, a, b);
}


BaldingNicholsRand::BaldingNicholsRand(double fixatingIndex, double frequency)
{
    setFixatingIndexAndFrequency(fixatingIndex, frequency);
}

std::string BaldingNicholsRand::name() const
{
    return "Balding-Nichols(" + toStringWithPrecision(getFixatingIndex()) + ", " + toStringWithPrecision(getFrequency()) + ")";
}

void BaldingNicholsRand::setFixatingIndexAndFrequency(double fixatingIndex, double frequency)
{
    F = fixatingIndex;
    if (F <= 0 || F >= 1)
        F = 0.5;

    p = frequency;
    if (p <= 0 || p >= 1)
        p = 0.5;

    double frac = (1.0 - F) / F, fracP = frac * p;
    BetaRand::setParameters(fracP, frac - fracP);
}
