#include "MarchenkoPasturRand.h"
#include "UniformRand.h"

MarchenkoPasturRand::MarchenkoPasturRand(double ratio, double scale)
{
    SetParameters(ratio, scale);
}

String MarchenkoPasturRand::Name() const
{
    return "Marchenko-Pastur(" + toStringWithPrecision(GetRatio()) + ", "
            + toStringWithPrecision(GetScale()) + ")";
}

void MarchenkoPasturRand::SetParameters(double ratio, double scale)
{
    if (ratio <= 0.0)
        throw std::invalid_argument("Marchenko-Pastur distribution: ratio parameter should be positive");
    if (scale <= 0.0)
        throw std::invalid_argument("Marchenko-Pastur distribution: scale should be positive");
    lambda = ratio;
    double sqrtLambda = std::sqrt(lambda);
    a = 1.0 - sqrtLambda;
    a *= a;
    b = 1.0 + sqrtLambda;
    b *= b;

    sigmaSq = scale;

    GENERATOR_ID genId = getIdOfUsedGenerator();
    if (genId == TINY_RATIO || genId == HUGE_RATIO) {
        BetaRV.SetShapes(1.5, 1.5);
        BetaRV.SetSupport(a, b);
        M = std::min(a / lambda, a);
    }
    else {
        BetaRV.SetShapes(0.5, 1.5);
        BetaRV.SetSupport(0, b);
        M = 0.25 * b / sqrtLambda;
        if (lambda > 1)
            M /= lambda * lambda;
    }
}

double MarchenkoPasturRand::f(const double &x) const
{
    if (x == 0.0)
        return (lambda < 1) ? 0.0 : INFINITY;
    double xSt = x / sigmaSq;
    if (xSt < a || xSt > b)
        return 0.0;
    double y = 0.5 * std::sqrt((b - xSt) * (xSt - a));
    y /= (lambda * xSt);
    y /= (M_PI * sigmaSq);
    return y;
}

double MarchenkoPasturRand::logf(const double &x) const
{
    if (x == 0.0)
        return (lambda < 1) ? -INFINITY : INFINITY;
    double xSt = x / sigmaSq;
    if (xSt < a || xSt > b)
        return -INFINITY;
    double y = 0.5 * std::log((b - xSt) * (xSt - a));
    y -= std::log(2 * M_PI * sigmaSq * lambda * xSt);
    return y;
}

double MarchenkoPasturRand::ccdfForLargeRatio(const double &x) const
{
    double y1 = 1.0 - x + lambda;
    double lambdam1 = lambda - 1.0;
    double lambdap1 = lambda + 1.0;
    double temp = std::sqrt(4 * lambda - y1 * y1);
    if (temp != 0.0) {
        y1 /= temp;
        y1 = lambdap1 * RandMath::atan(y1);
    }
    else
        y1 = RandMath::sign(y1) * M_PI_2;
    double y2 = x * lambdap1;
    y2 -= lambdam1 * lambdam1;
    y2 /= temp * lambdam1;
    y2 = lambdam1 * RandMath::atan(y2);
    double y = M_PI - temp + y1 + y2;
    y /= M_PI * lambda;
    return 0.5 * y;
}

double MarchenkoPasturRand::cdfForSmallRatio(const double &x) const
{
    double y1 = 1.0 - x + lambda;
    double temp = std::sqrt(4 * lambda - y1 * y1);
    double lambdam1 = lambda - 1.0;
    double lambdap1 = lambda + 1.0;
    if (temp != 0.0) {
        y1 /= temp;
        y1 = lambdap1 * RandMath::atan(y1);
    }
    else
        y1 = RandMath::sign(y1) * M_PI_2;
    double y2 = 0.0;
    if (lambdam1 != 0)
    {
        y2 = x * lambdap1;
        y2 -= lambdam1 * lambdam1;
        y2 /= -temp * lambdam1;
        y2 = lambdam1 * RandMath::atan(y2);
    }
    double y = M_PI * lambda + temp - y1 + y2;
    y /= M_PI * lambda;
    return 0.5 * y;
}

double MarchenkoPasturRand::F(const double &x) const
{
    double xSt = x / sigmaSq;
    if (xSt < 0.0)
        return 0.0;
    if (xSt >= b)
        return 1.0;
    if (lambda > 1.0)
        return (xSt > a) ? 1.0 - ccdfForLargeRatio(xSt) : 1.0 - 1.0 / lambda;
    return (xSt > a) ? cdfForSmallRatio(xSt) : 0.0;
}

double MarchenkoPasturRand::S(const double &x) const
{
    double xSt = x / sigmaSq;
    if (xSt < 0.0)
        return 1.0;
    if (xSt >= b)
        return 0.0;
    if (lambda > 1.0)
        return (xSt > a) ? ccdfForLargeRatio(xSt) : 1.0 / lambda;
    return (xSt > a) ? 1.0 - cdfForSmallRatio(xSt) : 1.0;
}

MarchenkoPasturRand::GENERATOR_ID MarchenkoPasturRand::getIdOfUsedGenerator() const
{
    if (lambda < 0.3)
        return TINY_RATIO;
    if (lambda <= 1.0)
        return SMALL_RATIO;
    return (lambda > 3.3) ? HUGE_RATIO : LARGE_RATIO;
}

double MarchenkoPasturRand::variateForTinyRatio() const
{
    int iter = 0;
    do {
        double X = BetaRV.Variate();
        double U = UniformRand::StandardVariate();
        if (U < M / X)
            return X;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail due to some error
}

double MarchenkoPasturRand::variateForSmallRatio() const
{
    int iter = 0;
    do {
        double X = BetaRV.Variate();
        double U = UniformRand::StandardVariate();
        double ratio = M * (1.0 - a / X);
        if (U * U < ratio)
            return X;
    } while (++iter <= MAX_ITER_REJECTION);
    return NAN; /// fail due to some error
}

double MarchenkoPasturRand::variateForLargeRatio() const
{
    return (UniformRand::StandardVariate() > 1.0 / lambda) ? 0.0 : variateForSmallRatio();
}

double MarchenkoPasturRand::variateForHugeRatio() const
{
    return (UniformRand::StandardVariate() > 1.0 / lambda) ? 0.0 : variateForTinyRatio();
}

double MarchenkoPasturRand::Variate() const
{
    switch (getIdOfUsedGenerator()) {
    case TINY_RATIO:
        return sigmaSq * variateForTinyRatio();
    case SMALL_RATIO:
        return sigmaSq * variateForSmallRatio();
    case LARGE_RATIO:
        return sigmaSq * variateForLargeRatio();
    case HUGE_RATIO:
        return sigmaSq * variateForHugeRatio();
    default:
        return NAN;
    }
}

void MarchenkoPasturRand::Sample(std::vector<double> &outputData) const
{
    switch (getIdOfUsedGenerator()) {
    case TINY_RATIO:
        for (double & var : outputData)
            var = sigmaSq * variateForTinyRatio();
        break;
    case SMALL_RATIO:
        for (double & var : outputData)
            var = sigmaSq * variateForSmallRatio();
        break;
    case LARGE_RATIO:
        for (double & var : outputData)
            var = sigmaSq * variateForLargeRatio();
        break;
    case HUGE_RATIO:
        for (double & var : outputData)
            var = sigmaSq * variateForHugeRatio();
        break;
    default:
        return;
    }
}

double MarchenkoPasturRand::Moment(int n) const
{
    if (n < 0)
        return NAN;
    switch (n) {
    case 0:
        return 1.0;
    case 1:
        return sigmaSq;
    case 2:
        return sigmaSq * sigmaSq * lambda;
    default: {
        double sum = 0.0;
        for (int k = 0; k != n; ++k) {
            double addon = RandMath::binom(n - 1, k);
            addon *= addon;
            addon /= n - k;
            addon /= k + 1;
            addon *= n * std::pow(lambda, k);
            sum += addon;
        }
        return std::pow(sigmaSq, n) * sum;
        }
    }
}

double MarchenkoPasturRand::Mean() const
{
    return sigmaSq;
}

double MarchenkoPasturRand::Variance() const
{
    return sigmaSq * sigmaSq * lambda;
}

double MarchenkoPasturRand::Mode() const
{
    if (lambda > 1)
        return 0.0;
    double mode = lambda - 1.0;
    mode *= mode;
    mode /= lambda + 1.0;
    return sigmaSq * mode;
}

double MarchenkoPasturRand::Skewness() const
{
    double mu = Mean();
    double var = Variance();
    double skewness = Moment(3);
    skewness -= std::pow(mu, 3);
    skewness -= 3 * mu * var;
    skewness /= std::pow(var, 1.5);
    return skewness;
}

double MarchenkoPasturRand::ExcessKurtosis() const
{
    double mu = Mean();
    double var = Variance();
    double moment3 = Moment(3);
    double kurtosis = Moment(4);
    double muSq = mu * mu;
    kurtosis -= 4 * mu * moment3;
    kurtosis += 3 * muSq * muSq;
    kurtosis += 6 * muSq * var;
    kurtosis /= var * var;
    return kurtosis - 3.0;
}

double MarchenkoPasturRand::quantileImpl(double p) const
{
    return (p < 1.0 - 1.0 / lambda) ? 0.0 : ContinuousDistribution::quantileImpl(p);
}

double MarchenkoPasturRand::quantileImpl1m(double p) const
{
    return (p > 1.0 / lambda) ? 0.0 : ContinuousDistribution::quantileImpl1m(p);
}

std::complex<double> MarchenkoPasturRand::CFImpl(double t) const
{
    if (lambda < 1)
        return ContinuousDistribution::CFImpl(t);
    /// otherwise we have singularity at point 0
    if (lambda == 1) {
        /// we split integrand for real part on (cos(tx)-1)f(x) and f(x)
        double re = ExpectedValue([this, t] (double x)
        {
            return std::cos(t * x) - 1.0;
        }, 0, 4 * sigmaSq);
        double im = ExpectedValue([this, t] (double x)
        {
            return std::sin(t * x);
        }, 0, 4 * sigmaSq);
        return std::complex<double>(1.0 + re, im);
    }
    /// for λ > 1 we split integral on 2 parts: at point 0 and the rest
    double re = ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    }, sigmaSq * a, sigmaSq * b);
    double im = ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, sigmaSq * a, sigmaSq * b);
    return std::complex<double>(1.0 - 1.0 / lambda + re, im);
}

