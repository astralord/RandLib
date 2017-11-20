#include "MarchenkoPasturRand.h"
#include "UniformRand.h"

template < typename RealType >
MarchenkoPasturRand<RealType>::MarchenkoPasturRand(double ratio, double scale)
{
    SetParameters(ratio, scale);
}

template < typename RealType >
String MarchenkoPasturRand<RealType>::Name() const
{
    return "Marchenko-Pastur(" + this->toStringWithPrecision(GetRatio()) + ", "
            + this->toStringWithPrecision(GetScale()) + ")";
}

template < typename RealType >
void MarchenkoPasturRand<RealType>::SetParameters(double ratio, double scale)
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
    logLambda = std::log(lambda);

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

template < typename RealType >
double MarchenkoPasturRand<RealType>::f(const RealType &x) const
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

template < typename RealType >
double MarchenkoPasturRand<RealType>::logf(const RealType &x) const
{
    if (x == 0.0)
        return (lambda < 1) ? -INFINITY : INFINITY;
    double xSt = x / sigmaSq;
    if (xSt < a || xSt > b)
        return -INFINITY;
    double y = 0.5 * std::log((b - xSt) * (xSt - a));
    y -= M_LN2 + M_LNPI + logLambda + std::log(x);
    return y;
}

template < typename RealType >
double MarchenkoPasturRand<RealType>::ccdfForLargeRatio(const RealType &x) const
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

template < typename RealType >
double MarchenkoPasturRand<RealType>::cdfForSmallRatio(const RealType &x) const
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

template < typename RealType >
double MarchenkoPasturRand<RealType>::F(const RealType &x) const
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

template < typename RealType >
double MarchenkoPasturRand<RealType>::S(const RealType &x) const
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

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::variateForTinyRatio() const
{
    size_t iter = 0;
    do {
        double X = BetaRV.Variate();
        double U = UniformRand::StandardVariate(this->localRandGenerator);
        if (U < M / X)
            return X;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    return NAN; /// fail due to some error
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::variateForSmallRatio() const
{
    size_t iter = 0;
    do {
        double X = BetaRV.Variate();
        double U = UniformRand::StandardVariate(this->localRandGenerator);
        double ratio = M * (1.0 - a / X);
        if (U * U < ratio)
            return X;
    } while (++iter <= ProbabilityDistribution<RealType>::MAX_ITER_REJECTION);
    return NAN; /// fail due to some error
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::variateForLargeRatio() const
{
    return (UniformRand::StandardVariate(this->localRandGenerator) > 1.0 / lambda) ? 0.0 : variateForSmallRatio();
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::variateForHugeRatio() const
{
    return (UniformRand::StandardVariate(this->localRandGenerator) > 1.0 / lambda) ? 0.0 : variateForTinyRatio();
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::Variate() const
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

template < typename RealType >
void MarchenkoPasturRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    switch (getIdOfUsedGenerator()) {
    case TINY_RATIO:
        for (RealType & var : outputData)
            var = sigmaSq * variateForTinyRatio();
        break;
    case SMALL_RATIO:
        for (RealType & var : outputData)
            var = sigmaSq * variateForSmallRatio();
        break;
    case LARGE_RATIO:
        for (RealType & var : outputData)
            var = sigmaSq * variateForLargeRatio();
        break;
    case HUGE_RATIO:
        for (RealType & var : outputData)
            var = sigmaSq * variateForHugeRatio();
        break;
    default:
        return;
    }
}

template < typename RealType >
void MarchenkoPasturRand<RealType>::Reseed(unsigned long seed) const
{
    this->localRandGenerator.Reseed(seed);
    BetaRV.Reseed(seed + 1);
}

template < typename RealType >
long double MarchenkoPasturRand<RealType>::Moment(int n) const
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
        long double sum = 0.0;
        long double temp = RandMath::lfact(n) + RandMath::lfact(n - 1);
        long double nlogSigmaSq = n * std::log(sigmaSq);
        for (int k = 0; k != n; ++k) {
            long double term = temp;
            term -= 2 * RandMath::lfact(k);
            term -= RandMath::lfact(n - k);
            term -= RandMath::lfact(n - k - 1);
            term += k * logLambda;
            term += nlogSigmaSq;
            term = std::exp(term) / (k + 1);
            sum += term;
        }
        return sum;
        }
    }
}

template < typename RealType >
long double MarchenkoPasturRand<RealType>::Mean() const
{
    return sigmaSq;
}

template < typename RealType >
long double MarchenkoPasturRand<RealType>::Variance() const
{
    return sigmaSq * sigmaSq * lambda;
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::Mode() const
{
    if (lambda > 1)
        return 0.0;
    RealType mode = lambda - 1.0;
    mode *= mode;
    mode /= lambda + 1.0;
    return sigmaSq * mode;
}

template < typename RealType >
long double MarchenkoPasturRand<RealType>::Skewness() const
{
    long double mu = Mean();
    long double var = Variance();
    long double skewness = Moment(3);
    skewness -= std::pow(mu, 3);
    skewness -= 3 * mu * var;
    skewness /= std::pow(var, 1.5);
    return skewness;
}

template < typename RealType >
long double MarchenkoPasturRand<RealType>::ExcessKurtosis() const
{
    long double mu = Mean();
    long double var = Variance();
    long double moment3 = Moment(3);
    long double kurtosis = Moment(4);
    long double muSq = mu * mu;
    kurtosis -= 4 * mu * moment3;
    kurtosis += 3 * muSq * muSq;
    kurtosis += 6 * muSq * var;
    kurtosis /= var * var;
    return kurtosis - 3.0;
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::quantileImpl(double p) const
{
    return (p < 1.0 - 1.0 / lambda) ? 0.0 : ContinuousDistribution<RealType>::quantileImpl(p);
}

template < typename RealType >
RealType MarchenkoPasturRand<RealType>::quantileImpl1m(double p) const
{
    return (p > 1.0 / lambda) ? 0.0 : ContinuousDistribution<RealType>::quantileImpl1m(p);
}

template < typename RealType >
std::complex<double> MarchenkoPasturRand<RealType>::CFImpl(double t) const
{
    if (lambda < 1)
        return ContinuousDistribution<RealType>::CFImpl(t);
    /// otherwise we have singularity at point 0
    if (lambda == 1) {
        /// we split integrand for real part on (cos(tx)-1)f(x) and f(x)
        double re = this->ExpectedValue([this, t] (double x)
        {
            return std::cos(t * x) - 1.0;
        }, 0, 4 * sigmaSq);
        double im = this->ExpectedValue([this, t] (double x)
        {
            return std::sin(t * x);
        }, 0, 4 * sigmaSq);
        return std::complex<double>(1.0 + re, im);
    }
    /// for λ > 1 we split integral on 2 parts: at point 0 and the rest
    double re = this->ExpectedValue([this, t] (double x)
    {
        return std::cos(t * x);
    }, sigmaSq * a, sigmaSq * b);
    double im = this->ExpectedValue([this, t] (double x)
    {
        return std::sin(t * x);
    }, sigmaSq * a, sigmaSq * b);
    return std::complex<double>(1.0 - 1.0 / lambda + re, im);
}

