#include "UniformRand.h"
#include "../BasicRandGenerator.h"
#include <algorithm>
#include <limits>

template < typename RealType >
UniformRand<RealType>::UniformRand(double minValue, double maxValue) :
    BetaDistribution<RealType>(1, 1, minValue, maxValue)
{
}

template < typename RealType >
String UniformRand<RealType>::Name() const
{
    return "Uniform(" + this->toStringWithPrecision(MinValue()) + ", " + this->toStringWithPrecision(MaxValue()) + ")";
}

template < typename RealType >
double UniformRand<RealType>::f(const RealType & x) const
{
    return (x < this->a || x > this->b) ? 0.0 : this->bmaInv;
}

template < typename RealType >
double UniformRand<RealType>::logf(const RealType & x) const
{
    return (x < this->a || x > this->b) ? -INFINITY : -this->logbma;
}

template < typename RealType >
double UniformRand<RealType>::F(const RealType & x) const
{
    if (x < this->a)
        return 0.0;
    return (x > this->b) ? 1.0 : this->bmaInv * (x - this->a);
}

template < typename RealType >
double UniformRand<RealType>::S(const RealType & x) const
{
    if (x < this->a)
        return 1.0;
    return (x > this->b) ? 0.0 : this->bmaInv * (this->b - x);
}

template < typename RealType >
RealType UniformRand<RealType>::Variate() const
{
    return this->a + StandardVariate(this->localRandGenerator) * this->bma;
}

template < typename RealType >
RealType UniformRand<RealType>::StandardVariate(RandGenerator &randGenerator)
{
    /// (k + 1/2) / 2^b for b random bits, where b is one less than the mantissa width:
    /// the result lies in (0, 1) and its largest value 1 - 2^-(b+1) is representable,
    /// so rounding can never produce 1
    constexpr int bits = std::min(std::numeric_limits<RealType>::digits - 1, 63);
    RealType k = static_cast<RealType>(randGenerator.Bits(bits));
    return std::ldexp(k + static_cast<RealType>(0.5), -bits);
}

template < typename RealType >
RealType UniformRand<RealType>::StandardVariateClosed(RandGenerator &randGenerator)
{
    constexpr int bits = std::min(std::numeric_limits<RealType>::digits, 63);
    RealType k = static_cast<RealType>(randGenerator.Bits(bits));
    return k / (std::ldexp(static_cast<RealType>(1), bits) - 1);
}

template < typename RealType >
RealType UniformRand<RealType>::StandardVariateHalfClosed(RandGenerator &randGenerator)
{
    constexpr int bits = std::min(std::numeric_limits<RealType>::digits, 63);
    RealType k = static_cast<RealType>(randGenerator.Bits(bits));
    return std::ldexp(k, -bits);
}

template < typename RealType >
void UniformRand<RealType>::Sample(std::vector<RealType> &outputData) const
{
    for (RealType & var : outputData)
        var = this->Variate();
}

template < typename RealType >
long double UniformRand<RealType>::Mean() const
{
    return 0.5 * (this->b + this->a);
}

template < typename RealType >
long double UniformRand<RealType>::Variance() const
{
    return this->bma * this->bma / 12;
}

template < typename RealType >
RealType UniformRand<RealType>::Median() const
{
    return 0.5 * (this->b + this->a);
}

template < typename RealType >
RealType UniformRand<RealType>::Mode() const
{
    /// this can be any value in [a, b]
    return 0.5 * (this->b + this->a);
}

template < typename RealType >
long double UniformRand<RealType>::Skewness() const
{
    return 0.0;
}

template < typename RealType >
long double UniformRand<RealType>::ExcessKurtosis() const
{
    return -1.2;
}

template < typename RealType >
RealType UniformRand<RealType>::quantileImpl(double p) const
{
    return this->a + this->bma * p;
}

template < typename RealType >
RealType UniformRand<RealType>::quantileImpl1m(double p) const
{
    return this->b - this->bma * p;
}

template < typename RealType >
std::complex<double> UniformRand<RealType>::CFImpl(double t) const
{
    double cosX = std::cos(t * this->b), sinX = std::sin(t * this->b);
    double cosY = std::cos(t * this->a), sinY = std::sin(t * this->a);
    std::complex<double> numerator(cosX - cosY, sinX - sinY);
    std::complex<double> denominator(0, t * this->bma);
    return numerator / denominator;
}

template < typename RealType >
long double UniformRand<RealType>::Entropy() const
{
    return (this->b == this->a) ? -INFINITY : std::log(this->bma);
}

template < typename RealType >
double UniformRand<RealType>::LikelihoodFunction(const std::vector<RealType> &sample) const
{
    bool sampleIsInsideInterval = this->allElementsAreNotSmallerThan(this->a, sample) && this->allElementsAreNotGreaterThan(this->b, sample);
    return sampleIsInsideInterval ? std::pow(this->bma, -sample.size()) : 0.0;
}

template < typename RealType >
double UniformRand<RealType>::LogLikelihoodFunction(const std::vector<RealType> &sample) const
{
    bool sampleIsInsideInterval = this->allElementsAreNotSmallerThan(this->a, sample) && this->allElementsAreNotGreaterThan(this->b, sample);
    int sample_size = sample.size();
    return sampleIsInsideInterval ? -sample_size * this->logbma : -INFINITY;
}

template < typename RealType >
constexpr char UniformRand<RealType>::TOO_LARGE_A[];
template < typename RealType >
constexpr char UniformRand<RealType>::TOO_SMALL_B[];

template < typename RealType >
void UniformRand<RealType>::FitMinimum(const std::vector<RealType> &sample, bool unbiased)
{
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b)));
    RealType minVar = *std::min_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[min] = b - n / (n + 1) * (this->b - this->a)
        RealType minVarAdj = (minVar * (n + 1) - this->b) / n;
        if (!this->allElementsAreNotSmallerThan(minVarAdj, sample))
            throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, TOO_LARGE_A + this->toStringWithPrecision(minVarAdj)));
        SetSupport(minVarAdj, this->b);
    }
    else {
        SetSupport(minVar, this->b);
    }
}

template < typename RealType >
void UniformRand<RealType>::FitMaximum(const std::vector<RealType> &sample, bool unbiased)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    RealType maxVar = *std::max_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[max] = (this->b - this->a) * n / (n + 1) + a
        RealType maxVarAdj = (maxVar * (n + 1) - this->a) / n;
        if (!this->allElementsAreNotGreaterThan(maxVarAdj, sample))
            throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, TOO_SMALL_B + this->toStringWithPrecision(maxVarAdj)));
        SetSupport(this->a, maxVarAdj);
    }
    else {
        SetSupport(this->a, maxVar);
    }
}

template < typename RealType >
void UniformRand<RealType>::Fit(const std::vector<RealType> &sample, bool unbiased)
{
    double minVar = *std::min_element(sample.begin(), sample.end());
    double maxVar = *std::max_element(sample.begin(), sample.end());
    if (unbiased == true) {
        int n = sample.size();
        /// E[min] = b - n / (n + 1) * (this->b - this->a)
        RealType minVarAdj = (minVar * n - maxVar) / (n - 1);
        /// E[max] = (this->b - this->a) * n / (n + 1) + a
        RealType maxVarAdj = (maxVar * n - minVar) / (n - 1);
        if (!this->allElementsAreNotSmallerThan(minVarAdj, sample))
            throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, TOO_LARGE_A + this->toStringWithPrecision(minVarAdj)));
        if (!this->allElementsAreNotGreaterThan(maxVarAdj, sample))
            throw std::runtime_error(this->fitErrorDescription(this->WRONG_RETURN, TOO_SMALL_B + this->toStringWithPrecision(maxVarAdj)));
        SetSupport(minVarAdj, maxVarAdj);
    }
    else {
        SetSupport(minVar, maxVar);
    }
}

template < typename RealType >
ParetoRand<RealType> UniformRand<RealType>::FitMaximumBayes(const std::vector<RealType> &sample, const ParetoRand<RealType> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw std::invalid_argument(this->fitErrorDescription(this->WRONG_SAMPLE, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a)));
    double maxVar = *std::max_element(sample.begin(), sample.end());
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newScale = std::max(priorDistribution.GetScale(), maxVar - this->a);
    ParetoRand<RealType> posteriorDistribution(newShape, newScale);
    double theta = MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean();
    SetSupport(this->a, this->a + theta);
    return posteriorDistribution;
}

template class UniformRand<float>;
template class UniformRand<double>;
template class UniformRand<long double>;
