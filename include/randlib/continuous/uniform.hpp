#pragma once

#include <span>

#include <algorithm>
#include <limits>

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/continuous/beta.hpp>

namespace randlib
{

/**
 * @brief The UniformRand class <BR>
 * Uniform continuous distribution
 *
 * f(x | a, b) = 1 / (b - a) for a < x < b
 *
 * Notation: X ~ U(a, b)
 *
 * Related distributions: <BR>
 * X ~ B(1, 1, a, b) <BR>
 * (X - a) / (b - a) ~ IH(1)
 */
template <typename RealType>
class UniformRand final : public BetaDistribution<RealType>
{
public:
    UniformRand(double minValue = 0, double maxValue = 1);
    [[nodiscard]] std::string Name() const override;

    using BetaDistribution<RealType>::SetSupport;

    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] RealType MinValue() const override { return this->a; }
    [[nodiscard]] RealType MaxValue() const override { return this->b; }

    [[nodiscard]] double f(const RealType & x) const override;
    [[nodiscard]] double logf(const RealType & x) const override;
    [[nodiscard]] double F(const RealType & x) const override;
    [[nodiscard]] double S(const RealType & x) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] RealType operator()(G &g) const;
protected:
    RealType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    /**
     * @fn StandardVariate
     * @param g
     * @return a random number on interval (0,1) if no preprocessors are specified
     */
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariate(G &g);

    /**
     * @fn StandardVariateClosed
     * @param g
     * @return a random number on interval [0,1]
     */
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariateClosed(G &g);

    /**
     * @fn StandardVariateHalfClosed
     * @param g
     * @return a random number on interval [0,1)
     */
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] static RealType StandardVariateHalfClosed(G &g);

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
    [[nodiscard]] RealType Median() const override;
    [[nodiscard]] RealType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    RealType quantileImpl(double p) const override;
    RealType quantileImpl1m(double p) const override;

    std::complex<double> CFImpl(double t) const override;

    static constexpr char TOO_LARGE_A[] = "Minimum element of the sample is smaller than lower boundary returned by method: ";
    static constexpr char TOO_SMALL_B[] = "Maximum element of the sample is greater than upper boundary returned by method: ";

public:
    [[nodiscard]] long double Entropy() const;

    [[nodiscard]] double LikelihoodFunction(std::span<const RealType> sample) const override;
    [[nodiscard]] double LogLikelihoodFunction(std::span<const RealType> sample) const override;

    /**
     * @fn FitMinimum
     * fit minimum with maximum-likelihood estimator if unbiased == false,
     * fit minimum using UMVU estimator otherwise
     * @param sample
     */
    void FitMinimum(std::span<const RealType> sample, bool unbiased = false);
    /**
     * @fn FitMaximum
     * fit maximum with maximum-likelihood estimator if unbiased == false,
     * fit maximum using UMVU estimator otherwise
     * @param sample
     */
    void FitMaximum(std::span<const RealType> sample, bool unbiased = false);
    /**
     * @fn Fit
     * fit support with maximum-likelihood estimator if unbiased == false,
     * fit support using UMVU estimator otherwise
     * @param sample
     */
    void Fit(std::span<const RealType> sample, bool unbiased = false);

    /**
     * @fn FitMaximumBayes
     * fit maximum parameter, using Bayesian estimation
     * @param sample
     * @param priorDistribution
     * @param MAP if true, use MAP estimator
     * @return posterior distribution
     */
    ParetoRand<RealType> FitMaximumBayes(std::span<const RealType> sample, const ParetoRand<RealType> &priorDistribution, bool MAP = false);
};

} // namespace randlib

#include <randlib/continuous/pareto.hpp>

namespace randlib
{

template < typename RealType >
UniformRand<RealType>::UniformRand(double minValue, double maxValue) :
    BetaDistribution<RealType>(1, 1, minValue, maxValue)
{
}

template < typename RealType >
std::string UniformRand<RealType>::Name() const
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
template <std::uniform_random_bit_generator G>
RealType UniformRand<RealType>::operator()(G &g) const
{
  return this->a + StandardVariate(g) * this->bma;
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType UniformRand<RealType>::StandardVariate(G &g)
{
  return detail::uniform01<RealType>(g);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType UniformRand<RealType>::StandardVariateClosed(G &g)
{
  return detail::uniform01Closed<RealType>(g);
}

template < typename RealType >
template <std::uniform_random_bit_generator G>
RealType UniformRand<RealType>::StandardVariateHalfClosed(G &g)
{
  return detail::uniform01HalfClosed<RealType>(g);
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
double UniformRand<RealType>::LikelihoodFunction(std::span<const RealType> sample) const
{
    bool sampleIsInsideInterval = this->allElementsAreNotSmallerThan(this->a, sample) && this->allElementsAreNotGreaterThan(this->b, sample);
  return sampleIsInsideInterval ? std::pow(this->bma, -sample.size()) : 0.0;
}

template < typename RealType >
double UniformRand<RealType>::LogLikelihoodFunction(std::span<const RealType> sample) const
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
void UniformRand<RealType>::FitMinimum(std::span<const RealType> sample, bool unbiased)
{
    if (!this->allElementsAreNotGreaterThan(this->b, sample))
        throw this->fitError(FitErrorCode::WrongSample, this->UPPER_LIMIT_VIOLATION + this->toStringWithPrecision(this->b));
    RealType minVar = *std::min_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[min] = b - n / (n + 1) * (this->b - this->a)
        RealType minVarAdj = (minVar * (n + 1) - this->b) / n;
        if (!this->allElementsAreNotSmallerThan(minVarAdj, sample))
            throw this->fitError(FitErrorCode::WrongReturn, TOO_LARGE_A + this->toStringWithPrecision(minVarAdj));
        SetSupport(minVarAdj, this->b);
    }
    else {
        SetSupport(minVar, this->b);
    }
}

template < typename RealType >
void UniformRand<RealType>::FitMaximum(std::span<const RealType> sample, bool unbiased)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    RealType maxVar = *std::max_element(sample.begin(), sample.end());

    if (unbiased == true) {
        int n = sample.size();
        /// E[max] = (this->b - this->a) * n / (n + 1) + a
        RealType maxVarAdj = (maxVar * (n + 1) - this->a) / n;
        if (!this->allElementsAreNotGreaterThan(maxVarAdj, sample))
            throw this->fitError(FitErrorCode::WrongReturn, TOO_SMALL_B + this->toStringWithPrecision(maxVarAdj));
        SetSupport(this->a, maxVarAdj);
    }
    else {
        SetSupport(this->a, maxVar);
    }
}

template < typename RealType >
void UniformRand<RealType>::Fit(std::span<const RealType> sample, bool unbiased)
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
            throw this->fitError(FitErrorCode::WrongReturn, TOO_LARGE_A + this->toStringWithPrecision(minVarAdj));
        if (!this->allElementsAreNotGreaterThan(maxVarAdj, sample))
            throw this->fitError(FitErrorCode::WrongReturn, TOO_SMALL_B + this->toStringWithPrecision(maxVarAdj));
        SetSupport(minVarAdj, maxVarAdj);
    }
    else {
        SetSupport(minVar, maxVar);
    }
}

template < typename RealType >
ParetoRand<RealType> UniformRand<RealType>::FitMaximumBayes(std::span<const RealType> sample, const ParetoRand<RealType> &priorDistribution, bool MAP)
{
    if (!this->allElementsAreNotSmallerThan(this->a, sample))
        throw this->fitError(FitErrorCode::WrongSample, this->LOWER_LIMIT_VIOLATION + this->toStringWithPrecision(this->a));
    double maxVar = *std::max_element(sample.begin(), sample.end());
    int n = sample.size();
    double newShape = priorDistribution.GetShape() + n;
    double newScale = std::max(priorDistribution.GetScale(), maxVar - this->a);
    ParetoRand<RealType> posteriorDistribution(newShape, newScale);
    double theta = MAP ? posteriorDistribution.Mode() : posteriorDistribution.Mean();
    SetSupport(this->a, this->a + theta);
  return posteriorDistribution;
}

} // namespace randlib
