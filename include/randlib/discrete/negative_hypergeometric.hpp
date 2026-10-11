#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The NegativeHyperGeometricRand class <BR>
 * Negative hypergeometric distribution
 *
 * Notation: X ~ NHG(N, M, m)
 */
template <typename IntType>
class NegativeHyperGeometricRand final : public DiscreteDistribution<IntType>
{
    IntType N = 1; ///< size of population
    IntType M = 1; ///< total amount of successes
    IntType m = 1; ///< limiting number of successes
    double pmfCoef = 1; ///< C(N, M)
    double p0 = 1; ///< M / N

public:
    NegativeHyperGeometricRand(IntType totalSize = 1, IntType totalSuccessesNum = 1, IntType limitSuccessesNum = 1);
    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return 0; }
    [[nodiscard]] IntType MaxValue() const override { return N - M; }

    void SetParameters(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum);
    [[nodiscard]] inline IntType GetTotalSize() { return N; }
    [[nodiscard]] inline IntType GetTotalSuccessesNum() { return M; }
    [[nodiscard]] inline IntType GetLimitSuccessesNum() { return m; }

    [[nodiscard]] double P(const IntType & k) const override;
    [[nodiscard]] double logP(const IntType & k) const override;
    [[nodiscard]] double F(const IntType & k) const override;
    template <std::uniform_random_bit_generator G>
    [[nodiscard]] IntType operator()(G &g) const;
protected:
    IntType variate(detail::EngineRef g) const override { return (*this)(g); }
public:

    [[nodiscard]] long double Mean() const override;
    [[nodiscard]] long double Variance() const override;
};

} // namespace randlib

#include <randlib/discrete/bernoulli.hpp>

namespace randlib
{

template < typename IntType >
NegativeHyperGeometricRand<IntType>::NegativeHyperGeometricRand(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum)
{
    SetParameters(totalSize, totalSuccessesNum, limitSuccessesNum);
}

template < typename IntType >
std::string NegativeHyperGeometricRand<IntType>::Name() const
{
    return "Negative hypergeometric(" + this->toStringWithPrecision(N) + ", "
                                      + this->toStringWithPrecision(M) + ", "
                                      + this->toStringWithPrecision(m) + ")";
}

template < typename IntType >
void NegativeHyperGeometricRand<IntType>::SetParameters(IntType totalSize, IntType totalSuccessesNum, IntType limitSuccessesNum)
{
    if (!(totalSize > 0 && totalSuccessesNum > 0 && limitSuccessesNum > 0))
        throw std::invalid_argument("Negative-HyperGeometric distribution: all parameters should be positive");
    if (!(totalSuccessesNum <= totalSize))
        throw std::invalid_argument("Negative-HyperGeometric distribution: total size shouldn't be smaller than total successes number");
    if (!(limitSuccessesNum <= totalSuccessesNum))
        throw std::invalid_argument("Negative-HyperGeometric distribution: total successes number shouldn't be smaller than limit successes number");

    N = totalSize;
    M = totalSuccessesNum;
    m = limitSuccessesNum;

    p0 = static_cast<double>(M) / N;
    pmfCoef = math::lfact(M);
    pmfCoef += math::lfact(N - M);
    pmfCoef -= math::lfact(m - 1);
    pmfCoef -= math::lfact(M - m);
    pmfCoef -= math::lfact(N);
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::P(const IntType & k) const
{
    return (k < MinValue() || k > MaxValue()) ? 0.0 : std::exp(logP(k));
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::logP(const IntType & k) const
{
    if (k < MinValue() || k > MaxValue())
        return -INFINITY;
    double p = math::lfact(k + m - 1);
    p += math::lfact(N - m - k);
    p -= math::lfact(k);
    p -= math::lfact(N - M - k);
    return p + pmfCoef;
}

template < typename IntType >
double NegativeHyperGeometricRand<IntType>::F(const IntType &k) const
{
    // relation with hypergeometric distribution can be used here instead
    if (k < MinValue())
        return 0.0;
    IntType maxVal = MaxValue();
    if (k >= maxVal)
        return 1.0;
    if (k <= 0.5 * maxVal) {
        /// sum P(X = i) going forward until k
        double sum = 0;
        for (IntType i = 0; i <= k; ++i)
            sum += P(i);
        return sum;
    }
    /// going backwards is faster
    double sum = 1.0;
    for (IntType i = k + 1; i <= maxVal; ++i)
        sum -= P(i);
    return sum;
}

template < typename IntType >
template <std::uniform_random_bit_generator G>
IntType NegativeHyperGeometricRand<IntType>::operator()(G &g) const
{
    double p = p0;
    IntType successesNum = 0;
    IntType num = 0;
    while (successesNum < m) {
        ++num;
        if (BernoulliRand::Variate(p, g) && ++successesNum == num - N + M)
            return N - M;
        p = M - successesNum;
        p /= N - num;
    }
    return num - successesNum;
}

template < typename IntType >
long double NegativeHyperGeometricRand<IntType>::Mean() const
{
    double mean = m;
    mean *= N - M;
    mean /= M + 1;
    return mean;
}

template < typename IntType >
long double NegativeHyperGeometricRand<IntType>::Variance() const
{
    double Mp1 = M + 1;
    double var = 1 - m / Mp1;
    var *= N + 1;
    var *= N - M;
    var /= Mp1 * (Mp1 + 1);
    return m * var;
}

} // namespace randlib
