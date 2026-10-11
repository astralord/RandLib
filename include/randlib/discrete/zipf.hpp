#pragma once

#include <randlib/math/constants.hpp>
#include <randlib/math/functions.hpp>
#include <randlib/core/errors.hpp>
#include <randlib/core/discrete.hpp>

namespace randlib
{

/**
 * @brief The ZipfRand class <BR>
 * Zipf distribution
 *
 * Notation: X ~ Zipf(s, n)
 */
template <typename IntType>
class ZipfRand final : public DiscreteDistribution<IntType>
{
    double s = 2; ///< exponent
    IntType n = 1; ///< number
    double invHarmonicNumber = 1; /// 1 / H(s, n)

    static constexpr int tableSize = 16;
    int hashedVarNum = 1;
    double table[tableSize];

public:
    ZipfRand(double exponent = 2.0, IntType number = 1);

    [[nodiscard]] std::string Name() const override;
    [[nodiscard]] Support SupportType() const override { return Support::Finite; }
    [[nodiscard]] IntType MinValue() const override { return 1; }
    [[nodiscard]] IntType MaxValue() const override { return n; }

    void SetParameters(double exponent, IntType number);
    [[nodiscard]] inline double GetExponent() const noexcept { return s; }
    [[nodiscard]] inline IntType GetNumber() const noexcept { return n; }

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
    [[nodiscard]] IntType Mode() const override;
    [[nodiscard]] long double Skewness() const override;
    [[nodiscard]] long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

} // namespace randlib

#include <randlib/continuous/uniform.hpp>

namespace randlib
{

template < typename IntType >
ZipfRand<IntType>::ZipfRand(double exponent, IntType number)
{
    SetParameters(exponent, number);
}

template < typename IntType >
std::string ZipfRand<IntType>::Name() const
{
  return "Zipf(" + this->toStringWithPrecision(GetExponent()) + ", "
                   + this->toStringWithPrecision(GetNumber()) + ")";
}

template < typename IntType >
void ZipfRand<IntType>::SetParameters(double exponent, IntType number)
{
    if (!(exponent > 1.0))
        throw std::invalid_argument("Zipf distribution: exponent should be greater than 1, but it's equal to "
                                    + std::to_string(exponent));
    if (!(number > 0))
        throw std::invalid_argument("Zipf distribution: number should be positive, but it's equal to "
                                    + std::to_string(number));
    s = exponent;
    n = number;

    invHarmonicNumber = 1.0 / math::harmonicNumber(s, n);

    // WARNING: we calculate pow here and in invHarmonic number
    hashedVarNum = n > tableSize ? tableSize : n;
    table[0] = 1.0;
    for (int i = 1; i < hashedVarNum - 1; ++i)
        table[i] = table[i - 1] + std::pow(i + 1, -s);
    if (hashedVarNum == n)
        table[hashedVarNum - 1] = 1.0 / invHarmonicNumber;
    else
        table[hashedVarNum - 1] = table[hashedVarNum - 2] + std::pow(hashedVarNum, -s);
    for (int i = 0; i < hashedVarNum; ++i)
        table[i] *= invHarmonicNumber;
}

template < typename IntType >
double ZipfRand<IntType>::P(const IntType &k) const
{
  return (k < 1 || k > n) ? 0.0 : std::pow(k, -s) * invHarmonicNumber;
}

template < typename IntType >
double ZipfRand<IntType>::logP(const IntType & k) const
{
  return (k < 1 || k > n) ? -INFINITY : -s * std::log(k) + std::log(invHarmonicNumber); // can be hashed
}

template < typename IntType >
double ZipfRand<IntType>::F(const IntType & k) const
{
    if (k < 1.0)
      return 0.0;
    if (k >= n)
      return 1.0;
  return math::harmonicNumber(s, k) * invHarmonicNumber;
}

template < typename IntType >
template <std::uniform_random_bit_generator G>
IntType ZipfRand<IntType>::operator()(G &g) const
{
    double U = detail::uniform01<double>(g);
    int k = 1;
    /// if we didn't manage to hash values for such U
    if (U > table[hashedVarNum - 1]) {
        k = hashedVarNum;
        double sum = table[hashedVarNum - 1];
        do {
            ++k;
            sum += std::pow(k, -s) * invHarmonicNumber;
        } while (sum < U);
    }
    else {
        while (k <= hashedVarNum && table[k - 1] < U)
            ++k;
    }
  return k;
}

template < typename IntType >
long double ZipfRand<IntType>::Mean() const
{
  return math::harmonicNumber(s - 1, n) * invHarmonicNumber;
}

template < typename IntType >
long double ZipfRand<IntType>::Variance() const
{
    double numerator = math::harmonicNumber(s - 1, n);
    numerator *= numerator;
    numerator = math::harmonicNumber(s - 2, n) * math::harmonicNumber(s, n) - numerator;
  return numerator * invHarmonicNumber * invHarmonicNumber;
}

template < typename IntType >
IntType ZipfRand<IntType>::Mode() const
{
  return 1;
}

template < typename IntType >
std::complex<double> ZipfRand<IntType>::CFImpl(double t) const
{
    std::complex<double> sum(0.0, 0.0);
    for (int i = 1; i <= n; ++i)
    {
        std::complex<double> addon(-s * std::log(i), i * t);
        sum += std::exp(addon);
    }
  return invHarmonicNumber * sum;
}

template < typename IntType >
long double ZipfRand<IntType>::Skewness() const
{
    long double harmonic0 = 1.0 / invHarmonicNumber;
    long double harmonic1 = math::harmonicNumber(s - 1, n);
    long double harmonic2 = math::harmonicNumber(s - 2, n);
    long double harmonic3 = math::harmonicNumber(s - 3, n);
    long double first = harmonic3 * harmonic0 * harmonic0;
    long double harmonic2prod0 = harmonic2 * harmonic0;
    long double second = -3 * harmonic1 * harmonic2prod0;
    long double harmonic1Sq = harmonic1 * harmonic1;
    long double third = 2 * harmonic1 * harmonic1Sq;
    long double numerator = first + second + third;
    long double denominator = harmonic2prod0 - harmonic1Sq;
    denominator *= std::sqrt(denominator);
  return numerator / denominator;
}

template < typename IntType >
long double ZipfRand<IntType>::ExcessKurtosis() const
{
    long double harmonic0 = 1.0 / invHarmonicNumber;
    long double harmonic1 = math::harmonicNumber(s - 1, n);
    long double harmonic2 = math::harmonicNumber(s - 2, n);
    long double harmonic3 = math::harmonicNumber(s - 3, n);
    long double harmonic4 = math::harmonicNumber(s - 4, n);
    long double harmonic2prod0 = harmonic2 * harmonic0;
    long double harmonic0Sq = harmonic0 * harmonic0;
    long double harmonic1Sq = harmonic1 * harmonic1;
    long double denominator = harmonic2prod0 - harmonic1Sq;
    denominator *= denominator;
    long double numerator = harmonic0Sq * harmonic0 * harmonic4;
    numerator -= 4 * harmonic0Sq * harmonic1 * harmonic3;
    numerator += 6 * harmonic2prod0 * harmonic1Sq;
    numerator -= 3 * harmonic1Sq * harmonic1Sq;
  return numerator / denominator - 3;
}

} // namespace randlib
