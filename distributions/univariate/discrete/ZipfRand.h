#ifndef ZIPFRAND_H
#define ZIPFRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The ZipfRand class <BR>
 * Zipf distribution
 *
 * Notation: X ~ Zipf(s, n)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT ZipfRand : public DiscreteDistribution<IntType>
{
    double s = 2; ///< exponent
    IntType n = 1; ///< number
    double invHarmonicNumber = 1; /// 1 / H(s, n)

    static constexpr int tableSize = 16;
    int hashedVarNum = 1;
    double table[tableSize];

public:
    ZipfRand(double exponent = 2.0, IntType number = 1);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    IntType MinValue() const override { return 1; }
    IntType MaxValue() const override { return n; }

    void SetParameters(double exponent, IntType number);
    inline double GetExponent() const { return s; }
    inline IntType GetNumber() const { return n; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    IntType Variate() const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // ZIPFRAND_H
