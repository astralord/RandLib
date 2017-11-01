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
class RANDLIBSHARED_EXPORT ZipfRand : public DiscreteDistribution
{
    double s = 1; ///< exponent
    int n = 1; ///< number
    double invHarmonicNumber = 1; /// 1 / H(s, n)

    static constexpr int tableSize = 16;
    int hashedVarNum = 1;
    double table[tableSize];

public:
    ZipfRand(double exponent, int number);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return n; }

    void SetParameters(double exponent, int number);
    inline double GetExponent() const { return s; }
    inline int GetNumber() const { return n; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    int Variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // ZIPFRAND_H
