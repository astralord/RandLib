#ifndef ZIPFRAND_H
#define ZIPFRAND_H

#include "DiscreteDistribution.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The ZipfRand class
 * Zipf distribution
 * X ~ Zipf(s, n)
 */
class RANDLIBSHARED_EXPORT ZipfRand : public DiscreteDistribution
{
    double s;
    int n;

    double invHarmonicNumber; /// 1 / harmonic_number

    static constexpr int tableSize = 16;
    int hashedVarNum;
    double table[tableSize];

public:
    ZipfRand(double exponent, int number);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return n; }

    void setParameters(double exponent, int number);
    inline double getExponent() const { return s; }
    inline int getNumber() const { return n; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;
};

#endif // ZIPFRAND_H
