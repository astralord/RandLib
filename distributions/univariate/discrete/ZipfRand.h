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
    std::string name() override;

    void setParameters(double exponent, int number);
    inline double getExponent() { return s; }
    inline size_t getNumber() { return n; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
    
    std::complex<double> CF(double t) const override;
     
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // ZIPFRAND_H
