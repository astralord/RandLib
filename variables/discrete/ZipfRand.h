#ifndef ZIPFRAND_H
#define ZIPFRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The ZipfRand class
 */
class RANDLIBSHARED_EXPORT ZipfRand : public DiscreteRand
{
    double s;
    int N;

    double invHarmonicNumber; /// 1 / harmonic_number

public:
    ZipfRand(double exponent, int number);
    std::string name() override;

    void setParameters(double exponent, int number);
    inline double getExponent() { return s; }
    inline size_t getNumber() { return N; }

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
