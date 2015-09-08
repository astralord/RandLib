#ifndef ZIPFRAND_H
#define ZIPFRAND_H

#include "DiscreteRand.h"
#include "../continuous/UniformRand.h"

/**
 * @brief The ZipfRand class
 */
class ZipfRand : public DiscreteRand<int>
{
    double s;
    size_t N;

    double invHarmonicNumber; /// 1 / harmonic_number

public:
    ZipfRand(double exponent, size_t number);
    virtual std::string name() override;

    void setParameters(double exponent, size_t number);
    inline double getExponent() { return s; }
    inline size_t getNumber() { return N; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return RandMath::harmonicNumber(s - 1, N) * invHarmonicNumber; }
    double Var() const override {
        double numerator = RandMath::harmonicNumber(s - 1, N);
        numerator *= numerator;
        numerator = RandMath::harmonicNumber(s - 2, N) * RandMath::harmonicNumber(s, N) - numerator;
        return numerator * invHarmonicNumber * invHarmonicNumber;
    }
};

#endif // ZIPFRAND_H
