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

    double denominator; /// 1 / harmonic_number

    // TODO: replace to RandMath
    static double harmonicNumber(double exponent, size_t number);
public:
    ZipfRand(double exponent, size_t number);
    virtual std::string name() override;

    void setParameters(double exponent, size_t number);
    inline double getExponent() { return s; }
    inline size_t getNumber() { return N; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return harmonicNumber(s - 1, N) * denominator; }
    double Var() const override {
        double numerator = harmonicNumber(s - 1, N);
        numerator *= numerator;
        numerator = harmonicNumber(s - 2, N) * harmonicNumber(s, N) - numerator;
        return numerator * denominator * denominator;
    }
};

#endif // ZIPFRAND_H
