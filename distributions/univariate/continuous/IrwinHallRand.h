#ifndef IRWINHALLRAND_H
#define IRWINHALLRAND_H

#include "UniformRand.h"

/**
 * @brief The IrwinHallRand class
 * Irwin-Hall distribution
 * X ~ IH(n)
 *
 * f(x|n) = 0.5 / (n - 1)! * sum_{k=0}^n (-1)^k * C(n,k) * (x - k) ^ (n - 1) * sgn(x - k)
 *
 * X ~ Y_1 + Y_2 + ... + Y_n, where Y_i ~ U(0,1)
 */
class RANDLIBSHARED_EXPORT IrwinHallRand : public ContinuousDistribution
{
    int n;
    UniformRand U;
    double pdfCoef, cdfCoef;
public:
    explicit IrwinHallRand(int number);
    std::string name() override;

    void setNumber(int number);
    inline int getNumber();

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // IRWINHALLRAND_H
