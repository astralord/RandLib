#ifndef IRWINHALLRAND_H
#define IRWINHALLRAND_H

#include "UniformRand.h"

/**
 * @brief The IrwinHallRand class
 *
 * f(x|n) = 0.5 / (n - 1)! * sum_{k=0}^n (-1)^k * C(n,k) * (x - k) ^ (n - 1) * sgn(x - k)
 *
 * The sum of a number of independent random variables, each having a uniform distribution on [0, 1]
 */
class RANDLIBSHARED_EXPORT IrwinHallRand : public ContinuousRand
{
    UniformRand U;
    double pdfCoef, cdfCoef;
    int n;
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

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // IRWINHALLRAND_H
