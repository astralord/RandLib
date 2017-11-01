#ifndef IRWINHALLRAND_H
#define IRWINHALLRAND_H

#include "UniformRand.h"

/**
 * @brief The IrwinHallRand class <BR>
 * Irwin-Hall distribution
 *
 * f(x | n) = 0.5 / (n - 1)! * sum_{k=0}^n (-1)^k * C(n,k) * (x - k) ^ (n - 1) * sign(x - k)
 *
 * Notation: X ~ IH(n)
 *
 * Related distributions: <BR>
 * X ~ Y_1 + Y_2 + ... + Y_n, where Y_i ~ U(0,1)
 */
class RANDLIBSHARED_EXPORT IrwinHallRand : public ContinuousDistribution
{
    int n = 1; ///< parameter of the distribution
    UniformRand U{};
public:
    explicit IrwinHallRand(size_t number);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    double MinValue() const override { return 0; }
    double MaxValue() const override { return n; }

    void SetNumber(int number);
    inline int GetNumber() const { return n; }

    double f(const double & x) const override;
    double logf(const double & x) const override;
    double F(const double & x) const override;
    double Variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // IRWINHALLRAND_H
