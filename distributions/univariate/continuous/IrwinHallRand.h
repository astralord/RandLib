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
template < typename RealType = double >
class RANDLIBSHARED_EXPORT IrwinHallRand : public ContinuousDistribution<RealType>
{
    int n = 1; ///< parameter of the distribution
    UniformRand<RealType> U{};

public:
    explicit IrwinHallRand(size_t number);

    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return FINITE_T; }
    RealType MinValue() const override { return 0; }
    RealType MaxValue() const override { return n; }

    void SetNumber(int number);
    inline int GetNumber() const { return n; }

    double f(const RealType & x) const override;
    double logf(const RealType & x) const override;
    double F(const RealType & x) const override;
    RealType Variate() const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    RealType Median() const override;
    RealType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // IRWINHALLRAND_H
