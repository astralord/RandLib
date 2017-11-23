#ifndef SKELLAMRAND_H
#define SKELLAMRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"

/**
 * @brief The SkellamRand class <BR>
 * Skellam distribution
 *
 * Notation: X ~ Skellam(μ1, μ2)
 *
 * Related distributions: <BR>
 * If Y ~ Po(μ1) and Z ~ Po(μ2) then Y - Z ~ Skellam(μ1, μ2)
 */
template < typename IntType = int >
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution<IntType>
{
    double mu1 = 1; ///< first rate μ1
    double mu2 = 1; ///< second rate μ2
    double logMu1 = 0; ///< log(μ1)
    double logMu2 = 0; ///< log(μ2)
    double sqrtMu1 = 1; ///< √μ1
    double sqrtMu2 = 1; ///< √μ2

    PoissonRand<IntType> X{}, Y{};

public:
    SkellamRand(double rate1 = 1.0, double rate2 = 1.0);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    IntType MinValue() const override { return std::numeric_limits<IntType>::lowest(); }
    IntType MaxValue() const override { return std::numeric_limits<IntType>::max(); }

    void SetRates(double rate1, double rate2);
    inline double GetFirstRate() const { return mu1; }
    inline double GetSecondRate() const { return mu2; }

    double P(const IntType & k) const override;
    double logP(const IntType & k) const override;
    double F(const IntType & k) const override;
    double S(const IntType & k) const override;
    IntType Variate() const override;
    void Sample(std::vector<IntType> &outputData) const override;
    void Reseed(unsigned long seed) const override;

    long double Mean() const override;
    long double Variance() const override;
    IntType Median() const override;
    IntType Mode() const override;
    long double Skewness() const override;
    long double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // SKELLAMRAND_H
