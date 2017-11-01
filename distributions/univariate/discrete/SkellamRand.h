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
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution
{
    double mu1 = 1; ///< first rate μ1
    double mu2 = 1; ///< second rate μ2
    double logMu1 = 0; ///< log(μ1)
    double logMu2 = 0; ///< log(μ2)
    double sqrtMu1 = 1; ///< √μ1
    double sqrtMu2 = 1; ///< √μ2

    PoissonRand X{}, Y{};

public:
    SkellamRand(double rate1, double rate2);
    String Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    int MinValue() const override { return INT_MIN; }
    int MaxValue() const override { return INT_MAX; }

    void SetRates(double rate1, double rate2);
    inline double GetFirstRate() const { return mu1; }
    inline double GetSecondRate() const { return mu2; }

    double P(const int & k) const override;
    double logP(const int & k) const override;
    double F(const int & k) const override;
    double S(const int & k) const override;
    int Variate() const override;
    void Sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    std::complex<double> CFImpl(double t) const override;
};

#endif // SKELLAMRAND_H
