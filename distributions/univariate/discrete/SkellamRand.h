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
 * Related distributions:
 * If Y ~ Po(μ1) and Z ~ Po(μ2) then Y - Z ~ Skellam(μ1, μ2)
 */
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution
{
    double mu1, mu2;
    double logMu1, logMu2;
    double sqrtMu1, sqrtMu2;

    PoissonRand X, Y;

public:
    SkellamRand(double mean1, double mean2);
    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    int MinValue() const override { return INT_MIN; }
    int MaxValue() const override { return INT_MAX; }

    void SetMeans(double mean1, double mean2);
    inline double GetFirstMean() const { return mu1; }
    inline double GetSecondMean() const { return mu2; }

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
