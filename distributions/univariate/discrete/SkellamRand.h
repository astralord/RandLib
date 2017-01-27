#ifndef SKELLAMRAND_H
#define SKELLAMRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"

/**
 * @brief The SkellamRand class
 * Skellam distribution
 *
 * Notation: X ~ Skellam(μ_1, μ_2)
 *
 * Related distributions:
 * If Y ~ Poisson(μ_1) and Z ~ Poisson(μ_2) then Y - Z ~ Skellam(μ_1, μ_2)
 */
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution
{
    double mu1, mu2;
    double pmfCoef1; /// 2 * sqrt(mu1 * mu2)
    double pmfCoef2; /// log(mu1 / mu2)

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

    double P(int k) const override;
    double F(int k) const override;
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
