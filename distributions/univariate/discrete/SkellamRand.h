#ifndef SKELLAMRAND_H
#define SKELLAMRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"

/**
 * @brief The SkellamRand class
 * Skellam distribution
 * X ~ Skellam(μ_1, μ_2)
 *
 * If Y ~ Poisson(μ_1) and Z ~ Poisson(μ_2) then Y - Z ~ Skellam(μ_1, μ_2)
 */
// TODO: define it as a compound distribution
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution
{
    double mu1, mu2;
    double pmfCoef1; /// exp(-(mu1 + mu2))
    double pmfCoef2; /// sqrt(mu1 / mu2)
    double pmfCoef3; /// 2 * sqrt(mu1 * mu2)

    PoissonRand X, Y;

public:
    SkellamRand(double mean1, double mean2);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    int MinValue() const override { return INT_MIN; }
    int MaxValue() const override { return INT_MAX; }

    void setMeans(double mean1, double mean2);
    inline double getFirstMean() const { return mu1; }
    inline double getSecondMean() const { return mu2; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    void sample(std::vector<int> &outputData) const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // SKELLAMRAND_H
