#ifndef SKELLAMRAND_H
#define SKELLAMRAND_H

#include "DiscreteDistribution.h"
#include "PoissonRand.h"

/**
 * @brief The SkellamRand class
 * Skellam distribution
 * X ~ Skellam(\mu_1, \mu_2)
 *
 * If Y ~ Poisson(\mu_1) and Z ~ Poisson(\mu_2) then Y - Z ~ Skellam(\mu_1, \mu_2)
 */
class RANDLIBSHARED_EXPORT SkellamRand : public DiscreteDistribution
{
    double mu1, mu2;
    double pmfCoef1; /// exp(-(mu1 + mu2))
    double pmfCoef2; /// sqrt(mu1 / mu2)
    double pmfCoef3; /// 2 * sqrt(mu1 * mu2)

    PoissonRand X, Y;

public:
    SkellamRand(double mean1, double mean2);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void setMeans(double mean1, double mean2);
    inline double getFirstMean() { return mu1; }
    inline double getSecondMean() { return mu2; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // SKELLAMRAND_H
