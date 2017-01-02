#ifndef GUMBELRAND_H
#define GUMBELRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GumbelRand class
 * Gumbel distribution
 */
class RANDLIBSHARED_EXPORT GumbelRand : public ContinuousDistribution
{
    double mu, beta;
    double betaInv;
public:
    GumbelRand(double location, double scale);

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void SetLocation(double location);
    void SetScale(double scale);
    inline double GetLocation() const { return mu; }
    inline double GetScale() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;
    static double Variate(double location, double scale);

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

private:
    double quantileImpl(double p) const override;
    double quantileImpl1m(double p) const override;

public:
    double Entropy() const;
};

#endif // GUMBELRAND_H
