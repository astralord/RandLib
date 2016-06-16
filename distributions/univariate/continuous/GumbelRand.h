#ifndef GUMBELRAND_H
#define GUMBELRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The GumbelRand class
 */
class RANDLIBSHARED_EXPORT GumbelRand : public ContinuousDistribution
{
    double mu, beta;
    double betaInv;
public:
    GumbelRand(double location, double scale);
    std::string name() override;
    SUPPORT_TYPE supportType() const override { return INFINITE_T; }
    double MinValue() const override { return -INFINITY; }
    double MaxValue() const override { return INFINITY; }

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return beta; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;
    static double variate(double location, double scale);

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // GUMBELRAND_H
