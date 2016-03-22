#ifndef FRECHETRAND_H
#define FRECHETRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The FrechetRand class
 */
class RANDLIBSHARED_EXPORT FrechetRand : public ContinuousDistribution
{
    double alpha, s, m;
    double alphaInv; /// 1.0 / alpha

public:
    FrechetRand(double shape, double scale, double location);
    std::string name() override;

    void setParameters(double shape, double scale, double location);
    inline double getShape() const { return alpha; }
    inline double getScale() const { return s; }
    inline double getLocation() const { return m; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    double Entropy() const;
};

#endif // FRECHETRAND_H
