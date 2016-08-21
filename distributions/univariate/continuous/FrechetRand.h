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

    std::string Name() const override;
    SUPPORT_TYPE SupportType() const override { return RIGHTSEMIFINITE_T; }
    double MinValue() const override { return m; }
    double MaxValue() const override { return INFINITY; }

    void SetParameters(double shape, double scale, double location);
    inline double GetShape() const { return alpha; }
    inline double GetScale() const { return s; }
    inline double GetLocation() const { return m; }

    double f(double x) const override;
    double F(double x) const override;
    double Variate() const override;

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

#endif // FRECHETRAND_H
