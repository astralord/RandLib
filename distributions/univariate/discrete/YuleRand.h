#ifndef YULERAND_H
#define YULERAND_H

#include "DiscreteDistribution.h"
#include "GeometricRand.h"
#include "../continuous/ExponentialRand.h"
#include "../continuous/ParetoRand.h"

/**
 * @brief The YuleRand class
 * Yule distribution
 * X ~ Yule(ρ)
 */
class RANDLIBSHARED_EXPORT YuleRand : public DiscreteDistribution
{
    double ro;
    double lgamma1pRo;
    
    ParetoRand X;
public:
    explicit YuleRand(double shape);
    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return RIGHTSEMIFINITE_T; }
    int MinValue() const override { return 1; }
    int MaxValue() const override { return INT_MAX; }

    void setShape(double shape);
    inline double getShape() const { return ro; }

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;
    static int variate(double shape);

    double Mean() const override;
    double Variance() const override;

    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // YULERAND_H
