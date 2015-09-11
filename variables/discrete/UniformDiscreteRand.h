#ifndef UNIFORMDISCRETERAND_H
#define UNIFORMDISCRETERAND_H

#include "DiscreteRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The UniformDiscreteRand class
 */
class RANDLIBSHARED_EXPORT UniformDiscreteRand : public DiscreteRand<int>
{
    int n, a, b;

public:
    UniformDiscreteRand(int minValue = 0, int maxValue = 1);
    virtual std::string name() override;

    void setBoundaries(int minValue, int maxValue);
    inline int getMinValue() const { return a; }
    inline int getMaxValue() const { return b; }

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override;
    double Var() const override;

    double quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    inline double Entropy() { return std::log(n); }
};

#endif // UNIFORMDISCRETERAND_H
