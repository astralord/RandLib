#ifndef RADEMACHERRAND_H
#define RADEMACHERRAND_H

#include "DiscreteRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The RademacherRand class
 */
class RANDLIBSHARED_EXPORT RademacherRand : public DiscreteRand<int>
{
public:
    RademacherRand();
    virtual void setName() override;

    virtual double P(int k) const override;
    virtual double F(double x) const override;
    virtual double variate() const override;

    double E() const override { return 0; }
    double Var() const override { return 1; }

    static constexpr double Median() { return 0; }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessiveKurtosis() { return -2; }

    static constexpr double Entropy() { return M_LN2; }
};

#endif // RADEMACHERRAND_H
