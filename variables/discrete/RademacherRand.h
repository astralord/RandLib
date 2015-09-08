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
    virtual std::string name() override;

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double E() const override { return 0; }
    double Var() const override { return 1; }

    static constexpr double Median() { return 0; }
    static constexpr double Skewness() { return 0; }
    static constexpr double ExcessiveKurtosis() { return -2; }

    std::complex<double> CF(double t) const override;
    
    static constexpr double Entropy() { return M_LN2; }
};

#endif // RADEMACHERRAND_H
