#ifndef RADEMACHERRAND_H
#define RADEMACHERRAND_H

#include "DiscreteRand.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The RademacherRand class
 */
class RANDLIBSHARED_EXPORT RademacherRand : public DiscreteRand
{
public:
    RademacherRand();
    virtual std::string name() override;

    double P(int k) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override { return 0; }
    double Variance() const override { return 1; }

    std::complex<double> CF(double t) const override;

    double Quantile(double p) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    static constexpr double Entropy() { return M_LN2; }
};

#endif // RADEMACHERRAND_H
