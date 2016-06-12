#ifndef RADEMACHERRAND_H
#define RADEMACHERRAND_H

#include "DiscreteDistribution.h"
#include "../BasicRandGenerator.h"

/**
 * @brief The RademacherRand class
 */
class RANDLIBSHARED_EXPORT RademacherRand : public DiscreteDistribution
{
public:
    RademacherRand();
    std::string name() override;

    double P(int k) const override;
    double F(int k) const override;
    int variate() const override;

    double Mean() const override { return 0; }
    double Variance() const override { return 1; }

    std::complex<double> CF(double t) const override;

    double Median() const override;
    int Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
    
    double Entropy();
};

#endif // RADEMACHERRAND_H
