#ifndef MAXWELLBOLTZMANNRAND_H
#define MAXWELLBOLTZMANNRAND_H

#include "ChiSquaredRand.h"

/**
 * @brief The MaxwellBoltzmannRand class
 */
class RANDLIBSHARED_EXPORT MaxwellBoltzmannRand : public ContinuousDistribution
{
    double a;

    ChiSquaredRand C;

public:
    explicit MaxwellBoltzmannRand(double scale);
    std::string name() override;

    void setScale(double scale);
    inline double getScale() const { return a; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // MAXWELLBOLTZMANNRAND_H
