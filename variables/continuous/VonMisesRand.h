#ifndef VONMISESRAND_H
#define VONMISESRAND_H

#include "ContinuousRand.h"

/**
 * @brief The VonMisesRand class
 */
class RANDLIBSHARED_EXPORT VonMisesRand : public ContinuousRand
{
    double mu, k;
    double I0kInv; /// 1.0 / I_0(k)
public:
    VonMisesRand(double location, double concentration);
    virtual std::string name() override;

    void setLocation(double location);
    void setConcentration(double concentration);
    inline double getLocation() const { return mu; }
    inline double getConcentration() const { return k; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    double Median() const override;
    double Mode() const override;
};

#endif // VONMISESRAND_H
