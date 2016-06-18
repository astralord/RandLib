#ifndef VONMISESRAND_H
#define VONMISESRAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The VonMisesRand class
 */
class RANDLIBSHARED_EXPORT VonMisesRand : public ContinuousDistribution
{
    double mu, k;
    double I0kInv; /// 1.0 / I_0(k)
    double s; /// generator coefficient

public:
    VonMisesRand(double location, double concentration);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return mu - M_PI; }
    double MaxValue() const override { return mu + M_PI; }

    void setLocation(double location);
    void setConcentration(double concentration);
    inline double getLocation() const { return mu; }
    inline double getConcentration() const { return k; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
};

#endif // VONMISESRAND_H
