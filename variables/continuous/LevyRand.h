#ifndef LEVYRAND_H
#define LEVYRAND_H

#include "ContinuousRand.h"
#include "NormalRand.h"

/**
 * @brief The LevyRand class
 * X ~ Levy(mu, c)
 * X ~ Stable(0.5, 1, mu, c)
 */
class RANDLIBSHARED_EXPORT LevyRand : public ContinuousRand
{
    double mu, c;
    double sqrtc_2pi; /// sqrtc_2pi = sqrt(0.5 * c / pi)
    NormalRand X;

public:
    LevyRand(double location = 0, double scale = 1);
    std::string name() override;

    void setLocation(double location);
    void setScale(double scale);
    inline double getLocation() const { return mu; }
    inline double getScale() const { return c; }

    double f(double x) const override;
    double F(double x) const override;
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
    
    std::complex<double> CF(double t) const override;
    double Quantile(double p) const override;
    
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};

#endif // LEVYRAND_H
