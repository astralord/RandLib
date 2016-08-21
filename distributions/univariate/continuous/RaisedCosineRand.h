#ifndef RAISEDCOSINERAND_H
#define RAISEDCOSINERAND_H

#include "ContinuousDistribution.h"

/**
 * @brief The RaisedCosineRand class
 */
class RANDLIBSHARED_EXPORT RaisedCosineRand : public ContinuousDistribution
{
    double mu, s;
    double s_pi; /// M_1_PI * s
public:
    RaisedCosineRand(double location, double scale);

    std::string name() const override;
    SUPPORT_TYPE supportType() const override { return FINITE_T; }
    double MinValue() const override { return mu - s; }
    double MaxValue() const override { return mu + s; }

    void setLocation(double location);
    inline double getLocation() const { return mu; }
    void setScale(double scale);
    inline double getScale() const { return s; }

    double f(double x) const override;
    double F(double x) const override;

    static double standardVariate();
    double variate() const override;

    double Mean() const override;
    double Variance() const override;
    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;

    std::complex<double> CF(double t) const override;
};


/**
 * @brief The RaabGreenRand class
 * X ~ Raised cosine (0.0, M_PI)
 */
class RANDLIBSHARED_EXPORT RaabGreenRand : public RaisedCosineRand
{
public:
    RaabGreenRand() : RaisedCosineRand(0.0, M_PI) {}
    std::string name() const override;
};


#endif // RAISEDCOSINERAND_H
