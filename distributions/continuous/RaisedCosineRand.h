#ifndef RAISEDCOSINERAND_H
#define RAISEDCOSINERAND_H

#include "ContinuousRand.h"

/**
 * @brief The RaisedCosineRand class
 */
class RANDLIBSHARED_EXPORT RaisedCosineRand : public ContinuousRand
{
    double mu, s;
    double s_pi; /// M_1_PI * s
public:
    RaisedCosineRand(double location, double scale);
    std::string name() override;

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

    std::complex<double> CF(double t) const override;

    double Median() const override;
    double Mode() const override;
    double Skewness() const override;
    double ExcessKurtosis() const override;
};


/**
 * @brief The RaabGreenRand class
 * X ~ Raised cosine (0.0, M_PI)
 */
class RANDLIBSHARED_EXPORT RaabGreenRand : public RaisedCosineRand
{
public:
    RaabGreenRand() : RaisedCosineRand(0.0, M_PI) {}
    std::string name() override;
};


#endif // RAISEDCOSINERAND_H
